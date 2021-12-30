function [X, out]= spopt(X, fun, opts, varargin)
%%-------------------------------------------------------------------------
% spopt is a solver for optimization on the symplectic Stiefel manifold:
%
%   min f(X), s.t., X'*J_{2n}*X = J_{2k},
%
%   where X \in R^{2n,2k} and J_{2k} = [zeros(k,k),eye(k);-eye(k),zeros(k,k)]
%
%   ------- main iterative update: X(t) -------
%   P = I - J*X*inv(X'*X)*X'*J' + 0.5*rho*X*X'; (or) P = (I-X*J*X'*J')(I-X*J*X'*J')' + 0.5*rho*X*X';
%   U = [-P*G, X*J];    V = [X*J, -P*G];
%   X(t) = X + t*U*inv(I+0.5*t*V'*J'*U)*V'*J*X
%   ---------------------------------------------
%
%   Input:
%           X --- 2n by 2k matrix such that X'*J_{2n}*X = J_{2k}
%
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                   [X, out]= spopt(X0, @fun, opts, data1, data2);
%
%        opts --- option structure with fields:
%                 record      0: no print out
%                 mxitr       max number of iterations
%                 xtol        stop control for ||X_k - X_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                 metric      Riemannian metric: 1 for canonical-like, else for Euclidean
%                 retr        retraction map: 1 for Cayley, 2 for quasi-geodesic, else for SR decomposition
%
%   Output:
%           X --- solution
%         Out --- output information
%                 feaX (final feasibility)     --  feaXs (iterative history)
%                 nrmG (final gradient norm)   --  kkts  (iterative history)
%                 fval (final function value)  --  fvals (iterative history)
%                 itr  (final iteration number)--  times (time line for iterations)
% -------------------------------------
% Reference:
%   Bin Gao, Nguyen Thanh Son, P.-A. Absil, Tatjana Stykel
%   1. Riemannian optimization on the symplectic Stiefel manifold
%   2. Geometry of the symplectic Stiefel manifold endowed with the Euclidean metric
%   3. Symplectic model order reduction of Hamiltonian systems: a Riemannian optimization approach
% Author: Bin Gao (https://www.gaobin.cc)
%   Version 0.1 ... 2019/11
%   Version 0.2 ... 2020/03: add Eucliean metric
%   Version 1.0 ... 2020/06: release at github: https://github.com/opt-gaobin/spopt
%   Version 1.1 ... 2021/02: add Cayley retraction to Euclidean metric
%   Version 1.2 ... 2021/03: add SR decomposition-based retraction by Nguyen Thanh Son (https://sites.google.com/view/ntson)
%--------------------------------------------------------------------------
%% default setting
if nargin < 2; error('[X, out]= spopt(X0, @fun, opts)'); end
if nargin < 3; opts = []; end

% problem size
if isempty(X)
    error('input X is an empty matrix');
else
    [nn, kk] = size(X);
    n = round(nn/2); k = round(kk/2);
end

% stop and print
if ~isfield(opts, 'xtol');      opts.xtol = 1e-6;      end
if ~isfield(opts, 'gtol');      opts.gtol = 1e-6;      end
if ~isfield(opts, 'ftol');      opts.ftol = 1e-12;     end
if ~isfield(opts, 'mxitr');     opts.mxitr  = 1000;    end
if ~isfield(opts, 'record');    opts.record = 0;       end

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'stepsizeRule');  opts.stepsizeRule = 1;   end % different strategies for stepsize
if ~isfield(opts, 'stepsize');  opts.stepsize = 1e-3;  end % initial stepsize (t)
if ~isfield(opts, 'maxgamma');  opts.maxgamma = 1e5;   end % maximal stepsize
if ~isfield(opts, 'mingamma');  opts.mingamma = 1e-15; end % nimimal stepsize
if ~isfield(opts, 'beta');      opts.beta     = 1e-4;  end % linear-search condition
if ~isfield(opts, 'delta');     opts.delta    = 0.1;   end % backtracking parameter
if ~isfield(opts, 'alpha');     opts.alpha    = 0.85;  end % non-monotone parameter
if ~isfield(opts, 'nt');        opts.nt       = 5;     end % max number of linear-search steps

% parameter for Riemannian optimization
if ~isfield(opts, 'retr');      opts.retr   = 1;       end % 1: Cayley, 2: quasi-geodesic, else: SR decomposition
if ~isfield(opts, 'pg');        opts.pg     = 1;       end % different choices of Riemannian gradient (different complement of X)
if ~isfield(opts, 'metric');    opts.metric = 1;       end % 1: Canonical-like, o.w.: Euclidean
if ~isfield(opts, 'rho')                                   % parameter of canonical-like metric (rho)
    if opts.pg == 1
        opts.rho = 0.5;
    else
        opts.rho = 1;
    end
end

% copy parameters
xtol     = opts.xtol;     gtol     = opts.gtol;      ftol     = opts.ftol;      nt       = opts.nt;  
beta     = opts.beta;     delta    = opts.delta;     alpha    = opts.alpha;     crit     = ones(nt, 3);
stepsize = opts.stepsize; gamma    = stepsize;       maxgamma = opts.maxgamma;  mingamma = opts.mingamma;       
retr     = opts.retr;     metric   = opts.metric;    rho      = opts.rho;       pg       = opts.pg;
record   = opts.record;  

% save metric and solver
if metric ~= 1; metricname = 'Euclidean'; else; metricname = 'Canonical-like'; end
if retr == 1
    retrname   = 'Cayley';
elseif retr == 2
    retrname   = 'quasi-geodesic';
else
    retrname = 'SR decomposition';
end
%% ------------------------------------------------------------------------
% Initialization
J2k = [zeros(k) eye(k);-eye(k) zeros(k)];
% evaluate function and gradient info.
[F,  G] = feval(fun, X , varargin{:});  out.nfe = 1;
% preparations for the first update
XJ = [-X(:,k+1:end) X(:,1:k)]; JX = [X(n+1:end,:); -X(1:n,:)];
if metric == 1
    GX = 0.5*rho*G'*X;
    if k < n
        if pg == 1
            XX = X'*X; [invXXXJG, ~] = linsolve(XX, JX'*G);
            PG = G - JX*invXXXJG + X*GX';
        else
            XJG = XJ'*G; JXXJG = JX*XJG;
            PG = G - JXXJG- XJ*(JX'*G)+ XJ*(JX'*JXXJG) + X*GX';
        end
    else
        PG = X*GX';
    end
else
    XX = X'*X; JXG = JX'*G; skewXJG = JXG-JXG';
end

if retr == 1
    invH = true; eye2n = eye(2*n); if k < n/2; invH = false; eye2k = eye(2*k); eye4k = eye(4*k); end
    if metric == 1
        if invH
            PGXJ = -PG*(XJ)'; H = PGXJ + PGXJ'; HJ = [-H(:,n+1:end), H(:,1:n)]; RJX = H*JX;
        else
            % direct way to get VJU and VJX (eq.(5.7) in the SIOPT paper)
            U =  [-PG, XJ]; V = [XJ, -PG];	VJU = V'*[-U(n+1:end,:); U(1:n,:)];
            VJX = V'*JX;
            
            % economic but may lose feasibility (eq.(5.8) in the SIOPT paper)
            % U =  [-PG, XJ]; PGJPG = [PG(n+1:end,:); -PG(1:n,:)]'*PG; VJU = [GX' J2k';PGJPG -GX];
            % VJX = [eye2k; GX(:,k+1:end) -GX(:,1:k)];
        end
    else
        % activate only if "Euclidean + Cayley" ----------------
        Omega = lyap(XX,-skewXJG); Omega = 0.5*(Omega - Omega'); W = -JXG + XX*Omega; W = -0.5*(W+W');
        PG = G - XJ*(0.5*JXG) - JX*Omega + XJ*(0.5*(XX*Omega));
        if invH
            PGXJ = -PG*(XJ)'; H = PGXJ + PGXJ'; HJ = [-H(:,n+1:end), H(:,1:n)]; RJX = H*JX;
        else
            U =  [-PG, XJ]; V = [XJ, -PG];	VJU = V'*[-U(n+1:end,:); U(1:n,:)];
            VJX = V'*JX;
        end
    end
elseif retr == 2
    eye2k = eye(2*k);
    if metric == 1
        W = [-GX(:,k+1:end) GX(:,1:k)]; W = W+W';
        % direct way to get W
        % XJG = XJ'*G; W = 0.5*rho*(XJG+XJG');
    else
        Omega = lyap(XX,-skewXJG); Omega = 0.5*(Omega - Omega'); W = -JXG + XX*Omega; W = -0.5*(W+W');
    end
else
    if metric ~= 1
        Omega = lyap(XX,-skewXJG); Omega = 0.5*(Omega - Omega');
    end
end
% compute initial error
XFeasi = norm(X'*JX - J2k,'fro');
if metric == 1; dtX = PG*(XJ'*JX) + XJ*(PG'*JX); else; dtX = G - JX*Omega; end
nrmG = norm(dtX, 'fro'); nrmG0 = nrmG; % initial gradient norm
% -------------------------------------------------------
% save history
out.fvals = []; out.fvals(1) = F;
out.kkts = [];  out.kkts(1) = nrmG;
out.feaXs = []; out.feaXs(1) = XFeasi;
out.times = []; out.times(1) = 0; runTime = tic;
% line-search parameter
Q = 1; Cval = F;
% print info.
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '------------- Riemannian Gradient Method with Line search -------------- \n');
    fprintf(fid, 'Solver setting... (%s metric, %s retraction) \n',metricname,retrname);
    fprintf(fid, '------\n');
    fprintf(fid, '%4s %9s %10s %9s %9s %9s %9s\n', 'Iter', 'stepsize', 'F(X)', 'nrmG', 'XDiff', 'FDiff', 'XFeasi');
    fprintf('%4d  %8s  %4.3e  %3.2e  %10s  %6s  %3.2e  %2s\n', ...
        0, '', F, nrmG, '', '', XFeasi, '');
end
%% ------------------------------------------------------------------------
% main iteration
for itr = 1 : opts.mxitr
    XP = X;     FP = F;   dtXP = dtX; % GP = G;
    
    % -------- scale step size by nonmonotone line-search --------
    stepsize = gamma; % initial trial stepsize
    nls = 1; deriv = beta*abs(iprod(G,dtXP)); % Riemannian metric ||dtXP||^2_X
    while 1
        % ----- retraction step -----
        if retr == 1
            % Cayley
            if invH
                [X, ~] = linsolve(eye2n - (0.5*stepsize)*HJ, XP + (0.5*stepsize)*RJX);
            else
                [aa, ~] = linsolve(eye4k + (0.5*stepsize)*VJU, VJX);
                X = XP + U*(stepsize*aa);
            end
        elseif retr == 2
            % quasi-geodesic
            U = [XP,-stepsize*dtX]; JWt = stepsize*[W(k+1:end,:); -W(1:k,:)];
            if nls == 1; JZJZ = [-dtX(:,k+1:end) dtX(:,1:k)]'*[dtX(n+1:end,:); -dtX(1:n,:)]; end
            H = [-JWt -(stepsize^2)*JZJZ; eye2k -JWt]; ExpH = expm(H);
            X = U*(ExpH(:,1:2*k)*expm(JWt));
            
            % direct computation based on quasi-geodesic curve
            % U = [XP,-dtX]; JW = [W(k+1:end,:); -W(1:k,:)];
            % if nls == 1; JZJZ = [-dtX(:,k+1:end) dtX(:,1:k)]'*[dtX(n+1:end,:); -dtX(1:n,:)]; end;
            % H = [-JW -JZJZ; eye2k -JW]; ExpH = expm(stepsize*H);
            % X = U*(ExpH(:,1:2*k)*expm(stepsize*JW));
        else
            % SR decomposition based
            X = MsGS(XP - stepsize*dtX);
        end
        
        % ----- evaluate function -----
        [F,G] = feval(fun, X, varargin{:});
        out.nfe = out.nfe + 1;
        
        % ----- line-search -----
        if F <= Cval - stepsize*deriv || nls >= nt
            break;
        end
        stepsize = delta*stepsize;          nls = nls+1;
    end
    
    % -------------------- prepare retraction --------------------
    XJ = [-X(:,k+1:end) X(:,1:k)]; JX = [X(n+1:end,:); -X(1:n,:)];
    if metric == 1
        GX = 0.5*rho*G'*X;
        if k < n
            if pg == 1
                XX = X'*X; [invXXXJG, ~] = linsolve(XX, JX'*G);
                PG = G - JX*invXXXJG + X*GX';
            else
                XJG = XJ'*G; JXXJG = JX*XJG;
                PG = G - JXXJG - XJ*(JX'*G) + XJ*(JX'*JXXJG) + X*GX';
            end
        else
            PG = X*GX';
        end
    else
        XX = X'*X; JXG = JX'*G; skewXJG = JXG-JXG';
    end
    
    if retr == 1
        if metric == 1
            if invH
                PGXJ = -PG*(XJ)';  H = PGXJ + PGXJ';  HJ = [-H(:,n+1:end) H(:,1:n)];  RJX = H*JX;
            else
                % direct way to get VJU and VJX (eq.(5.7) in the SIOPT paper)
                U =  [-PG, XJ]; V = [XJ, -PG];	VJU = V'*[-U(n+1:end,:); U(1:n,:)];
                VJX = V'*JX;
                
                % economic but may lose feasibility (eq.(5.8) in the SIOPT paper)
                % U =  [-PG, XJ]; PGJPG = [PG(n+1:end,:); -PG(1:n,:)]'*PG; VJU = [GX' J2k';PGJPG -GX];
                % VJX = [eye2k; GX(:,k+1:end) -GX(:,1:k)];
            end
        else
            % activate only if "Euclidean + Cayley" ----------------
            Omega = lyap(XX,-skewXJG); Omega = 0.5*(Omega - Omega'); W = -JXG + XX*Omega; W = -0.5*(W+W');
            PG = G - XJ*(0.5*JXG) - JX*Omega + XJ*(0.5*(XX*Omega));
            if invH
                PGXJ = -PG*(XJ)'; H = PGXJ + PGXJ'; HJ = [-H(:,n+1:end), H(:,1:n)]; RJX = H*JX;
            else
                U =  [-PG, XJ]; V = [XJ, -PG];	VJU = V'*[-U(n+1:end,:); U(1:n,:)];
                VJX = V'*JX;
            end
        end
    elseif retr == 2
        if metric == 1
            W = [-GX(:,k+1:end) GX(:,1:k)]; W = W+W';
            % direct way to get W
            % XJG = XJ'*G; W = 0.5*rho*(XJG+XJG');
        else
            Omega = lyap(XX,-skewXJG); Omega = 0.5*(Omega - Omega'); W = -JXG + XX*Omega; W = -0.5*(W+W');
        end
    else
        if metric ~= 1
            Omega = lyap(XX,-skewXJG); Omega = 0.5*(Omega - Omega');            
        end
    end
    
    % ---------------------- compute error ----------------------
    S = X - XP; XDiff = norm(S,'fro')/sqrt(n); FDiff = abs(FP-F)/(abs(FP)+1); XFeasi = norm(X'*JX - J2k,'fro');
    stepsizek = stepsize;
    
    if metric == 1; dtX = PG*(XJ'*JX) + XJ*(PG'*JX); else; dtX = G - JX*Omega; end; nrmG  = norm(dtX, 'fro');
    out.fvals(itr+1) = F; out.kkts(itr+1) = nrmG; out.feaXs(itr+1) = XFeasi; out.times(itr+1) = toc(runTime);
    % print history
    if (record >= 1 && mod(itr,15) == 0)
        fprintf(fid, '%4d  %3.2e  %4.3e  %3.2e  %3.2e  %3.2e  %3.2e  %2d\n', ...
            itr, stepsizek, F, nrmG, XDiff, FDiff, XFeasi, nls);
    end
    
    % -------------------- update step size ---------------------
    switch opts.stepsizeRule
        case 1 % Alternating BB
            %Y = G - GP;     SY = abs(iprod(S,Y));
            Y = dtX - dtXP;     SY = abs(iprod(S,Y));
            if mod(itr,2)==0
                gamma = (norm(S,'fro')^2)/SY;
            else
                gamma  = SY/(norm(Y,'fro')^2);
            end
        case 2 % BB1
            Y = dtX - dtXP;     SY = abs(iprod(S,Y));
            gamma = (norm(S,'fro')^2)/SY;
        case 3 % BB2
            Y = dtX - dtXP;     SY = abs(iprod(S,Y));
            gamma  = SY/(norm(Y,'fro')^2);
        case 4 % Nocedal & Wright, p59(3.60)
            gamma = abs(2*(FP - F)/iprod(G,dtX)); % Riemannian metric
            % gamma = abs(2*(FP - F)/nrmG^2); % Euclidean
    end
    gamma = max(min(gamma, maxgamma), mingamma);
    
    % ----------------------- stop criteria ----------------------
    % crit(itr,:) = [nrmG, XDiff, FDiff];
    % mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);
    if nrmG < gtol
        %if nrmG < gtol*nrmG0
        %if (XDiff < xtol && nrmG < gtol ) || FDiff < ftol
        %if (XDiff < xtol || nrmG < gtol ) || FDiff < ftol
        %if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol
        %if ( XDiff < xtol || FDiff < ftol ) || nrmG < gtol
        %if any(mcrit < [gtol, xtol, ftol])
        %if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol || all(mcrit(2:3) < 10*[xtol, ftol])
        out.msg = 'converge';
        break;
    end
    
    % --------------------- nonmonotone update --------------------
    Qp = Q; Q = alpha*Qp + 1; Cval = (alpha*Qp*Cval + F)/Q;
end
%% ------------------------------------------------------------------------
% output
if itr >= opts.mxitr; out.msg = 'exceed max iteration'; end
if record >= 1
    fprintf(fid, '%s at...\n',out.msg);
    fprintf(fid, '%4d  %3.2e  %4.3e  %3.2e  %3.2e  %3.2e  %3.2e  %2d\n', ...
        itr, stepsizek, F, nrmG, XDiff, FDiff, XFeasi, nls);
    fprintf(fid, '------------------------------------------------------------------------\n');
end
out.feaX = XFeasi;
out.nrmG = nrmG;
out.fval = F;
out.itr = itr;
end
%% ------------------------------------------------------------------------
% nest-function: inner product
function a = iprod(x,y)
a = real(sum(sum(x.*y)));
% a = real(sum(sum(conj(x).*y)));
end
