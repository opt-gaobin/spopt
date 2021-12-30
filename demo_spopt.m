function demo_spopt
%%-------------------------------------------------------------------------
% This demo shows how to call spopt to solve
%       min  f(X), s.t.  X'*J2n*X=J2k.
% where J2k = [zeros(k,k),eye(k);-eye(k),zeros(k,k)].
% -------------------------------------
% objective 1:    nearest symplectic matrix problem, f(X):= norm(X-A,'fro')^2
% objective 2:    trace minimization, f(X):= trace(X'*A*X)
% solver:       Cayley, quasi-geodesic, SR decomposition retraction (Canonical-like metric, Euclidean metric)
% output:       function infomation, iterative figures
% -------------------------------------
% Author: Bin Gao (https://www.gaobin.cc)
%   Version 1.0 ... 2020/06
%   Version 1.1 ... 2021/03: add objective function 2 and comparison with the SR decomposition-based retraction by Nguyen Thanh Son (https://sites.google.com/view/ntson)
%--------------------------------------------------------------------------
%% objective function 1
    function [F,G] = fun(X,A)
        F = norm(X-A,'fro')^2;
        G = 2*(X-A); % gradient
    end

% --- Problem generation ---
n = 200; k = 40; 

% --- scenario 1:
% A is randomly generated as a perturbation (in the normal space) of 
% a symplectic matrix, which means that the problem has the closed-form 
% solution, and achieves the optimal function value f_star = 0.
% J2n = [zeros(n) eye(n);-eye(n) zeros(n)]; J2k = [zeros(k) eye(k);-eye(k) zeros(k)];
% WA = randn(2*n,2*n); WA = WA'*WA+0.1*eye(2*n); EA = expm([WA(n+1:end,:); -WA(1:n,:)]);
% A = [EA(:,1:k) EA(:,n+1:n+k)]; 
% s = 1e-8; K = rand(2*k,2*k); K = K - K'; B = A; A = B + s*J2n*(B*K);

% --- scenario 2:
% A is totally randomly generated
A = randn(2*n,2*k); A = A/norm(A,2); 
%% --------------------------------------------------------------
% objective function 2
%     function [F,G] = fun(X,A)
%         AX = A*X;
%         F = X(:)'*AX(:);
%         G = 2*AX;
%     end
% 
% % --- Problem generation ---
% n = 200; k = 40;
% 
% Gauss.c = 1.2;
% Gauss.m = round(n/5);
% Gauss.d = -sqrt(Gauss.m);
% DD = diag([1:n 1:n]);
% rng default
% Q1 = randn(n,n); Q2 = randn(n,n);
% [U,~,~] = svd(Q1+sqrt(-1)*Q2);
% A = [real(U) -imag(U);imag(U) real(U)];
% d1 = [ones(1,Gauss.m-2) Gauss.c Gauss.c ones(1,n-Gauss.m)];
% d2 = [ones(1,Gauss.m-2) 1/Gauss.c 1/Gauss.c ones(1,n-Gauss.m)];
% L2 = zeros(n,n); L2(Gauss.m,Gauss.m-1) = Gauss.d;
% L2(Gauss.m-1,Gauss.m) = Gauss.d;
% Lmcd = [diag(d1) L2;zeros(n,n) diag(d2)]; %sympl. Gauss transformation type I
% A = A*Lmcd;
% A = 0.5*((A*DD)*A' + A*(DD*A'));
% A = 0.5*(A+A');
%% --------------------------------------------------------------
% --- parameters ---
opts.record = 1;
opts.mxitr  = 1000;
opts.gtol = 1e-6;

% --- generate initial guess ---
% type 1: "identity"
% X0 = zeros(2*n,2*k); X0(1:k,1:k) = eye(k); X0(n+1:n+k,k+1:end) = eye(k);
% type 2: random
W = randn(2*k,2*k); W = W'*W+0.1*eye(2*k); E = expm([W(k+1:end,:); -W(1:k,:)]);
X0 = [E(1:k,:);zeros(n-k,2*k);E(k+1:end,:);zeros(n-k,2*k)];
%% ------------------------------------------------------------------------
% call solver
opts.metric = 1;

% --- Cayley retraction ---
opts.retr = 1;
tic; [~, out1]= spopt(X0, @fun, opts, A); tsolve1 = toc;
% --- Quasi-geodesic ---
opts.retr = 2;
tic; [~, out2]= spopt(X0, @fun, opts, A); tsolve2 = toc;
% --- SR decomposition based ---
opts.retr = 3;
tic; [~, out3]= spopt(X0, @fun, opts, A); tsolve3 = toc;

fprintf('spopt-cay: obj: %7.6e, itr: %d, nrmG: %3.2e, nfe: %d, time: %f, |XT*JX-J|: %3.2e \n', ...
    out1.fval, out1.itr, out1.nrmG, out1.nfe, tsolve1, out1.feaX);
fprintf('spopt-geo: obj: %7.6e, itr: %d, nrmG: %3.2e, nfe: %d, time: %f, |XT*JX-J|: %3.2e \n', ...
    out2.fval, out2.itr, out2.nrmG, out2.nfe, tsolve2, out2.feaX);
fprintf('spopt-SR : obj: %7.6e, itr: %d, nrmG: %3.2e, nfe: %d, time: %f, |XT*JX-J|: %3.2e \n', ...
    out3.fval, out3.itr, out3.nrmG, out3.nfe, tsolve3, out3.feaX);
%% ------------------------------------------------------------------------
% figure
% --- function value ---
f_fval = figure(1); hold off;
semilogy(out1.times,out1.fvals,'r-', ...
    out2.times,out2.fvals,'b--', ...
    out3.times,out3.fvals,'k:', ...
    'linewidth',1.5);
set(gca,'fontsize',20);  % grid on
xlabel('time (s)','fontsize',20);
ylabel('function value','fontsize',20);
title(['Size: ',num2str(2*n),'\times',num2str(2*k)],'fontsize',20);
legend('Sp-Cayley','Quasi-geodesics','SR','Location','northeast');
% --- gradient ---
f_kkt = figure(2); hold off;
semilogy(out1.times, out1.kkts,'r-o', ...
     out2.times, out2.kkts,'b--+', ... 
     out3.times, out3.kkts,'k:>', ... 
     'linewidth',1.5);
set(gca,'fontsize',20);  % grid on
xlabel('time (s)','fontsize',20);
ylabel('||grad f||','fontsize',20);
title(['Size: ',num2str(2*n),'\times',num2str(2*k)],'fontsize',20);
legend('Sp-Cayley','Quasi-geodesics','SR','Location','northeast');
end

