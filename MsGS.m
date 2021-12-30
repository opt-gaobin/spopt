function [S,R] = MsGS(A,type,tol)
% This function runs the modified  symplectic Gram-Schmidt algorithm
% that produces S,R such that A = SR where V is symplectic and R is an
% permuted upper triangular-type matrix
% INPUT
% A: 2n x 2k at least full rank; There is still an oppotunity for an
% incurable break due to the isotropy.
% type: 'd': produces diagonal matrix H; 'f': produces full matrix H with 1
% on the main diagonal. Default value is 'd' (as the report)
% tol: the numerical tolerance. Default value is 1e-15
% OUTPUT
% S: 2n x 2k symplectic matrix such that A = SR, 

% Reference 
% [1] A. Salam (2005), On theoretical and numerical aspects of symplectic
% Gram-Schmidt-like algorithms, Numer. Algorithm. 39: 437-462.
% [2] A first-order optimization method on
% symplectic Stiefel manifold, in preparation.
%   Author
% N.T. Son, UCLouvain, 2019-11-14
if nargin < 2
    type = 'd'; 
end
if nargin < 3
    tol = 1e-15;
end
[n,k] = size(A); k = k/2;n = n/2;
P = zeros(2*k,2*k);
for i = 1:k
    P(2*i-1,i) = 1;P(2*i,k+i) = 1;
end
%A = A*P';
A = mulPT(A);
S = zeros(size(A));
R = zeros(2*k,2*k);
[S(:,1:2),R(1:2,1:2)] = ESRm(A(:,1:2),type,tol);
for j = 2:k
    Wj = A(:,2*j-1:2*j);
    for i = 1:j-1
        SJW = S(:,2*i-1:2*i)'*[Wj(n+1:end,:);-Wj(1:n,:)];
        R(2*i-1:2*i,2*j-1:2*j) = [-SJW(2,:);SJW(1,:)];
        %R(2*i-1:2*i,2*j-1:2*j) = J2'*(S(:,2*i-1:2*i)'*Jmul(Wj));
        Wj = Wj - S(:,2*i-1:2*i)*R(2*i-1:2*i,2*j-1:2*j);
    end
    [S(:,2*j-1:2*j),R(2*j-1:2*j,2*j-1:2*j)] = ESRm(Wj,type,tol);
end
%S = S*P;
S = mulP(S);
%R = P'*R*P;
R = PTmul(mulP(R));
R = sparse(R);
end
function MPT = mulPT(M)
MPT = zeros(size(M));
kM = size(M,2)/2;
for j = 1:kM
    MPT(:,2*j-1) = M(:,j);
    MPT(:,2*j) = M(:,kM+j);
end
end

function MP = mulP(M)
MP = zeros(size(M));
kM = size(M,2)/2;
for j = 1:kM
     MP(:,j) = M(:,2*j-1);
     MP(:,kM+j) = M(:,2*j);
end
end

function PTM = PTmul(M)
PTM = zeros(size(M));
kM = size(M,1)/2;
for j = 1:kM
    PTM(j,:) = M(2*j-1,:);
    PTM(kM+j,:) = M(2*j,:);
end
end




