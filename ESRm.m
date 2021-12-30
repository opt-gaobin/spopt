function [V,H] = ESRm(A,type,tol)
% This function runs the modified elementary SR factorization so that is
% produced H = I when A is already symplectic. In general, H does not need
% to be upper triangular.
% INPUT
% A: 2n x 2 non-isotropic matrix, i.e. A1'*J2*A2 =\ 0.
% type: 'd': produces diagonal matrix H; 'f': produces full matrix H with 1
% on the main diagonal
% tol: the numerical tolerance. Referred value ~ 1e-13 t0 1e-10
% OUTPUT
% V: 2n x 2 symplectic matrix such that A = VH, where H is upper triangle
% Reference 
% [1] A. Salam (2005), On theoretical and numerical aspects of symplectic
% Gram-Schmidt-like algorithms, Numer. Algorithm. 39: 437-462.
% [2] A first-order optimization method on
% symplectic Stiefel manifold, in preparation.
%   Author
% N.T. Son, UC Louvain, 2019-10-18
% N.T. Son, UC Louvain, 2019-11-12
n = size(A,1)/2;
%J = [sparse(n,n) speye(n,n);-speye(n,n) sparse(n,n)];
if nargin < 2
    type = 'd'; 
end
if nargin < 3
    tol = 1e-15;
end
%alp = A(:,2)'*Jmul(A(:,1)); 
alp = A(1:n,2)'*A(n+1:end,1) - A(n+1:end,2)'*A(1:n,1);
V = zeros(size(A));
switch type
    case 'f'
        H = eye(2,2);
        if abs(alp) <= tol
            error('input matrix is isotropic')
        elseif abs(alp+1) <= tol
            %disp('A is already symplectic')
            V = A;
        elseif alp+1 >= tol
            H(1,2) = sqrt(alp+1);
            H(2,1) = H(1,2);
            V(:,1) = (H(2,1)*A(:,2)-A(:,1))/(H(1,2)*H(2,1)-1);
            V(:,2) = (H(1,2)*A(:,1)-A(:,2))/(H(1,2)*H(2,1)-1);
        else
            H(1,2) = sqrt(-alp-1);
            H(2,1) = -H(1,2);
            V(:,1) = (H(2,1)*A(:,2)-A(:,1))/(H(1,2)*H(2,1)-1);
            V(:,2) = (H(1,2)*A(:,1)-A(:,2))/(H(1,2)*H(2,1)-1);
        end
    case 'd'
        H = zeros(2,2);
        if abs(alp) <= tol
            error('input matrix is isotropic')
        elseif abs(alp+1) <= tol
            %disp('A is already symplectic')
            V = A; H(1,1) = 1; H(2,2) = 1;
        else
            H(1,1) = sqrt(abs(alp)); H(2,2) = sign(-alp)*H(1,1);
            V(:,1) = A(:,1)/H(1,1); V(:,2) = A(:,2)/H(2,2);
        end
    otherwise
        error('type must be either f or d')
end
end
            
            


