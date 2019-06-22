function [NegLnLike,Psi,U]=likelihood(x)
% Calculates the negative of the concentrated ln 每 likelihood
%
%Inputs:
% x 每 vetor of log(theta) parameters
%
%Global variables used:
% ModelInfo.X 每 n x k matrix of sample locations
% ModelInfo.y 每 n x 1 vector of observed data
%
%Outputs:
% NegLnLike 每 concentrated ln 每 likelihood for minimizing
% Psi 每 correlation matrix
% U 每 Cholesky factorization of correlation matrix
global ModelInfo
X=ModelInfo.X;
y=ModelInfo.y;
theta=10.^x;
n=size(X,1);
one=ones(n,1);
% Pre每allocate memory
Psi=zeros(n,n);
% Build upper half of correlation matrix
for i=1:n
    for j=i+1:n
        Psi(i,j)=exp(-sum(theta.*(X(i,:)-X(j,:)).^2));
    end
end
% Add upper and lower halves and diagonal of ones plus
% small number to reduce ill conditioning

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi=Psi+Psi'+eye(n)+eye(n).*eps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cholesky factorization
[U,p]=chol(Psi);
% Use penalty if ill每conditioned
if p > 0
    NegLnLike=1e4;
else
    % Sum lns of diagonal to find ln(det(Psi))
    LnDetPsi=2*sum(log(abs(diag(U))));
    % Use back每substitution of Cholesky instead of inverse
    mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
    SigmaSqr=((y-one*mu)'*(U\(U'\ (y-one*mu))))/n;
    NegLnLike=(n/2)*log(SigmaSqr)+0.5*LnDetPsi;
end
end