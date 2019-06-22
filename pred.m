function [f, MSE] = pred( x )
% Calculates a Kriging prediction at x
%
%Inputs:
% x �C 1 x k vetor of design variables
%
%Global variables used:
% ModelInfo.X �C n x k matrix of sample locations
% ModelInfo.y �C n x 1 vector of observed data
% ModelInfo.Theta �C 1 x k vector of log(theta)
% ModelInfo.U �C n x n Cholesky factorization of Psi
%
%Outputs:
% f �C scalar Kriging prediction
global ModelInfo
% Extract variables from data structure
% slower, but makes code easier to follow
X=ModelInfo.X;
y=ModelInfo.y;
theta=10.^ModelInfo.Theta;
U=ModelInfo.U;
% Calculate number of sample points
n=size(X,1);
% Vector of ones
one=ones(n,1);
% Calculate mu
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
% Initialize psi to vector of ones
psi=ones(n,1);
% Fill psi vector
for i=1:n
    psi(i)=exp(-sum(theta.*abs(X(i,:)-x).^2));
end
% Calculate prediction
f=mu+psi'*(U\(U'\(y-one*mu)));
MSE=((y-one*mu)'*(U\(U'\ (y-one*mu))))/n;
end