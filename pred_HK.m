function [pred_HK, RMSE] = pred_HK( x )

% calculate the predicted value by Hierarchial Kriging 

%input:
%     site to be predicted
%output:
%     prediction

global ModelInfo_HK
X=ModelInfo_HK.Xe;
U=ModelInfo_HK.U;
F=ModelInfo_HK.F;
V_HK=ModelInfo_HK.V_HK;
theta=ModelInfo_HK.theta;
belta0=ModelInfo_HK.belta0;
sigmaSqr=ModelInfo_HK.sigmaSqr;
n=ModelInfo_HK.ne;
R=U'*U;

r=ones(n,1,'double');
for i= 1:n
    r(i)=cor_function(theta,x,X(i,:));
end
pred_lf=pred(x);
pred_HK=belta0*pred_lf+r'*V_HK;
RMSE=sqrt(abs(sigmaSqr*(1-r'*R^(-1)*r+(r'*R^(-1)*F)-pred_lf)^2*...
    (F'*R^(-1)*F)^(-1)));

end