function training_HK(X,Y)

%tunning the hypermeter of HK model

%input:
%     high fidelity data
global ModelInfo_HK
ModelInfo_HK.Xe=X;
ModelInfo_HK.Ye=Y;
k=size(X,2);
ModelInfo_HK.k=k;
n=size(X,1);
ModelInfo_HK.ne=n;

theta_Lower=-3*ones(1,k);
theta_Upper=2*ones(1,k);
theta=ga(@likelihood_HK,k,[],[],[],[],theta_Lower,theta_Upper);
[~,U,belta0, sigmaSqr, F]= likelihood_HK(theta);
ModelInfo_HK.theta=theta;
ModelInfo_HK.V_HK=(U\(U'\(Y-belta0*F)));
ModelInfo_HK.U=U;
ModelInfo_HK.belta0=belta0;
ModelInfo_HK.sigmaSqr=sigmaSqr;
ModelInfo_HK.F=F;

end