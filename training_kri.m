function  training_kri(S,Y)
%���ƿ����ģ�͵ĳ�����theta
global ModelInfo
ModelInfo.X=S;
ModelInfo.y=Y;
k=size(ModelInfo.X,2);
% Set upper and lower bounds for search of log theta
UpperTheta=ones(1,k).*2;
LowerTheta=ones(1,k).*(-3);
% Run GA search of likelihood
[ModelInfo.Theta,~]=ga(@likelihood,k,[],[],[],[], LowerTheta,UpperTheta);
%ga()���������Ŵ��㷨����Ŀ�꺯������Сֵ
% Put Cholesky factorization of Psi, into ModelInfo
[~,ModelInfo.Psi,ModelInfo.U]=likelihood(ModelInfo.Theta);       
end

