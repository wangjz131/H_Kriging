function [ NegLnlike, U ,belta0 ,sigmaSqr, F] = likelihood_HK( theta )

%likelihood function of Hierarchial Kriging method

%input:  
%     theta (hyperparameter of the coordinates)
%     sites and responses of high fidelity code
%output:
%     likelihood object function value

global ModelInfo_HK;
x=ModelInfo_HK. Xe;
y=ModelInfo_HK. Ye;

n=size(x,1);
R=zeros(n,n,'double');
F=zeros(n,1,'double');
for i=1:n
    for j=i+1:n
        R(i,j)=cor_function(theta,x(i,:),x(j,:));    
    end
    F(i)=pred(x(i,:));
end

R=R+R'+eye(n)+eps*eye(n);%correlation matrix
[U,p]=chol(R);

if p>0
    NegLnlike=10^4;
else
    LnDetR=2*sum(log(abs(diag(U))));
    belta0=(F'*(U\(U'\y)))/(F'*(U\(U'\F)));
    sigmaSqr=(y-belta0*F)'*(U\(U'\(y-belta0*F)))/n;
    NegLnlike=LnDetR+n*log(sigmaSqr);
end

end