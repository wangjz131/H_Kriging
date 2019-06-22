function r=cor_function(theta,x1,x2)

%correlation model
%cubic spline function
theta=10.^theta;
k=size(x1,2);
x1=x1-x2;
r=1;
certh_K=theta.*abs(x1);
for i = 1:k
    if certh_K(i)>=1
        r=0;
        break;
    elseif certh_K(i)<1&&certh_K(i)>0.2
        r=r*1.25*(1-certh_K(i))^3;
    elseif certh_K(i)>=0&&certh_K(i)<=0.2
        r=r*(1-15*certh_K(i)^2+30*certh_K(i)^3);
    end
end

end