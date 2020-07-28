function [gam]=Gamma(ndim,x_a,x,near,q,lam0,dlam,t)

lam = lam0 + t*dlam';
sum2=0;
for id=1:ndim
    sum2=sum2 + lam(id)*(x(id)-x_a(near,id));
end

temp=q'.*exp(sum2);
Z=sum(temp);
%p_a=temp/Z;
gam=log(Z);
