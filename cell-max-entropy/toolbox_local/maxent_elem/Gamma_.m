function [gam,dgam,hgam,p_a,Z]=Gamma_(ndim,x_a,x,lam,near,q)

sum2=0;
for id=1:ndim
    sum2=sum2 + lam(id)*(x(id)-x_a(near,id));
end

temp=q'.*exp(sum2);
Z=sum(temp);
p_a=temp/Z;
gam=log(Z);

%CHECK IT!
p_a(q<eps)=q(q<eps);


dgam = zeros(1,ndim);
for id=1:ndim
    dgam(id)=sum( (x(id)-x_a(near,id)).*p_a(:) );
end

hgam = zeros(ndim,ndim);
for id=1:ndim
    for jd=1:ndim
        hgam(id,jd)=sum( p_a(:).* ...
            (x(id)-x_a(near,id)).*(x(jd)-x_a(near,jd)) )  ...
            - dgam(id)*dgam(jd);
    end
end
