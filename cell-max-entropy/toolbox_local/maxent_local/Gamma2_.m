function [R,J,p_a]=Gamma2_(ndim,exp_beta,dx,lam)
%dx  : nxD array with the difference vector between x and x_a
%lam : 1xD array
%p_a : nx1 array
%beta: 1xn array

sum2= lam*dx';
temp= exp_beta .* exp(sum2');

% sum1= sum(dx.^2,2);
% sum2= lam*dx';
% temp= exp(-beta'.*sum1 + sum2');
Z   = sum(temp);
p_a = temp/Z;
% gam = log(Z);

R = dx'*p_a;  %residual vector

J=zeros(ndim);       %jacobian matrix
for id=1:ndim
  for jd=1:ndim
    J(id,jd)=sum( (p_a.*dx(:,id)).*dx(:,jd) );
  end
end
J = J - R*R';
