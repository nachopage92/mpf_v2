function [gam] = Gamma2_ls(t,exp_beta,dx,lam,dlam)
%function [gam] = Gamma2_ls(t,exp_beta,dx,lam,dlam)
% t        : scalar (line search parameter)
% exp_beta : part of the shape functions related with beta
% dx       : nxD array with the difference vector between x and x_a
% lam      : 1xD array
% dlam     : 1xD array 

lam = lam + t * dlam;
sum2= lam*dx';
temp= exp_beta .* exp(sum2');
Z     = sum(temp);
gam   = log(Z);
