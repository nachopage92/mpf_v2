function [W, dW] = Wfunction_noconvex(mesh,x_a,x_n,x_s,gamma,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% =======
%
% Inputs:
% x_n : [x1 y1; x2 y2; ...; xn yn], n vertices of the polygon in CCW (Counter Clock Wise)
%       surronging x_a; polygon is open (n vertices) or closed (n+1 vertices; repeat 1st vertex)
% x_s : [ xs(1) ys(1); ... ; xs(m) ys(m) ], the m sample points for which we want
%       to compute the levet set function and its derivatives
%
% Outputs:
% D  : level set function = [D_1; ...; D_m] (m x 1 array)
% dD : gradient of level set function = [dD_1; ...; dD_m] (m x 2 array)
%
%
% D. Millan UPC-BarcelonaTech, February 2014
%
% REFERENCE:
% [1] C.A. Duarte, D.-J. Kim, D.M. Quaresma, Arbitrarily smooth generalized finite element
%     approximations, CMAME, 2006


N    = size(mesh,1);  % N edges
m    = size(x_s,1);
rho  = zeros(N,1); % N edges
drho = zeros(N,2); % N edges

rhoA = zeros(2,1); % N edges

W    = zeros(m,1);
dW   = zeros(m,2);

mpow = 2;

coeff = ((1-2^gamma)/log(beta));
for k = 1:m
  x = x_s(k,:);
  
  % Compute the nonnegative edge weight and its derivative: cloud boundary function.
  [rho,drho] = cloudweight(mesh,x_a,x_n,x,gamma,coeff,rho,drho);
  
  rhoA = cloudweightA(mesh(7:8,:),x_a,x_n,gamma,coeff,rhoA);
  
  % Compute the Weighting function and its gradient
  [w_k,dw_k]  = formW(coeff,rho(1:6),drho(1:6,:));
  
  % weight function is the square of the R-function 
  srho2 = rho(7)^2 + rho(8)^2;
  if srho2 < eps
    w_r  = 0;
    dw_r = [0 0];
  else
    srhoA2 = rhoA(1)^2 + rhoA(2)^2;
    denA   = (rhoA(1) + rhoA(2) + sqrt(srhoA2)) * (srhoA2)^(mpow/2);

    w_r   = (rho(7) + rho(8) + sqrt(srho2))*(srho2)^(mpow/2)/denA;
    
    f1    = rho(7);
    df1   = rho(7)*drho(7,:);
    f2    = rho(8);
    df2   = rho(8)*drho(8,:);
    dw_r  = (df1 + df2 + (f1*df1 + f2*df2)/sqrt(srho2))*(srho2)^(mpow/2)/denA;
    dw_r  = dw_r + w_r * mpow * (f1*df1+f2*df2)/srho2;   
  end
  
  if norm(w_k*dw_r)>100
      figure(30);clf
      plot(x_n(:,1),x_n(:,2),'r.')
      hold on
      plot(x_a(1),x_a(1),'rs','markersize',14,'linewidth',2)
      for i=1:N
        ii = mesh(i,:);
        plot(x_n(ii,1),x_n(ii,2),'ro-','markersize',8,'linewidth',1)
      end
      plot(x(1),x(2),'cx','markersize',14,'linewidth',3)
      hold off
      axis equal
      pause
  end
  
  
  W(k)    = w_k*w_r;
  dW(k,:) = dw_k*w_r+ w_k*dw_r;
  
  if (sum(isnan(dW(k,:))) > 0) % NaN values at the vertices of the polygon
    %     mesh
    %     rho
    %     drho
    %     w
    %     dw
    %     disp('NaN values at the vertices of the polygon')
    dW(k,isnan(dW(k,:))) = 0;
    %     pause
  end
  
end

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,drho] = cloudweight(mesh,x_a,x_n,x,gamma,coeff0,rho,drho)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative
% of the polygon using the cloud boundary function Eq (5) [1]
%
N    = size(mesh,1);  % N edges

coeff = coeff0^(1/gamma);

for j = 1:N
  
  %The edges ar in CCW (Counter Clock Wise) so we define x1 and x2 in CW (clock wise)
  x1 = x_n(mesh(j,2),:);
  x2 = x_n(mesh(j,1),:);
  
  b_aj = (x1+x2)*0.5;
  d_aj = x_a - b_aj;
  
  %normal to edge j
  s    = x2 - x1;
  n_aj = [s(2) -s(1)]/norm(s);
  
  n_aj = coeff/abs(n_aj*d_aj')*n_aj;
  xi_j = n_aj*(x-b_aj)';
  
  % cloud boundary function and its derivative for the edge j
  %drho(j,:) = gamma*xi_j^(-gamma-1)*rho(j)*n_aj;
  if xi_j > 0
    rho(j)    = exp(-1/(xi_j^gamma));
    drho(j,:) = gamma*xi_j^(-gamma-1)*n_aj;
  else
    rho(j)    = 0;
    drho(j,:) = [0 0];
  end
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhoA = cloudweightA(mesh,x_a,x_n,gamma,coeff0,rhoA)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative
% of the polygon using the cloud boundary function Eq (5) [1]
%
N     = size(mesh,1);  % N edges
coeff = coeff0^(1/gamma);

for j = 1:N
  
  %The edges ar in CCW (Counter Clock Wise) so we define x1 and x2 in CW (clock wise)
  x1 = x_n(mesh(j,2),:);
  x2 = x_n(mesh(j,1),:);
  
  b_aj = (x1+x2)*0.5;
  d_aj = x_a - b_aj;
  
  %normal to edge j
  s    = x2 - x1;
  n_aj = [s(2) -s(1)]/norm(s);
  
  n_aj = coeff/abs(n_aj*d_aj')*n_aj;
  
  xi_aj= n_aj*(x_a-b_aj)';
  
  % cloud boundary function and its derivative for the edge j
  %drho(j,:) = gamma*xi_j^(-gamma-1)*rho(j)*n_aj;
  if xi_aj > 0
    rhoA(j)    = exp(-1/(xi_aj^gamma));
  else
    rhoA(j)    = 0;
  end
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formW(coeff0,rho,drho)
% Purpose
% =======
% Compute the R function and its gradient
%
N    = length(rho);

coeff = N/coeff0;

w  = exp(coeff)*prod(rho);
dw = w*sum(drho);
return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotrho0(mesh,x_a,x_n,x_s,x,rho,dw)
% Purpose
% =======
% Checking zero values of rho
%
N    = size(mesh,1);
ids  = 1:N;
ids0 = ids(rho<eps);
nrh0  = length(ids0);
if nrh0>2
  figure(30);clf
  plot(x_n(:,1),x_n(:,2),'r.')
  hold on
  plot(x_a(1),x_a(1),'rs','markersize',14,'linewidth',2)
  for i=1:nrh0
    ii = mesh(ids0(i),:);
    plot(x_n(ii,1),x_n(ii,2),'ro-','markersize',8)
  end
  plot(x(1),x(2),'cx','markersize',14,'linewidth',2)
  hold off
  axis equal
  fprintf(1,'Rfunction_equivJP: rho have %d zero values at x=[%4.2f %4.2f]\n',nrh0,x(1),x(2));
  
  figure(31);clf
  plot(x_n(:,1),x_n(:,2),'ro-','markersize',12,'linewidth',2)
  hold on
  plot(x_s(:,1),x_s(:,2),'bx',x(1),x(2),'bs')
  quiver(x(1),x(2),dw(1),dw(2))
  hold off
  axis equal
  pause
end

return
end