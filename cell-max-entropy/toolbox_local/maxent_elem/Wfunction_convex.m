function [W, dW] = Wfunction_convex(mesh,x_a,x_n,x_s,spow,mpow,h)
% function [W, dW] = Wfunction_convex(mesh,x_a,x_n,x_s,spow,mpow)
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

W    = zeros(m,1);
dW   = zeros(m,2);

% rhoA  = zeros(N,1); % N edges
% rhoA  = rhoweightA(mesh,x_n,x_a,rhoA);
% grhoA = (1./rhoA)*[1 1];

h_a   = sqrt(2)*h;

%rhoA  = exp(rhoA.^(-gamma));
%coeff = prod(rhoA);

for k = 1:m
   x = x_s(k,:);
  
  % Compute the nonnegative edge weight and its derivative: cloud boundary function.
  [rho,drho] = rhoweight(mesh,x_n,x,rho,drho);
  
  % Compute the Weighting function and its gradient
  %[w_k,dw_k]  = formW(coeff,gamma,rho,drho);

  %---------------
  %rho = rho./rhoA;
  %grho = spow*(rho.^(-spow-1).*exp(-rho.^(-spow)))./rhoA*[1 1];
  %drho = grho.*drho;
  %rho  = exp(-rho.^(-spow));
  
  %[w_k,dw_k]  = formR(rho,drho,mpow);
  
  %W(k)    = w_k;
  %dW(k,:) = dw_k;

  %-----------------
  %-- CASE 1
  %   rho  = rho./rhoA;
  %   drho = grhoA.*drho;
  %   [w_k,dw_k]  = formR(rho,drho,mpow);
  %
  %   W(k)    = w_k^spow;
  %   dW(k,:) = spow*w_k^(spow-1)*dw_k;
  
  %-- CASE 2
  [w_k,dw_k]  = formR(rho,drho,mpow);
  q_a     = exp(-(x-x_a)*(x-x_a)'/(h_a*h_a));
  dq_a    = q_a*(-2*(x-x_a))/(h_a*h_a);
  W(k)    = w_k^spow * q_a;
  dW(k,:) = spow*w_k^(spow-1)*dw_k * q_a + w_k^spow * dq_a;

  %-- CASE 3
  %   W(k)    = exp(-w_k^(-spow));
  %   dW(k,:) = exp(-w_k^(-spow))*spow*w_k^(-spow-1)*dw_k;

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
function [rho,drho] = rhoweight(mesh,x_n,x,rho,drho)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative
% of the polygon using the triangle inequality
%
N    = size(mesh,1);  % N edges
for j = 1:N
  %The edges ar in CCW (Counter Clock Wise) so we define x1 and x2 in CW (clock wise)
  x1 = x_n(mesh(j,1),:);
  x2 = x_n(mesh(j,2),:);
  
  s1 = x - x1;

  s  = x2 - x1;
  d  = norm(s);
  id = 1/d;
  
  xc = (x1+x2)*0.5;
  f  = (s1(1)*s(2) - s1(2)*s(1))*id;
  
  df = [s(2), -s(1)]*id; %in WWC it points in the outer direction
  t  = ( (0.5*d)^2 - sum((x-xc).^2) )*id;
  dt = -2*(x-xc)*id;
  
  val = sqrt(t*t+f^4);
  
  h_j = sqrt(f^2 + 0.25*(val-t)^2 );
  
  if h_j<2*eps  %Spacing of floating point numbers: 2^(-52) ~ 2.22e-16
    rho(j)    = 0;
    drho(j,:) = df;
  else
    rho(j)    = h_j;
    drho(j,:) = (f*df + 0.25*(val-t) * ((t*dt+2*f^3*df)/val - dt) )/h_j;
  end
end


return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho] = rhoweightA(mesh,x_n,x,rho)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative
% of the polygon using the triangle inequality
%
N    = size(mesh,1);  % N edges
for j = 1:N
  %The edges ar in CCW (Counter Clock Wise) so we define x1 and x2 in CW (clock wise)
  x1 = x_n(mesh(j,1),:);
  x2 = x_n(mesh(j,2),:);
  
%   s1 = x - x1;
% 
%   s  = x2 - x1;
%   d  = norm(s);
%   id = 1/d;
%   
%   f  = (s1(1)*s(2) - s1(2)*s(1))*id;
  
  xc = (x1+x2)*0.5;
  %t  = ( (0.5*d)^2 - sum((x-xc).^2) )*id;
  
  %val = sqrt(t*t+f^4);
  
  %h_j = sqrt(f^2 + 0.25*(val-t)^2 );

  %rho(j) = max([norm(x-x1) norm(x-x2)]);
  rho(j) = norm(x-xc);
end


return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formW(coeff0,gamma,rho,drho)
% Purpose
% =======
% Compute the R function and its gradient
%

erho = exp(-rho.^(-gamma));
w    = coeff0*prod(erho);

grho = rho.^(-gamma-1)*[1 1];
dw   = w*gamma*sum(grho.*drho);
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formR(rho,drho,mpow)
% Purpose
% =======
% Compute the R function and its gradient
%

N    = length(rho);
ids  = 1:N;
ids  = ids(rho<eps);
nrho0= length(ids);

if nrho0>0
  w = 0;
  if nrho0 == 1
    dw = drho(ids,:);
  else
    dw = mean(drho(ids,:));
  end
else
  ws  = sum(1./(rho.^mpow));
  w   = ws^(-1/mpow);
  dw  = ( (rho.^(-mpow-1))' * drho ) * w/ws;
end

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