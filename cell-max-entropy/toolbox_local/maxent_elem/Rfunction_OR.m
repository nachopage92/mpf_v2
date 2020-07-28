function [W, dW] = Rfunction_OR(mesh,x_a,x_n,x_s,gamma,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% =======
% Evaluate level set function using R-functions (conjunctions): rho^1 ^
% rho^2 . . . . ^ rho^n-1 where rho^i is the edge-weight for the i-th edge
% USAGE: [D, dD] = Rfunction(v,p,in)
%
% Inputs:
% v  : [x1 y1; x2 y2; ...; xn yn], n vertices of the polygon in ccw; polygon
%      is open (n vertices) or closed (n+1 vertices; repeat 1st vertex)
% p  : [ px(1) py(1); ... ; px(m) py(m) ], the m points for which we want
%      to compute the levet set function and its derivatives
%
% Outputs:
% D  : level set function = [D_1; ...; D_m] (m x 1 array)
% dD : gradient of level set function = [dD_1; ...; dD_m] (m x 2 array)
%
% N. Sukumar: UC Davis, February 2014
% D. Millan UPC-BarcelonaTech, February 2014 

N    = size(mesh,1);  % N edges
m    = size(x_s,1);
rho  = zeros(N-1,1); % N-1 edges
drho = zeros(N-1,2); % N-1 edges

W    = zeros(m,1);
dW   = zeros(m,2);

mpow = 2;

coeff0 = ((1-2^gamma)/log(beta));

for k = 1:m
  x = x_s(k,:);
  
  % Compute the nonnegative edge weight and its derivative
  [rho,drho] = rhoweight(mesh,x_n,x,rho,drho);
  
  %[rho,drho] = cloudweight(mesh,x_a,x_n,x,gamma,coeff0,rho,drho);
  
  % weight function is the square of the R-function
  [w_k,dw_k] = formR(rho,drho,mpow);
  
  W(k)    =  w_k;
  dW(k,:) = dw_k;
  
  if (sum(isnan(dW(k,:))) > 0) % NaN values at the vertices of the polygon
    dW(k,:) = [0 0];
  end
  
end

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,drho] = rhoweight(mesh,v,x,rho,drho)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative of the polygon
%
N    = size(mesh,1);  % N edges
for i = 1:N
  v1 = v(mesh(i,1),:);
  v2 = v(mesh(i,2),:);
  
  s1 = x - v1;
  %s2 = x - v2;
  s  = v2 - v1;
  d  = norm(s);
  id = 1/d;
  
  xc = (v1+v2)*0.5;
  f  = (s1(1)*s(2) - s1(2)*s(1))*id;
  
  df = [s(2), -s(1)]*id;
  t  = ( (0.5*d)^2 - sum((x-xc).^2) )*id;
  dt = -2*(x-xc)*id;
  
  val = sqrt(t*t+f^4);
  
  rho (i)   = sqrt(f^2 + 0.25*(val-t)^2 );
  
  if abs(rho(i))<2*eps  %Spacing of floating point numbers: 2^(-52) ~ 2.22e-16
    drho(i,:) = -df;
  else
    drho(i,:) = (f*df + 0.25*(val-t) * ((t*dt+2*f^3*df)/val - dt) )/rho(i);
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
    drho(j,:) = gamma*xi_j^(-gamma-1)*rho(j)*n_aj;
  else
    rho(j)    = 0;
    drho(j,:) = [0 0];
  end
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formR(rho,drho,mpow)
% Purpose
% =======
% Compute the R function (conjuntion) and its gradient. 
% 
% NOTE: This operation is NOT associative which means that the constructed field depends
%       on the order in which the individual segments are joined.

N = size(rho,1);
w  = rho(1);
dw = drho(1,:);
for j = 2:N
  f  = rho(j);
  wp = w^2 + f^2;
  w  = (w + rho(j) + sqrt(wp))*wp^(mpow/2);

  
  df = drho(j,:);
  dw = (dw + df + (w*dw+f*df)/sqrt(wp))*(wp)^(mpow/2) + w*mpow*(w*dw+f*df)/wp;
end

return
end
