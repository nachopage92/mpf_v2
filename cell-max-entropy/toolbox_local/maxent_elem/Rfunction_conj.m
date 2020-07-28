function [D, dD] = Rfunction_conj(mesh,v,p,apow,bpow)
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
m    = size(p,1);
rho  = zeros(N-1,1); % N-1 edges
drho = zeros(N-1,2); % N-1 edges

D    = zeros(m,1);
dD   = zeros(m,2);

for j = 1:m
  x = p(j,:);
  
  % Compute the nonnegative edge weight and its derivative
  [rho,drho] = rhoweight(mesh,v,x,rho,drho);
  
  % weight function is the square of the R-function
  [w,dw] = formR(rho,drho,bpow);
  
  D(j) =  w^apow;
  dD(j,:) = apow.*w^(apow-1)*dw;
  
  if (sum(isnan(dD(j,:))) > 0) % NaN values at the vertices of the polygon
    dD(j,:) = [0 0];
  end
  
end

return
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [rho,drho] = rhoweight(mesh,v,x,rho,drho)
% % Purpose
% % =======
% % Compute the nonnegative edge weight and its derivative for the i-th edge
% % of the polygon using the triangle inequality
% %
% N    = size(mesh,1);  % N edges
% for i = 1:N
%   v1 = v(mesh(i,1),:);
%   v2 = v(mesh(i,2),:);
%   s1   = v1 - x; 
%   s2   = v2 - x; 
%   s    = v2 - v1;
%   
%   rho(i)    = norm(s1) + norm(s2) - norm(s);
%   drho(i,:) = - s1/norm(s1) - s2/norm(s2);
% end
% return
% end

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
  if abs(f)<eps
    rho (i)   = 0;
    drho(i,:) = [0 0];
  else
    df = [s(2), -s(1)]*id;
    t  = ( (0.5*d)^2 - sum((x-xc).^2) )*id;
    dt = -2*(x-xc)*id;

    val = sqrt(t^2+f^4);

    rho (i)   = sqrt(f^2 + 0.25*(val-t)^2 );
    drho(i,:) = (f*df + 0.25*(val-t) * ((t*dt+2*f^3*df)/val - dt) )/rho(i);
  end
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formR(rho,drho,bpow)
% Purpose
% =======
% Compute the R function (conjuntion) and its gradient. 
% 
% NOTE: This operation is NOT associative which means that the constructed field depends
%       on the order in which the individual segments are joined.

N = size(rho,1);
w  = rho(1);
dw = drho(1,:);
for i = 2:N
%   wp = sqrt(w^2 + rho(i)^2);
%   w  = w + rho(i) - wp;
%   dw = dw + drho(i,:) - ( w*dw + rho(i)*drho(i,:) ) / wp;

  wp = (w^bpow + rho(i)^bpow);
  w  = w + rho(i) - wp^(1/bpow);
  dw = dw + drho(i,:) - wp^(1/bpow-1) * ( w^(bpow-1)*dw + rho(i)^(bpow-1)*drho(i,:) );
end

return
end
