function [D, dD] = Rfunction(mesh,x0,p,apower,bpower)
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
% in : boolean flag (0 for points that are not inside; 1 for interior points)
%
% Outputs:
% D  : level set function = [D_1; ...; D_m] (m x 1 array)
% dD : gradient of level set function = [dD_1; ...; dD_m] (m x 2 array)
%
% N. Sukumar: UC Davis, February 2014 
%

N    = size(mesh,1);  % N edges
m    = size(p,1);
rho  = zeros(N,1); % N edges
drho = zeros(N,2); % N edges
%rho_p  = zeros(N,1); % N edges
%drho_p = zeros(N,2); % N edges
%rho_m  = zeros(N,1); % N edges
%drho_m = zeros(N,2); % N edges

D    = zeros(m,1);
dD   = zeros(m,2);

for j = 1:m
 

  x = p(j,:);

%  h = 1e-7;
  
  for i = 1:N
      v1 = x0(mesh(i,1),:);
      v2 = x0(mesh(i,2),:);
      [rho(i),drho(i,:)] = rhoweight(v1,v2,x);
%      x_p = x + [0 h];
%      x_m = x - [0 h];
%      [rho_p(i),drho_p(i,:)] = rhoweight(v1,v2,x_p);
%      [rho_m(i),drho_m(i,:)] = rhoweight(v1,v2,x_m);
%      disp('Numerical')
%      disp((rho_p(i)-rho_m(i))/2/h)
%      disp('Analytical')
%      disp(drho(i,2))
%      pause
  end


  [w,dw] = formR(rho,drho,bpower);
%  [w_p,dw_p] = formR(rho_p,drho_p,bpower);
%  [w_m,dw_m] = formR(rho_m,drho_m,bpower);
%  disp('Numerical')
%  disp((w_p-w_m)/2/h)
%  disp('Analytical')
%  disp(dw(1))
%  pause

  
  D(j) =  w^apower;
  dD(j,:) = apower*w^(apower-1)*dw;
  
%  tmp_p = w_p^apower;
%  tmp_m = w_m^apower; 
%  disp('Numerical')
%  disp((tmp_p-tmp_m)/2/h)
%  disp('Analytical')
%  disp(dD(j,2))
%  pause
  
  

  if (sum(isnan(dD(j,:))) > 0) % NaN values at the vertices of the polygon
    %disp('message 13')
    dD(j,:) = [0 0];
  end

end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,drho] = rhoweight(v1,v2,x)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative for the i-th edge 
% of the polygon using the triangle inequality
%

s1 = x - v1; s2 = x - v2; s = v2 - v1; d = norm(s);
xc = (v1+v2)/2;
f  = (s1(1)*s(2) - s1(2)*s(1))/d;
df = 1/d*[s(2), -s(1)];
t  = 1/d*( (d/2)^2 - sum((x-xc).^2) );
dt = -2/d*(x-xc);
rho = sqrt(f^2 + (sqrt(t^2+f^4)-t)^2 / 4 );
drho = 1/rho*(f*df + (sqrt(t^2+f^4)-t)/4* ((t*dt+2*f^3*df)/sqrt(t^2+f^4) - dt) );
%s1 = v1 - x; s2 = v2 - x; s = v2 - v1;
%rho = norm(s1) + norm(s2) - norm(s);
%drho = - s1/norm(s1) - s2/norm(s2);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formR(rho,drho,p)
% Purpose
% =======
% Compute the R function and its gradient
%
M = size(rho,1);

w=0;
dw=[0,0];
for i = 1:M
    w = w + (1/rho(i))^p;
    dw = dw + rho(i)^(-p-1)*drho(i,:);
end
dw = dw * w^(-1/p-1);
w  = w^(-1/p);

%w  = rho(1);
%dw = drho(1);
%for i = 2:M
%  w  = w + rho(i) - sqrt(w^2 + rho(i)^2);   
%  dw = dw + drho(i) - ( w*dw + rho(i)*drho(i) ) / sqrt(w^2 + rho(i)^2);   
%end

return
end
