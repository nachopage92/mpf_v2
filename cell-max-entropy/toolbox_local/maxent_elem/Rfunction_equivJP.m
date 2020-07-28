function [D, dD] = Rfunction_equivJP(mesh,v,p,spow,mpow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% =======
% Evaluate level set function using R-functions (equivalence): rho^1 ^ rho^2 ^
% rho^3 . . . . ^ rho^n-1 where rho^i is the edge-weight for the i-th edge
% USAGE: [D, dD] = Rfunction_equiv(v,p,in)
%
% Inputs:
% v  : [x1 y1; x2 y2; ...; xn yn], n vertices of the polygon in ccw (Counter Clock Wise); 
%      polygon is open (n vertices) or closed (n+1 vertices; repeat 1st vertex)
% p  : [ px(1) py(1); ... ; px(m) py(m) ], the m points for which we want 
%      to compute the levet set function and its derivatives
%
% Outputs:
% D  : level set function = [D_1; ...; D_m] (m x 1 array)
% dD : gradient of level set function = [dD_1; ...; dD_m] (m x 2 array)
%
% N. Sukumar: UC Davis, M. Arroyo UPC-BarcelonaTech
% D. Millan UPC-BarcelonaTech, February 2014 

N    = size(mesh,1);  % N edges
m    = size(p,1);
rho  = zeros(N,1); % N edges
drho = zeros(N,2); % N edges


ids_q1= unique(mesh(:,1));
ids_q2= unique(mesh(:,2));
ids_q = intersect(ids_q1,ids_q2);
v_q  = v(ids_q,:);
Nq   = size(ids_q,1);  % N points

qjp  = zeros(Nq,1); % N edges
dqjp = zeros(Nq,2); % N edges


D    = zeros(m,1);
dD   = zeros(m,2);

for j = 1:m
   x = p(j,:);
  
  % Compute the nonnegative edge weight and its derivative
  [rho,drho] = rhoweight(mesh,v,x,rho,drho);
  
  if isempty(ids_q)==0
    %normalized field measuring the distance to the joining point: qjp
    [qjp,dqjp] = joinweight(v_q,x,qjp,dqjp);
  
    % Compute the R function and its gradient
    [w,dw]  = formR_JP(rho,drho,qjp,dqjp,mpow);
  else
    % Compute the R function and its gradient
    [w,dw]  = formR(rho,drho,mpow);
  end
  
  plotrho0(mesh,v,p,x,rho,dw);
  
  D(j)    = w^spow;
  dD(j,:) = spow*w^(spow-1)*dw;

  if (sum(isnan(dD(j,:))) > 0) % NaN values at the vertices of the polygon
%     mesh
%     rho
%     drho
%     w
%     dw
%     disp('NaN values at the vertices of the polygon')
    dD(j,isnan(dD(j,:))) = 0;
%     pause
  end

end

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,drho] = rhoweight(mesh,v,x,rho,drho)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative
% of the polygon using the triangle inequality
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
%   if abs(f)<2*eps  %Spacing of floating point numbers: 2^(-52) ~ 2.22e-16
%     rho (i)   = 0;
%     drho(i,:) = [s(2), -s(1)]*id;
%   else
    df = [s(2), -s(1)]*id;
    t  = ( (0.5*d)^2 - sum((x-xc).^2) )*id;
    dt = -2*(x-xc)*id;

    val = sqrt(t^2+f^4);
    %val = sqrt(t^2);
%    val = sqrt(t^4+f^2);

    rho (i)   = sqrt(f^2 + 0.25*(val-t)^2 );
     
    if abs(rho(i))<2*eps  %Spacing of floating point numbers: 2^(-52) ~ 2.22e-16
      drho(i,:) = -df;
    else
      drho(i,:) = (f*df + 0.25*(val-t) * ((t*dt+2*f^3*df)/val - dt) )/rho(i);
      %drho(i,:) = (f*df + 0.25*(val-t) * ((t*dt+f*df/4)/val - dt) )/rho(i);
    end
%   end
end
%fprintf(1,'\td=%4.2f  f=%4.2f  t=%4.2f  rho[%d]=%4.2f\n', d, f, t, N, rho(N));

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qjp,dqjp] = joinweight(v_q,x,qjp,dqjp)
% Purpose
% =======
% Compute the nonnegative edge weight and its derivative
% of the polygon using the triangle inequality
%
Nq    = size(v_q,1);  % N points
for i = 1:Nq
  vi    = v_q(i,:);
  dxv   = norm(x-vi);
  
  qjp(i)= dxv;
  if dxv < 2*eps
    dqjp(i,:) = [0 0];
  else
    dqjp(i,:) = (x-vi)/dxv;
  end
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formR(rho,drho,p)
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
  ws  = sum(1./(rho.^p));
  w   = ws^(-1/p);
  dw  = ( (rho.^(-p-1))' * drho ) * w/ws;
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,dw] = formR_JP(rho,drho,qjp,dqjp,p)
% Purpose
% =======
% Compute the R function and its gradient
%

% N    = length(rho);
% ids  = 1:N;
% ids  = ids(rho<eps);
% nrho0= length(ids);
% 
% if nrho0>0
%   w = 0;
%   if nrho0 == 1
%     dw = drho(ids,:);
%   else
%     dw = mean(drho(ids,:));
%   end
% else
%   if (sum(1./(rho.^p)) - sum(1./qjp.^p.*exp(-qjp.^p))) <0
%     sum(1./(rho.^p))
%     sum(1./qjp.^p.*exp(-qjp.^p))
%     pause
%   end
%   
%   ws  = sum(1./(rho.^p)) - sum(1./qjp.^p.*exp(-qjp.^p));
%   w   = ws^(-1/p);
%   dw  = ( (rho.^(-p-1))' * drho - ...
%         ( (qjp.^(-p-1) + qjp.^(-1)).*exp(-qjp.^p))' * dqjp) * w/ws;
% end
  ws  = sum(1./(rho.^p)) - sum(1./qjp.^p);
  w   = ws^(-1/p);
  dw  = ( (rho.^(-p-1))' * drho - ...
          (qjp.^(-p-1))' * dqjp ) * w/ws;
 
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotrho0(mesh,v,p,x,rho,dw)
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
  plot(v(:,1),v(:,2),'r.')
  hold on
  for i=1:nrh0
    ii = mesh(ids0(i),:);
    plot(v(ii,1),v(ii,2),'ro-','markersize',8)
  end
  plot(x(1),x(2),'cx','markersize',14,'linewidth',2)
  hold off
  axis equal
  fprintf(1,'Rfunction_equivJP: rho have %d zero values at x=[%4.2f %4.2f]\n',nrh0,x(1),x(2));
  
  figure(31);clf
  plot(v(:,1),v(:,2),'ro-','markersize',12,'linewidth',2)
  hold on
  plot(p(:,1),p(:,2),'bx',x(1),x(2),'bs')
  quiver(x(1),x(2),dw(1),dw(2))
  hold off
  axis equal
  pause
end

return
end