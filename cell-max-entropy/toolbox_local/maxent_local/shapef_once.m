function [p_a dp_a hp_a outLME] = shapef_once(dim,TolNR,n_near,beta,x_a,x)
% Calculation of the local maximum-entropy shape functions with Newton's
% method as described in section 4.2 of [1]
%
% INPUT
% =====
% dim  :   spacial dimension
% x_a  :   nearest nodes to x (n_near, dim)
% beta :   value of the thermalization parameter at each node point
% x    :   sample point where the shape functions are evaluated
% TolNR:   tolerance for Newton's iterations
%
% OUTPUT
% ======
% p_a :    p_a contains the values of the shape functions corresponding
%          to the closest neighbor nodes at the sample point x
% dp_a:    dp_a contains the values of the dim spacial derivatives of
%          the shape functions corresponding to the closest neighbor nodes at
%          the sample point x (gradient)
% hp_a:    hp_a contains the values of the dim second order spacial derivatives of
%          the shape functions corresponding to the closest neighbor nodes at
%          the sample point x (hessian)
%
% outLME.lam   = lam    : Lagrange multiplier
% outLME.err   = err    : error identifier
% outLME.niter = niter  : number of Newton-Raphson iterations
% 
% Reference:
% [1] Marino Arroyo and Michael Ortiz, "Local maximum-entropy approximation
%     schemes: a seamless bridge between finite elements and meshfree methods",
%     International Journal for Numerical Methods in Engineering, 65:2167-2202 (2006).

warning('off','All');
lam = zeros(1,dim);
R   = 10*ones(1,dim);
val = 1;
verb= 'off';
err = 0; %error id

%Newton iteration
niter=0;
%vector difference from the x_sample to the nears points in x_a
dx  = repmat(x,n_near,1)-x_a;

sum1= sum(dx.^2,2);
exp_beta= exp(-beta'.*sum1);

step = norm(beta'.*sum1)/norm(dx);

while (norm(R)>TolNR && val>TolNR)
  [R,J,p_a] = Gamma2_(dim,exp_beta,dx,lam);
  if isnan(rcond(J))
    err = 1;
    break;
  end
  dlam = -(J\R)';
  
  val = norm(dlam);
  if isnan(val)
    disp('DLAM is NaN - - LME')
    dlam = lam*0.01;
  end
  if val>step 
      dlam = 0.5*step*dlam/val;
  end
  lam  = lam + dlam;
  niter= niter+1;
  if (niter>100)
    err = 2;  
    break;
  end
end

%check for error
if err == 1
  disp(det(J))
  disp(norm(R))
  disp(norm(lam))
  if strcmp(verb,'on')
    disp('LME :: Newton Failed 1 -- condition of X in the 1-norm is NaN');
  end
elseif err == 2
  if strcmp(verb,'on')
    disp('LME :: Newton Failed 2 -- maximum number of iterations reached (100)');
  end
else
  if strcmp(verb,'on') && niter>10
    fprintf(1,'LME :: Iterations reached: %d\n',niter);
  end
end

outLME.lam   = lam;
outLME.err   = err;
outLME.niter = niter;

if err > 0
  p_a  = [];
  dp_a = [];
  hp_a = [];
  return
end

if nargout < 2
  return
end

%% Spacial Gradients
iJx=J\dx';

pa_o =  p_a*ones(1,dim);
dp_a = -pa_o.*iJx';

%% % Spatial Hessian matrix
if nargout < 3
  return
end

V  = dx*iJx + ones(n_near,n_near);
K = {zeros(n_near,1)};
for ia=1:n_near
  K{ia} = -iJx(:,ia)*dp_a(ia,:);
end


%% Hessian matrix
hp_a = {zeros(dim,dim)};
for ia=1:n_near
  GJ= zeros(dim,dim); %Contraction iJ*G*.grad(p_a*)*iJ
  for ib=1:n_near
    GJ = GJ + K{ib}*V(ia,ib);
  end % ib
  hp_a{ia} = K{ia} - p_a(ia)*GJ;   
end %ia