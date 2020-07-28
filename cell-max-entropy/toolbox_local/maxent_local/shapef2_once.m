function [p_a dp_a hp_a outLME] = shapef2_once(dim,TolNR,n_near,beta,x_a,x)
% Calculation of the local maximum-entropy shape functions with Newton's
% method as described in section 4.2 of [1]
%
% INPUT
% =====
% dim:     spacial dimension
% TolNR:   tolerance for Newton's iterations
% x_a:      nearest nodes to x (n_near, dim)
% beta:     value of the thermalization parameter at each node point (1,n_near)
% x:        sample point where the shape functions are evaluated (1,dim)
%
% OUTPUT
% ======
% p_a:      p_a contains the values of the shape functions corresponding
%           to the closest neighbor nodes at the sample point x
% dp_a:     dp_a contains the values of the spatial derivatives of
%           the shape functions at the sample point x (gradient)
% hp_a:     hp_a contains the values of the second order spatial derivatives of
%           the shape functions at the sample point x (hessian)
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
verb= 'off';
err = 0; %error id
lam = zeros(1,dim);

R   = 10*ones(1,dim);
% dlam= 10*ones(1,dim);

%Newton iteration
niter=0;
%vector difference from the x_sample to the nears points in x_a
dx  = repmat(x,n_near,1)-x_a;

sum1    = sum(dx.^2,2);
exp_beta= exp(-beta'.*sum1);


while ( norm(R)>TolNR )
  [R,J,p_a] = Gamma2_(dim,exp_beta,dx,lam);
  if isnan(rcond(J))
    err = 1;
    break;
  end
  dlam = -(J\R)';
  lam  = lam + dlam;
  niter= niter+1;
  if (niter>100)
    err = 2;
    break;
  end
end

%check for error
if err == 1
  if strcmp(verb,'on')
    disp('LME2 :: Newton Failed 1 -- condition of X in the 1-norm is NaN');
  end
elseif err == 2
  if strcmp(verb,'on')
    disp('LME2 :: Newton Failed 2 -- maximum number of iterations reached (100)');
  end
else
  if strcmp(verb,'on') && niter>10
    fprintf(1,'LME2 :: Iterations reached: %d\n',niter);
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
beta_o= ones(dim,1)*beta;
pa_o  = (p_a*ones(1,dim))';
rb    = beta_o.*pa_o.*dx';
Jb = rb*dx;
rb = sum(rb,2);

rb = 2*rb;
Jb = 2*Jb;
id = eye(dim);

DL = (Jb-id)/J;
%DL = (Jb-id)*iJ;
Mx_a = 2*beta_o'.*dx - dx*DL';
gF_a = ones(n_near,1)*rb' - Mx_a;
dp_a = pa_o'.*gF_a;


%% Spatial Hessian matrix
if nargout < 3
  return
end

iJx= J\dx';
V  = dx*iJx + ones(n_near,n_near);
K  = {zeros(n_near,1)};
for ia=1:n_near
  K{ia} = p_a(ia)*Mx_a(ia,:)'*Mx_a(ia,:);
end

% % Hessian matrix computation
% hp_a = zeros((dim*dim+dim)/2,1);
hp_a = {zeros(dim,dim)};
bp   = beta*p_a;
Rb   = rb*rb';
for ia=1:n_near
  GJ  = zeros(dim,dim); %Contraction G*.grad(p_a*)
  for ib=1:n_near
    GJ = GJ + K{ib}*V(ia,ib);
  end % ib
  A = rb*iJx(:,ia)';
  
  hp_a{ia}(:,:) = p_a(ia)*( gF_a(ia,:)'*gF_a(ia,:) + ...
    2*(bp-beta(ia))*id - GJ + ...
    (rb'*iJx(:,ia))*id + ...
    A + A' + Rb );
end %ia