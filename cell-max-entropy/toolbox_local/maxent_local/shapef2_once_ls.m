function [p_a dp_a hp_a outLME] = shapef2_once_ls(dim,TolNR,n_near,beta,x_a,x)
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
verb = 'off';
err  = 0; %error id
lam = zeros(1,dim);
R   = 10*ones(1,dim);

%Newton iteration
niter=0;
%vector difference from the x_sample to the nears points in x_a

dx       = repmat(x,n_near,1)-x_a;
sum1     = sum(dx.^2,2);
exp_beta = exp(-beta'.*sum1);
step     = norm(beta'.*sum1)/norm(dx);

% options for the line search
opts=optimset('TolX',1.e-08,'MaxIter',5000);


% fprintf(1,'step=%f\n',step);
% graf=[];
while ( norm(R)>TolNR )
  [R,J,p_a] = Gamma2_(dim,exp_beta,dx,lam);
  if isnan(rcond(J))
    err = 1;
    break;
  end
  dlam = -(J\R)';
  
  val = norm(dlam);
  if isnan(val)
    dlam = lam*0.01;
  end
  if val>step
    dlam = 0.5*step*dlam/val;
  end

  % line search
  t = fminbnd(@(t) Gamma2_ls(t,exp_beta,dx,lam,dlam),0,1,opts);
  lam  = lam + t*dlam;
  niter= niter+1;
%   graf(niter)=val;

  if (niter>1000)
    err = 2;
    break;
  end
end

% check for error
if err == 1
  p_a = [];
  if strcmp(verb,'on')
    disp('LME2:LS::Newton Failed 1 -- condition of X in the 1-norm is NaN');
  end
elseif err == 2
  p_a = [];
  if strcmp(verb,'on')
      disp('LME2:LS::Newton Failed 2 -- maximum number of iterations reached (1000)');
  end
else
  if strcmp(verb,'on') && niter>20
    fprintf(1,'LME2:LS::Iterations reached: %d\n',niter);
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
dp_a = zeros(n_near,dim);

Mx_a = zeros(n_near,dim);
gF_a = zeros(n_near,dim);

rb = 0.0;
Jb = zeros(dim,dim);
for ia=1:n_near
  rb = rb + beta(ia)*p_a(ia)*dx(ia,:)';
  Jb = Jb + beta(ia)*p_a(ia)*dx(ia,:)'*dx(ia,:);
end
rb = 2*rb;
Jb = 2*Jb;
id = eye(dim);
% iJ = inv(J);
DL = (Jb-id)/J;

if dim==1
  for ia=1:n_near
    dp_a(ia,:)= p_a(ia)* ( rb - (2*beta(ia)-DL)*dx(ia) );
  end
else
  for ia=1:n_near
    Mx_a(ia,:)= (2*beta(ia)*id - DL)*dx(ia,:)';
    gF_a(ia,:)= rb' - Mx_a(ia,:);
    dp_a(ia,:)= p_a(ia)* gF_a(ia,:);
  end
end

%% Spatial Hessian matrix
if nargout < 3
  return
end
hp_a = {zeros(dim,dim)};

iJx= J/dx';
bp = 0.0;
V  = zeros(n_near, n_near);
K  = {zeros(n_near,1)};
for ia=1:n_near
  bp    = bp + beta(ia)*p_a(ia);
  V(ia,ia) = dx(ia,:)*iJx(:,ia) + 1;
  for ib=ia+1:n_near
    V(ia,ib) = dx(ia,:)*iJx(:,ib) + 1;
    V(ib,ia) = V(ia,ib);
  end
  K{ia} = p_a(ia)*Mx_a(ia,:)'*Mx_a(ia,:);
end

% Hessian matrix
Rb = rb*rb';
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