function output = wrapper_lme_e(x_n,x_s,options)
% INPUT:
%   x_n: node coordinates
%   x_s: node coordinates
%   options: parameters setting
%
% OUTPUT:
%   p_lme  : are the shape functions values for each sample point
%   dp_lme : gradient for the shape functions at the sample points
%   hp_lme : hessian
%   id_samp: samples where the basis functions have been computed OK

if nargout > 3
  error('The number of output is too big');
end
if isfield(options, 'verb')
  verb = options.verb;
else
  verb = 1;
end
if isfield(options, 'hess')
  hess = options.hess;
else
  hess = 0;
end
if isfield(options, 'grad')
  grad = options.grad;
else
  grad = 0;
end

t=cputime;

if isfield(options, 'dim')
  dim = options.dim;
else
  error('Variable DIM has not be defined in options');
end

sPts   = size(x_s,1);

%% samples adjacency structure is computed
beta_n  = options.beta_n;
e_nears = options.e_nears;
nElem   = length(e_nears);
gPts    = sPts/nElem;

%% LME Shape functions and spatial derivatives up to second order computation
p_lme  = {zeros(sPts,1)};
if grad==1
  dp_lme = {zeros(sPts,1)};
else
  dp_lme = {[]};
end
if hess==1
  hp_lme = {zeros(sPts,1)};
else
  hp_lme = {[]};
end

id_samp= zeros(1,sPts);

TolNR = options.TolNR;  %Newton-Raphson tolerance

k = 1;
for e=1:nElem
  nn_ids = e_nears{e}; %nearest neighbour indices list
  n_near = length(nn_ids);
  x_near = x_n(nn_ids,:);
  
  %beta_n have the values of beta_a for each neighbor of xsample(i)
  beta_a = beta_n(nn_ids);

  for g=1:gPts    
    if hess==1
      % shape functions, gradient and hessian
      [p_lme{k},dp_lme{k},hp_k] = shapef2_once(dim,TolNR,n_near,beta_a,x_near,x_s(k,:));
      
      %if there is an error then LineSearch is done
      if isempty(dp_lme{k})==1
        if verb==1
          disp('err :=> line search')
        end
        [p_lme{k},dp_lme{k},hp_k] = shapef2_once_ls(dim,TolNR,n_near,beta_a,x_near,x_s(k,:));
      end
      
      if dim==2
        hp_voigt = zeros(n_near,3);
        for ii=1:n_near
          hp_voigt(ii,1) =   hp_k{ii}(1,1);
          hp_voigt(ii,2) =   hp_k{ii}(2,2);
          hp_voigt(ii,3) = 2*hp_k{ii}(1,2);
        end
        hp_lme{k} = hp_voigt;
      else
        hp_voigt = zeros(n_near,1);
        for ii=1:n_near
          hp_voigt(ii,1) =   hp_k{ii}(1,1);
        end
        hp_lme{k} = hp_voigt;
      end
      
    elseif grad==1
      % shape functions, gradient
      [p_lme{k},dp_lme{k}] = shapef2_once(dim,TolNR,n_near,beta_a,x_near,x_s(k,:));
      
      %if there is an error then LineSearch is done
      if isempty(dp_lme{k})==1
        if verb==1
          disp('err :=> line search')
        end
        [p_lme{k},dp_lme{k}] = shapef2_once_ls(dim,TolNR,n_near,beta_a,x_near,x_s(k,:));
      end
    else
      % shape functions, gradient and hessian
      [p_lme{k}] = shapef2_once(dim,TolNR,n_near,beta_a,x_near,x_s(k,:));
      
      %if there is an error then LineSearch is done
      if isempty(p_lme{k})==1
        if verb==1
          disp('err :=> line search')
        end
        [p_lme{k}] = shapef2_once_ls(dim,TolNR,n_near,beta_a,x_near,x_s(k,:));
      end
    end
    
    if sum(isnan(p_lme{k})) == 0 && isempty(p_lme{k})==0
      id_samp(k)=k;
    end
    k = k + 1;
  end
end
id_samp(id_samp==0)=[];

if verb==1
  format compact
  fprintf(1,'\tCPUTIME LME: %f\n', cputime-t);
end

output.id_samp = id_samp;
output.p_samp  = p_lme;
output.dp_samp = dp_lme;
output.hp_samp = hp_lme;
