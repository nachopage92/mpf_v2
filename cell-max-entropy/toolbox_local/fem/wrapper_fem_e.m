function output = wrapper_fem_e(x_n,x_s,options)
% INPUT:
%   x_n: node coordinates
%   x_s: node coordinates
%   options: parameters setting
%
% OUTPUT:
%   N_fem  : are the shape functions values for each sample point
%   dN_fem : gradient for the shape functions at the sample points
%

if nargout > 3
  error('The number of output is too big');
end
if isfield(options, 'verb')
  verb = options.verb;
else
  verb = 1;
end
if isfield(options, 'grad')
  grad = options.grad;
else
  grad = 0;
end

t=cputime;

sPts   = size(x_s,1);

%% samples adjacency structure is computed
e_nears= options.e_nears;
nElem  = length(e_nears);
gPts   = sPts/nElem;

%% FEM Shape functions and spatial derivatives up to second order computation
N_fem  = {zeros(sPts,1)};
if grad==1
  dN_fem = {zeros(sPts,1)};
else
  dN_fem = {[]};
end

k=1;
for e=1:nElem
  nn_ids = e_nears{e}; %nearest neighbour indices list
  x_near = x_n(nn_ids,:);
  
  for g=1:gPts
    if grad==1
      % shape functions, gradient
      [N_fem{k},dN_fem{k}] = TriangleBasisFunctions(x_near,x_s(k,:));
    else
      % shape functions
      N_fem{k} = TriangleBasisFunctions(x_near,x);
    end
    k = k + 1;
  end
end

if verb==1
  format compact
  fprintf(1,'\tCPUTIME LME: %f\n', cputime-t);
end

output.p_samp  = N_fem;
output.dp_samp = dN_fem;
