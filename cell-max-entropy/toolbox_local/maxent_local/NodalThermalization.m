function [beta_n,range_n,h_n]= NodalThermalization(x_n, options, id_bd)
% function [beta_n range_n h_n]= NodalThermalization(x_n, options, id_bd)
% Computation of the beta parameter and the range for searching nearest
% neighbours
% INPUT:
%   x_n     : pointset, (NxDim)
%   options : this structure must be contain the main settings for computing
%             the nearest neighbors, such that gamma, knn...
%   id_bd   : indices for the nodes on the boundary (or close to it), these will have
%             a special treatment
if nargin<3
  id_bd = [];
end

Tol0 = options.Tol0;
nPts = size(x_n,1);

beta_n = zeros(1,nPts);
range_n= zeros(1,nPts);

if isfield(options, 'radius')
  radius = options.radius;
  beta_n = -log(Tol0)/(radius^2) * ones(1,nPts);
  range_n=  radius * ones(1,nPts);
elseif isfield(options, 'knn')
  gamma = options.gamma;
  if options.knn==1
    h_n = NodalSpacing(x_n, options, id_bd);
  else
    if isfield(options,'spacing')
      if length(options.spacing)==nPts
        h_n = options.spacing;
      else
        h_n = options.spacing * ones(nPts,1);
      end
    else
      error('SPACING is needed (knn=0)')
    end
  end

  % beta and range are computed
  for i=1:nPts
    h_i       = h_n(i);
    beta_n(i) = gamma/(h_i*h_i);
    range     = sqrt(-log(Tol0)/gamma)*h_i;
    range_n(i)= max(range, 1.1*h_i);
  end

else
  error('RADIUS or KNN search criteria is needed');
end
