function [xi, wg] = GaussLegendreCubature(order,dim)

if nargin==1
  dim = 1;
end

if dim==1
  n = max(ceil((order+1)/2),1);
  [xi,wg] = GaussLegendreCubature1D(n);
elseif dim==2
  [xi,wg] = GaussLegendreCubature2D(order);
else
  error('NOT IMPLEMENTED YET')
end