function [x_s,w_s] = MakeGLSamples1D(x_n,options)

nPts  = length(x_n);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Gauss-Legendre points and weights for a specific order are computed at
% the origin (reference)
order      = options.orderGL;
[pts,wgs]  = GaussLegendreCubature(order,1);
[~,sorted]= sort(x_n);
nel        = length(sorted)-1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Generation of samples for each node
gPts = length(pts);
sPts = nel*gPts;
x_s = zeros(sPts,1);
w_s = zeros(sPts,1);


k=1;
for iel=1:nel
  v1  = x_n(sorted(iel));
  len = x_n(sorted(iel+1)) - v1;
  for ig=1:gPts
     x_s(k) = v1 + (1-pts(ig))*len/2;
     w_s(k) = wgs(ig)*len/2;
     k = k+1;
  end
end

% figure(1)
% y_n = zeros(1, length(x_n));
% y_s = zeros(1, length(x_s));
% plot(x_n, y_n, 'ro', x_s, y_s, 'bx')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Pruning: the samples that fall outside of the domain are not included
if isfield(options, 'pruning')
  pruning = options.pruning;
else
  pruning = 0;
end

if pruning == 0
  return;
end

sPts= length(x_s);
ids = (1:sPts);

% limits calculation
range= options.range_n;
xMin = min(x_n)+pruning*range(1);
xMax = max(x_n)-pruning*range(nPts);

% final extraction
ids = ids( x_s<xMin | x_s>xMax );
x_s= x_s(ids);
w_s= w_s(ids);

% figure(2)
% y_n = zeros(1, length(x_n));
% y_s = zeros(1, length(x_s));
% plot(x_n, y_n, 'ro', x_s, y_s, 'bx')