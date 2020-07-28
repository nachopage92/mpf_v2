function plotLME2D(x_n,u_n, nSX, nSY, Lx, Ly, centre, options)
% function PlotLME2D(x_n,u_n, nSX, nSY, Lx, Ly, centre, options)
% this functios plots the LME's shape functions in 2D
% x_n: nodes coordinates
% u_n: scalar field for each node to be interpolated
% nSX: number of samples in X
% nSY: number of samples in Y
% Lx, Ly are the size of the rectangular region and center its center.

x_s = UniformGrid2D(nSX, nSY, Lx, Ly, centre);
hs    = Lx/(nSX-1);
x_samp= (0:hs:Lx)-Lx/2+centre(1);
y_samp= (0:hs:Ly)-Ly/2+centre(2);
sPts  = length(x_s);

U_mat = zeros(nSX,nSY);
id_smp= (1:sPts);

[p_s nn_s] = wrapper_lme2D(x_n,x_s,options);

u_s= zeros(sPts,1);
for k=1:sPts
  u_s(k) = u_n(nn_s{k})'*p_s{k};
end
k=1;
for i=1:nSX
  for j=1:nSY
    ij  = id_smp(k);
    U_mat(i,j)= u_s(ij);
    k=k+1;
  end
end

axis off
plot3(x_n(:,1),x_n(:,2),u_n,'ko','LineWidth',2,...
  'MarkerFaceColor','g','MarkerSize',6);
hold on
surf(x_samp,y_samp,U_mat, ... %'LineWidth',0.5,...
  'FaceColor','red','FaceLighting','phong')
view([10 40])
hold off
%material shiny
%lighting gouraud
lightangle(80, -40)
lightangle(-90, 60)


%% Wrapper to compute the LOCAL MAX ENTROPY of first order (only the shape functions)
function [p_lme s_nears] = wrapper_lme2D(x_n,x_s,options)
% INPUT:
%   x_n: node coordinates
%   x_s: samples coordinates
%   options: parameters setting
%
% OUTPUT:
%   p_lme  : are the shape functions values for each sample point
%   s_nears: samples adjacency structure

if isfield(options, 'dim')
  dim = options.dim;
else
  error('Variable DIM has not be defined in options');
end

n_s   = size(x_s,1);
if n_s==1
  x_s=[x_s;x_s];
  n_s = 2;
end

% samples adjacency structure is computed
beta_n  = options.beta_n;
range_n = options.range_n;
if isfield(options, 's_nears')
  s_nears = options.s_nears;
else
  s_nears = SamplesAdjacency(x_n,x_s,range_n);
end

% LME Shape functions and spatial derivatives up to second order computation
p_lme = {zeros(n_s,1)};

TolNR = options.TolNR;  %Newton-Raphson tolerance
for i=1:n_s
  nn_ids = s_nears{i}; %nearest neighbour indices list
  n_near = length(nn_ids);
  x_near = x_n(nn_ids,:);

  %beta_n have the values of beta_a for each neighbor of xsample(i)
  beta_a = beta_n(nn_ids);
  
  % shape functions, gradient and hessian
  if std(beta_a/mean(beta_a)) < 1.e-8
      p_lme{i} = shapef_once( dim,TolNR,n_near,beta_a,x_near,x_s(i,:));
  else
      p_lme{i} = shapef2_once(dim,TolNR,n_near,beta_a,x_near,x_s(i,:));
  end
end