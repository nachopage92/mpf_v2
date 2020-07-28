function [rhs,ind_bcD] = Beam_RHS_fem(x_nodes,options,parameters)
% function [rhs,ind_Dirichlet] = Beam_RHS_fem(x_nodes,options,parameters);
%
% This function computes the right hand side rhs
%
% Input
%    x_nodes   : node points
%    parameters: L (length), D (diameter), nu (Poisson coefficient), E
%                (Young modulus)
%    options   : lme options
%
% Output
%    rhs    : right hand side
%    ind_bcD: position of Dirichlet's nodes


% Number of nodes
nPts = length(x_nodes);
% The rhs is initialized
rhs    = zeros(2*nPts,1);

% Material parameters
P  = parameters.P;
D  = parameters.D;
L  = parameters.L;

I = D^3/12;

%% The BCs nodes are classified
ids   = (1:nPts);

% (x=0)
ind_bc1 = ids(abs(x_nodes(:,1))<1.e-6);
% (x=L)
ind_bc2 = ids(abs(x_nodes(:,1)-L)<1.e-6);
% (y=0)
ind_bc3 = ids(abs(x_nodes(:,2))<1.e-6);
%(x=0,y=0)
ind_bcD(1) = ids(abs(x_nodes(:,1))<1.e-6 & abs(x_nodes(:,2))<1.e-6);
%(x=0,y=D/2)
ind_bcD(2) = ids(abs(x_nodes(:,1))<1.e-6 & abs(x_nodes(:,2)-0.5*D)<1.e-6);

% options for 1D integration and 1D lme
opt_int1D.orderGL = 10;
opt_int1D.pruning = 0;
options1D      = options;
options1D.dim  = 1;
options1D.grad = 0;

%% Boundary nodes corresponding to x=0
y1D           = x_nodes(ind_bc1,2);
[y1D,ind_y1D] = sort(y1D);
ind_bc1       = ind_bc1(ind_y1D);
[y_s1D,w_s1D] = MakeGLSamples1D(y1D, opt_int1D);
options1D.beta_n  = options.beta_n(ind_bc1);
options1D.range_n = options.range_n(ind_bc1);

% adjacency structure with the nearest neighbors nodes to each sample point
%s_near1D = SamplesAdjacency(y1D,y_s1D,options1D.range_n);

sPts = length(y_s1D);
bcPts= length(y1D);
elem = [1:(bcPts-1);2:bcPts]';
nElem= size(elem,1);
gPts = sPts/nElem;
s_near1D = SamplesAdjacencyFEM(elem,gPts);

% Local-max entropy basis functions computation
options1D.s_near = s_near1D;
outLME  = wrapper_lme(y1D,y_s1D,options1D);
p_lme1D = outLME.p_samp;

rhs_x = zeros(bcPts,1);
rhs_y = zeros(bcPts,1);
for i=1:sPts
  y_s      = y_s1D(i);
  w_s      = w_s1D(i);
  i_nears  = s_near1D{i};
  p_nears  = p_lme1D{i};
  rhs_x(i_nears) = rhs_x(i_nears) + P*L/I*y_s*p_nears*w_s;
  rhs_y(i_nears) = rhs_y(i_nears) - P/(2*I)*(D^2/4-y_s^2)*p_nears*w_s;
end

rhs(2*ind_bc1-1) = rhs(2*ind_bc1-1) + rhs_x;
rhs(2*ind_bc1)   = rhs(2*ind_bc1)   + rhs_y;


%% Boundary nodes corresponding to x=L
y1D           = x_nodes(ind_bc2,2);
[y1D,ind_y1D] = sort(y1D);
ind_bc2       = ind_bc2(ind_y1D);
[y_s1D,w_s1D] = MakeGLSamples1D(y1D, opt_int1D);

options1D.beta_n  = options.beta_n(ind_bc2);
options1D.range_n = options.range_n(ind_bc2);

% adjacency structure with the nearest neighbors nodes to each sample point
%s_near1D = SamplesAdjacency(y1D,y_s1D,options1D.range_n);
sPts = length(y_s1D);
bcPts= length(y1D);

elem = [1:(bcPts-1);2:bcPts]';
nElem= size(elem,1);
gPts = sPts/nElem;

s_near1D = SamplesAdjacencyFEM(elem,gPts);


% Local-max entropy basis functions computation
options1D.s_near = s_near1D;
outLME  = wrapper_lme(y1D,y_s1D,options1D);
p_lme1D = outLME.p_samp;

rhs_y = zeros(bcPts,1);
for i=1:sPts
  y_s      = y_s1D(i);
  w_s      = w_s1D(i);
  i_nears  = s_near1D{i};
  p_nears  = p_lme1D{i};
  rhs_y(i_nears) = rhs_y(i_nears) + P/(2*I)*(D^2/4-y_s^2)*p_nears*w_s;
end

rhs(2*ind_bc2) = rhs(2*ind_bc2) + rhs_y;


%% Boundary nodes corresponding to y=0
x1D           = x_nodes(ind_bc3,1);
[x1D,ind_x1D] = sort(x1D);
ind_bc3       = ind_bc3(ind_x1D);
[x_s1D,w_s1D] = MakeGLSamples1D(x1D, opt_int1D);
options1D.beta_n  = options.beta_n(ind_bc3);
options1D.range_n = options.range_n(ind_bc3);

% adjacency structure with the nearest neighbors nodes to each sample point
%s_near1D = SamplesAdjacency(x1D,x_s1D,options1D.range_n);

sPts = length(x_s1D);
bcPts = length(x1D);

elem = [1:(bcPts-1);2:bcPts]';
nElem= size(elem,1);
gPts = sPts/nElem;

s_near1D = SamplesAdjacencyFEM(elem,gPts);


% Local-max entropy basis functions computation
options1D.s_near = s_near1D;
outLME  = wrapper_lme(x1D,x_s1D,options1D);
p_lme1D = outLME.p_samp;

rhs_x = zeros(bcPts,1);
for i=1:sPts
  w_s      = w_s1D(i);
  i_nears  = s_near1D{i};
  p_nears  = p_lme1D{i};
  rhs_x(i_nears) = rhs_x(i_nears) - P/(2*I)*D^2/4*p_nears*w_s;
end

rhs(2*ind_bc3-1) = rhs(2*ind_bc3-1) + rhs_x;



