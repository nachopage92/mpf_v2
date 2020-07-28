function [Dd fd]= DisplacementEBC(x_n,id_bd,options)
% function [Dd fd]= DisplacementEBC(x_n,id_bd,options)
%
% Input:
%    parameters:
%    options   : lme options (Tol0, TolNR, gamma)
%
% Output:
%    Dbc   : matrix to impose the Displacement Boundary Condition
%

nPts = length(x_n);

%% plate boundaries and modelization options

% options for 1D integration and 1D lme
options1D       = options;
options1D.dim   = 1;
optGL1D.orderGL = 6;
options1D.gamma = 4.8;%options.gammaLag;

options.dim     = 2;

dPts = length(id_bd.D);

%Stiffness matrix for enforcing the essential boundary conditions
Dd   = zeros(dPts,nPts);
fd   = zeros(dPts,1);

id_disp = setdiff(id_bd.D,id_bd.G);

x_d1D   = x_n(id_disp,1);
x_d1D   = sort(x_d1D);

[x_s1D w_s1D] = MakeGLSamples1D(x_d1D, optGL1D);
sPts = length(x_s1D);

% nodal thermalization and range search for the nodes on the boundary with LME options
% options1D
x_d1D   = x_n(id_bd.D,1);
x_d1D   = sort(x_d1D);
[beta_n1D range_n1D] = NodalThermalization(x_d1D, options1D);
options1D.beta_n     = beta_n1D;
options1D.range_n    = range_n1D;

%samples adjacency structure is computed
b_nears1D        = SamplesAdjacency(x_d1D,x_s1D,range_n1D);
options1D.s_near = b_nears1D;
options1D.grad   = 0;
options1D.hess   = 0;

outLME1D = wrapper_lme(x_d1D,x_s1D,options1D);
N_s1D = outLME1D.p_samp;
%END 1D ===

% the samples on the boundary are expressed in 2D global coordinates for computing the
% shape functions (zero for interior nodes), and gradients
x_s2D = [x_s1D  zeros(sPts,1)];

% plot(x_n(:,1),x_n(:,2),'ro',x_s2D(:,1),x_s2D(:,2),'b*',...
%   x_n(id_bd.D,1),x_n(id_bd.D,2),'k*')
% xlabel('X')
% ylabel('Y')
% axis equal
% pause

% LME's shape functions are evaluated in the samples on the boundary
options.adjac_s= SamplesAdjacency(x_n,x_s2D,options.range_n);
p_s2D     = wrapper_lme(x_n,x_s2D,options);
s_nears2D = options.adjac_s;

% 'cross' tangential vector (see Krysl-Belytschko 1996)
for k=1:sPts
  b_near = b_nears1D{k};      % nearest nodes on the boundary with rotational EBC
  k_near = s_nears2D{k};      % nearest nodes in the domain
  n_b    = length(b_near);
  n_k    = length(k_near);
  N_b    = N_s1D{k};
  p_k    = p_s2D{k};
  w_gk   = w_s1D(k);
  
  % Attention, this way to perform the ASSEMBLY of A is very bad, and for real runs it
  % needs to be improved for doing the common operations of FEM (gather,scatter)
  for ia=1:n_b
    i = b_near(ia);
    for jb=1:n_k
      j = k_near(jb);
      Dd(i,j) = Dd(i,j) - N_b(ia)*p_k(jb)*w_gk;
    end
    fd(i) = fd(i) - N_b(ia)*w_gk;
  end
end