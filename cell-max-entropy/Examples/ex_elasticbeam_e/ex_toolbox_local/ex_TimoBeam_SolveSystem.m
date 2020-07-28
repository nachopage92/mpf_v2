function [u] = ex_TimoBeam_SolveSystem(dp_lme,s_near,x_nodes,x_samples,w_samples,parameters,options)
% function [u] =
% ex_TimoBeam_SolveSystem(dp_lme,s_near,x_nodes,x_samples,w_samples,parameters,options)
%
% This function assembled the stiffness matrix K and the right hand side
%     rhs corresponding to the cantilever beam problem explained in
%     Timoshenko's book.
%
% Input:
%    p_lme     : shape functions
%    s_lme     : shape functions gradient
%    s_near   : list of neighbors
%    x_nodes   : node points
%    x_samples : sample point
%    w_samples : gauss weigth for each sample point
%    parameters: L (length), D (diameter), nu (Poisson coefficient), E
%                (Young modulus)
%    options   : lme options
%
% Output:
%    u     : vectorial displacement field
%


% Material parameters
E  = parameters.E;
nu = parameters.nu;

nPts = size(x_nodes,1);
sPts = size(x_samples,1);


%% ------------------------------------------------------------------------
% The right hand side rhs is computed
[rhs,ind_Dirichlet] = ex_TimoBeam_RHS(x_nodes,options,parameters);


%% ------------------------------------------------------------------------
%  The stiffness matrix is assembled

optAssem.E       = E;
optAssem.nu      = nu;
optAssem.nPts    = nPts;
optAssem.s_near  = s_near;
optAssem.dp_samp = dp_lme;
optAssem.w_samp  = w_samples;
K = StiffnessMatrix_Assembly(optAssem);

%% Dirichlet BCs are applied

%(x=0,y=0)  ux=0 uy=0
ind = ind_Dirichlet(1);
K(2*ind,:)         = 0;
K(:,2*ind)         = 0;
rhs(2*ind)         = 0;
K(2*ind,2*ind)     = 1;
K(2*ind-1,:)       = 0;
K(:,2*ind-1)       = 0;
rhs(2*ind-1)       = 0;
K(2*ind-1,2*ind-1) = 1;

%(x=0,y=D/2)  ux=0
ind = ind_Dirichlet(2);
K(2*ind-1,:)       = 0;
K(:,2*ind-1)       = 0;
rhs(2*ind-1)       = 0;
K(2*ind-1,2*ind-1) = 1;



%% ------------------------------------------------------------------------
% The system is solved
u=K\rhs;

