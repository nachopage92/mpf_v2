%% Timoshenko cantilever beam - cell-based max-entropy

clear all
close all

format short
addpath('ex_toolbox_local')

fprintf(1,'====================================================================\n');

%% CME options
%  Input parameters: CME options
TolNR = 1.e-14;  % Newton-Raphson tolerance fo the nonlinear max-ent problem

optCME.verb  = 0;       % 0:off    1:on
optCME.grad  = 1;       % Computation of the Gradient 0:OFF 1:ON
optCME.isConv= 1;       % flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)

optCME.nring = 3;       % number ring of neighbors to consider for each node: 1, 2, 3.
optCME.spow  = 3;       % w^spow, where w is the approximation to the distance function (R-function)
optCME.mpow  = 2;       % distance function derivative maximum degree of the approximation

%% Matrial Parameters
D  = 1;          % height
L  = 4*D;% length

P  = -1000;     % external load
E  =  1e07;     % Young modulus
nu =  0.300;    % Poisson coefficient: 0.3, 0.499 (quasi incomp)

% sigma = C_stiff*epsilon
C_stiff = E/(1+nu)/(1-2*nu)*[1-nu,   nu,          0 ;...
                               nu, 1-nu,          0 ;...
                                0,    0, (1-2*nu)/2];

%epsilon = S_stiff*sigma
S_stiff = (1+nu)/E*[1-nu,  -nu, 0;
                     -nu, 1-nu, 0;
                       0,    0, 2];
parameters.L  = L;
parameters.D  = D;
parameters.E  = E;
parameters.nu = nu;
parameters.P  = P;
Lx = L;
Ly = 0.5 * D;


disp('###############################################################################')
fprintf(1,'\tCME  NRing=%d  spow=%d  mpow=%d\n', optCME.nring, optCME.spow, optCME.mpow);

%% Convervence curves
if optCME.nring == 2
  N_  = [3 5 9 17];% 33 65];
else
  N_  = [5 9 17];% 33 65];
end
m      = length(N_);

m_id   = (m-2:m);      %index of the L2 error norm to perform the fit
h_n    = zeros(m,1);
Ndof   = zeros(m,1);

normL2 = zeros(m,1);
normE  = zeros(m,1);
nnzK   = zeros(m,1);

% -------------------------------------------------------------------
% Integration options
% Numerical integration, cubature order
% cubature order: 1 (1 GPts), 2 (3 GPts), 4 (6 GPts), 5 (7 GPts), 6 (12 GPts)
orderGL= 6;%10;
gPtsV  = 12;%25;
optGL.quality = 0.01;

disp('=============================================================')
%order 1 (1 gPts), 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
optGL.orderGL = orderGL;

fprintf(1,'\tNumerical integration Gauss-Legendre  order=%2d  gPts=%2d\n',orderGL, gPtsV);

for n=1:m
  disp('------------------------------------------------')
  
  if n>2
    optCME.TolNR = TolNR;
  else
    optCME.TolNR = 1.e-10;
  end
  % Node Points
  Ny = N_(n);
  Nx = N_(n)*8-3;
  hx = Lx/(Nx-1);
  hy = Ly/(Ny-1);
  fprintf(1,'Quasi Uniform grid %dx%d ---- h=[%5.3f %5.3f]\n', Nx, Ny, hx, hy);
  
  
  fileIN  = sprintf('data_mesh/timo_mesh_N%03d_ring%d.mat', Ny, optCME.nring);
  data    = load(fileIN);
  nPts    = length(data.x_n);
  
  h       = 0.5*(hx+hy);
  h_n(n)  = h;
  Ndof(n) = sqrt(nPts*2);
  
  fprintf(1,'\tn=%d  Ny=%2d   nPts=%4d   nodal spacing=%6.4f\n',n, Ny, nPts, h);
  
  x_n     = data.x_n;
  tri     = data.tri;
  ibnd    = data.ibnd;
  connecStr = data.connecStr;
  
  % Sample points and gauss weights
  optGL.tri = tri;
  [x_s,w_s] = MakeGLSamples2D(x_n, optGL);
  sPts  = size(x_s,1);  %total number of Gauss points
  nElem = size(tri,1);        %number of elements
  gPts  = sPts/nElem;         %number of Gauss points in an element
  
  fprintf(1,'\tnring=%d | s=%d  m=%d | order=%2d\n',optCME.nring,optCME.spow,optCME.mpow,optGL.orderGL);
  
  %% =======================================================================================
  %  CME: basis functions computation
  t=cputime;
  
  optCME.node_spacing  = h;
  optCME.tri           = tri;
  optCME.bd_segm       = connecStr.bd_segm;
  optCME.bd_segm_out   = connecStr.bd_segm_out;
  optCME.bd_nodes      = connecStr.bd_nodes;
  optCME.nodes_in_elem = connecStr.nodes_in_elem;
  
  outCME  = wrapper_cme_e(x_n,x_s,optCME);
  
  p_samp  = outCME.p_samp;
  dp_samp = outCME.dp_samp;
  e_nears = outCME.e_nears;
  
  fprintf(1,'\tcputime CME basis functions    BULK: %5.3f seconds\n', cputime-t);
  
  
  %% =======================================================================================
  % -------------------------------------------------------------------
  % System matrix ASSEMBLY
  t=cputime;
  % ------------------------------------------------------------------------
  %% The right hand side rhs is computed
  [rhs,ind_Dirichlet] = Beam_RHS_cme(x_n,optCME,parameters);
  
  % ------------------------------------------------------------------------
  %%  The stiffness matrix is assembled
  optAssem.E       = E;
  optAssem.nu      = nu;
  optAssem.Cmat    = C_stiff;
  optAssem.Smat    = S_stiff;
  optAssem.P       = P;
  optAssem.L       = L;
  optAssem.D       = D;
  optAssem.nPts    = nPts;
  optAssem.e_nears = e_nears;
  optAssem.p_samp  = p_samp;
  optAssem.dp_samp = dp_samp;
  optAssem.w_samp  = w_s;
  K = Beam_StiffnessMatrix_Assembly(optAssem);
  nnzK(n) = nnz(K);
  
  % ------------------------------------------------------------------------
  %% Dirichlet BCs are applied
  %(x=0,y=0)  ux=0 uy=0
  ind = ind_Dirichlet(1);
  K(2*ind,:)         = 0;
  K(:,2*ind)         = 0;
  K(2*ind,2*ind)     = 1;
  K(2*ind-1,:)       = 0;
  K(:,2*ind-1)       = 0;
  K(2*ind-1,2*ind-1) = 1;
  
  rhs(2*ind)         = 0;
  rhs(2*ind-1)       = 0;
  
  %(x=0,y=D/2)  ux=0
  ind = ind_Dirichlet(2);
  K(2*ind-1,:)       = 0;
  K(:,2*ind-1)       = 0;
  K(2*ind-1,2*ind-1) = 1;
  
  rhs(2*ind-1)       = 0;
  
  % ------------------------------------------------------------------------
  %% The displacement field is computed: system is solved
  u_h = K\rhs;
  
  fprintf(1,'\tcputime system: %4.3f\n', cputime-t);
  
  % ------------------------------------------------------------------------
  %% L2 norm and energy norm
  [normL2(n),normE(n)] = Beam_ErrorNorms(u_h,x_s,optAssem);
  fprintf(1,'\t|u-u_h| L2 norm=%9.3e    Energy L2 error=%9.3e\n',normL2(n),normE(n));
end

disp('---------------------------------------------------');

%% Save data                 -------------------------------------------------------------
nring = optCME.nring;
spow  = optCME.spow;
mpow  = optCME.mpow;

if nu == 0.3
  fileOUT = sprintf('normL2_nuA/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts);
else
  fileOUT = sprintf('normL2_nuB/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts);
end
fprintf(1,'Save output in file: %s\n',fileOUT);
%         save(fileOUT,'Ndof','N_','normL2','nnzK','normE','spow','mpow','nring','gPts','parameters');


%% Analysis of Results       -------------------------------------------------------------
if nu==0.3
  data_fem   = load('normL2_nuA/normL2_fem_GP03.mat');
else
  data_fem   = load('normL2_nuB/normL2_fem_GP03.mat');
end
Ndof_fem   = data_fem.Ndof;
normL2_fem = data_fem.normL2;
normE_fem  = data_fem.normE;

if nu==0.3
  data_lme   = load('normL2_nuA/normL2_lme_g0.8_GP12_Tol1e6.mat');
  data_cme   = load(sprintf('normL2_nuA/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts));
else
  data_lme   = load('normL2_nuB/normL2_lme_g0.8_GP12_Tol1e6.mat');
  data_cme   = load(sprintf('normL2_nuB/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts));
end
Ndof_lme   = data_lme.Ndof;
normL2_lme = data_lme.normL2;
normE_lme  = data_lme.normE;

Ndof_cme   = data_cme.Ndof;
normL2_cme = data_cme.normL2;
normE_cme  = data_cme.normE;

fit_fem_L2 = polyfit(log10(1./Ndof_fem(m_id)),log10(normL2_fem(m_id)),1);
fit_fem_E  = polyfit(log10(1./Ndof_fem(m_id)),log10(normE_fem(m_id)),1);
fit_lme_L2 = polyfit(log10(1./Ndof_lme(m_id)),log10(normL2_lme(m_id)),1);
fit_lme_E  = polyfit(log10(1./Ndof_lme(m_id)),log10(normE_lme(m_id)),1);
fit_L2     = polyfit(log10(1./Ndof(m_id)),log10(normL2(m_id)),1);
fit_E      = polyfit(log10(1./Ndof(m_id)),log10(normE(m_id)),1);

disp('- - - - - - - - - - - - - - -');
h_line=[0.01 0.1]*2;
err_line1=[0.00001 0.001]*2;
err_line2=[0.020 0.2]*2;
fit_line1=polyfit(log10(h_line),log10(err_line1),1);
fit_line2=polyfit(log10(h_line),log10(err_line2),1);


fprintf ( 1, 'FEM  normL2: m=%7.3g\n', fit_fem_L2(1));
fprintf ( 1, 'LME  normL2: m=%7.3g\n', fit_lme_L2(1));
fprintf ( 1, 'CME  normL2: m=%7.3g\n', fit_L2(1));
fprintf ( 1, 'ref slope1 : m=%7.3g\n', fit_line1(1));

fprintf ( 1, 'FEM  normE : m=%7.3g\n', fit_fem_E(1));
fprintf ( 1, 'LME  normE : m=%7.3g\n', fit_lme_E(1));
fprintf ( 1, 'CME  normE : m=%7.3g\n', fit_E(1));
fprintf ( 1, 'ref slope2 : m=%7.3g\n', fit_line2(1));

format shortg
% disp([normL2 normL2_lme(2:4) normE normE_lme(2:4)])

% Plots ================================

% Error in L2 norm
figure(1);clf
loglog(1./Ndof_fem,normL2_fem,'go-','markersize',12,'markerfacecolor',[0,0.5,0],'linewidth',2)
hold on
loglog(1./Ndof_lme,normL2_lme,'bo--','markersize',12,'markerfacecolor',[0,0.4,0.8],'linewidth',2)
loglog(1./Ndof_cme,normL2_cme,'ro-','markersize',12,'markerfacecolor',[1,0.4,0.4],'linewidth',2)
loglog(1./Ndof,normL2,'rs-','markersize',12,'linewidth',2)
loglog(h_line,err_line1,'k-','linewidth',2)
hold off
title('Error in L2 norm  |u-u_h|_2','fontsize',18,'fontweight','b')
xlabel('DOF^{-1/2}', 'fontsize',16,'fontweight','b')
ylabel('L2 error', 'fontsize',16,'fontweight','b')
legend(...
  'FEM P1',...
  'LME \gamma=0.8',...
  sprintf('CME r=%d s=%d m=%d GP%02d - OLD',nring,spow,mpow,gPts),...
  sprintf('CME r=%d s=%d m=%d GP%02d - NEW',nring,spow,mpow,gPts),...
  'ref slope=2')



% Error in Energy norm
figure(2);clf
loglog(1./Ndof_fem,normE_fem,'go-','markersize',12,'markerfacecolor',[0,0.5,0],'linewidth',2)
hold on
loglog(1./Ndof_lme,normE_lme,'bo--','markersize',12,'markerfacecolor',[0,0.4,0.8],'linewidth',2)
loglog(1./Ndof_cme,normE_cme,'ro-','markersize',12,'markerfacecolor',[1,0.4,0.4],'linewidth',2)
loglog(1./Ndof,normE,'rs-','markersize',12,'linewidth',2)
loglog(h_line,err_line2,'k-','linewidth',2)
hold off
title('Energy error in L2 norm','fontsize',18,'fontweight','b')
xlabel('DOF^{-1/2}', 'fontsize',16,'fontweight','b')
ylabel('Energy error', 'fontsize',16,'fontweight','b')
legend(...
  'FEM P1',...
  'LME \gamma=0.8',...
  sprintf('CME r=%d s=%d m=%d GP%02d - OLD',nring,spow,mpow,gPts),...
  sprintf('CME r=%d s=%d m=%d GP%02d - NEW',nring,spow,mpow,gPts),...
  'ref slope=1')

return
%% Node and Sample Points
figure(3)
plot(x_n(:,1),x_n(:,2),'or',...
  x_s(:,1),x_s(:,2),'.b');
axis equal
legend('Node points','Sample points');


% Deformed and undeformed (numerical and analytical) configurations
figure(4)
u_sol = Beam_AnalyticalSolution(x_n,100*P,D,L,E,nu);
u_sol = u_sol+x_n;
u_num = reshape(u_h,2,length(x_n))';
u_num = u_num+x_n;
plot(x_n(:,1),x_n(:,2),'.r')
axis equal
hold on
plot(u_num(:,1),u_num(:,2),'xb')
hold on
plot(u_sol(:,1),u_sol(:,2),'.g')
hold on
legend('Undeformed','Deformed (numerical)','Deformed (analytical)');
