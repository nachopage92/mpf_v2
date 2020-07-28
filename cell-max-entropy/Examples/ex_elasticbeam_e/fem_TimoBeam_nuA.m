%% Timoshenko cantilever beam - Finite Element Method (P1: linear interpolants)
% plane strain state

clear all
close all

format short
addpath('ex_toolbox_local')


fprintf(1,'====================================================================\n');

%% FEM options
optFEM.dim   = 2;
optFEM.verb  = 0; %0:off
optFEM.grad  = 1;           % Computation of the Gradient 0:OFF 1:ON
optFEM.hess  = 0;           % Computation of the Hessian  0:OFF 1:ON
optFEM.knn   = 0;

% Nodal Integration options
%cubature order: 1 (1 GPts), 2 (3 GPts), 4 (6 GPts), 5 (7 GPts), 6 (12 GPts)
optGL.orderGL = 2;
fprintf(1,'orderGL=%2d\n', optGL.orderGL);
  
%% Matrial Parameters
D  = 1;         % height
L  = 4*D;% length

P  = -1000;     % external load
E  =  1e07;     % Young modulus
nu =   0.3;     % Poisson coefficient: 0.3, 0.499 (quasi incomp)

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

%% Convervence curves
N_     = [3 5 9 17 33];

m      = length(N_);
m_id   = (m-2:m);      %index of the L2 error norm to perform the fit
h_n    = zeros(m,1);
Ndof   = zeros(m,1);
normL2 = zeros(m,1);
normE  = zeros(m,1);
nnzK   = zeros(m,1);

for n=1:m
  disp('------------------------------------------------')
  % Node Points
  Ny = N_(n);
  hy = Ly/(Ny-1);

  Nx = N_(n)*8-3;
  hx = Lx/(Nx-1);
  fprintf(1,'Quasi Uniform grid %dx%d ---- h=[%5.3f %5.3f]\n', Nx, Ny, hx, hy);

  
  fileIN  = sprintf('data_mesh/timo_mesh_N%03d_ring2.mat', Ny);
  data    = load(fileIN);
  x_n     = data.x_n;
  tri     = data.tri;
  ibnd    = data.ibnd;
  
  nPts    = length(x_n);

  h       = 0.5*(hx+hy);
  h_n(n)  = h;
  Ndof(n) = sqrt(nPts*2);
  
  fprintf(1,'\tn=%d  Ny=%2d   nPts=%4d   nodal spacing=%6.4f\n',n, Ny, nPts, h);
  
  optFEM.spacing= h;
  optFEM.h_n    = ones(nPts,1)*h;

  % Sample points and gauss weights
  optGL.tri = tri;
  [x_s,w_s] = MakeGLSamples2D(x_n, optGL);
  sPts  = size(x_s,1);  %total number of Gauss points
  nElem = size(tri,1);        %number of elements
  gPts  = sPts/nElem;         %number of Gauss points in an element

 
  % The RHS computes on each boundary line LME as FEM
  optFEM.TolNR = 1.e-12;
  optFEM.knn   = 0;
  optFEM.Tol0  = 1.e-06;
  optFEM.gamma = 5.8;      % value of gamma to compute LME
  optFEM.spacing= h;
  optFEM.h_n    = ones(nPts,1)*h;
  [beta_n,range_n] = NodalThermalization(x_n, optFEM);
  optFEM.beta   = optFEM.gamma/(h*h);
  optFEM.beta_n = beta_n;
  optFEM.range_n= range_n;

  
  %% The basis functions are computed
  % adjacency structure with the nearest neighbors nodes to each sample point
  % nodal shape parameter
  t=cputime;
  e_nears = SamplesAdjacencyFEM(tri,1);
  
  % Local-max entropy basis functions computation
  optFEM.e_nears = e_nears;
  outFEM = wrapper_fem_e(x_n,x_s,optFEM);
  
  p_samp  = outFEM.p_samp;
  dp_samp = outFEM.dp_samp;
  fprintf(1,'\tcputime FEM    : %4.3f\n', cputime-t);

  t=cputime;
  % ------------------------------------------------------------------------
  %% The right hand side rhs is computed
  [rhs,ind_Dirichlet] = Beam_RHS_fem(x_n,optFEM,parameters);

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
end

%% Save data
if nu==0.3
  fileOUT = sprintf('normL2_nuA/normL2_fem_GP%02d.mat',gPts);
else
  fileOUT = sprintf('normL2_nuB/normL2_fem_GP%02d.mat',gPts);
end
%save(fileOUT,'Ndof','N_','normL2','normE','nnzK','gPts','parameters');

%return

%% Analysis of Results
fit_num  = polyfit(log10(h_n(m_id)),log10(normL2(m_id)),1);
disp('- - - - - - - - - - - - - - -');
fprintf ( 1, 'FEM  : m=%7.4g  b=%7.4g\n', fit_num(1),  fit_num(2));

% Plots
h_line=[0.01 0.1]*2;
err_line1=[0.0050 0.05]*2;
err_line2=[0.0001 0.01]*2;
fit_line1=polyfit(log10(h_line),log10(err_line1),1);
fit_line2=polyfit(log10(h_line),log10(err_line2),1);
fprintf ( 1, 'line : mref1=%7.4g   mref2=%7.4g\n', fit_line1(1),fit_line2(1));

figure(1);clf
loglog(1./Ndof,normL2,'ro-',h_line, err_line2,'k-', ...
  'LineWidth',2,'MarkerSize',6)
xlabel('DOF^{-1/2]','fontsize',14,'fontweight','b')
ylabel('L2 Error','fontsize',14,'fontweight','b')
legend('FEM','Location','Best','Orientation','horizontal')
title('Error in Norm L_2','fontsize',16,'fontweight','b')
set(gca,'XTick',0.02:0.02:0.1)
%xlim([0.001 0.11])

figure(2);clf
loglog(1./Ndof,normE,'ro-',h_line, err_line1,'k-', ...
  'LineWidth',2,'MarkerSize',6)
xlabel('DOF^{-1/2]','fontsize',14,'fontweight','b')
ylabel('Energy semi-norm','fontsize',14,'fontweight','b')
legend('FEM','Location','Best','Orientation','horizontal')
title('Error in Energy semi-norm','fontsize',16,'fontweight','b')
set(gca,'XTick',0.02:0.02:0.1)


%% Node and Sample Points
figure(3)
plot(x_n(:,1),x_n(:,2),'or', x_s(:,1),x_s(:,2),'.b');
axis equal
legend('Node points','Sample points');


% Deformed and undeformed (numerical and analytical) configurations
figure(4)
u_sol = Beam_AnalyticalSolution(x_n,P,D,L,E,nu);
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
