%% Timoshenko cantilever beam - local max-entropy

clear all
close all

format short
addpath('ex_toolbox_local')

fprintf(1,'====================================================================\n');

%% LME options
optLME.dim   = 2;
optLME.verb  = 0; %0:off
optLME.grad  = 1;           % Computation of the Gradient 0:OFF 1:ON
optLME.hess  = 0;           % Computation of the Hessian  0:OFF 1:ON
optLME.TolNR = 1.e-12;
optLME.knn   = 0;

%% Matrial Parameters
D  = 1;         % height
L  = 4*D;       % length

P  = -1000;     % external load
E  =  1e07;     % Young modulus
nu =   0.300;   % Poisson coefficient: 0.3, 0.499 (quasi incomp)

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

tpow  = 6;%[4,6,8];
Tol0V = 1e-6;%[1e-4,1e-6,1e-8];
gammaV= [1.8,4.8];
for gg=1:length(gammaV)
  optLME.gamma = gammaV(gg);      % value of gamma to compute LME
  
  % Nodal Integration options
  %cubature order: 1 (1 GPts), 2 (3 GPts), 4 (6 GPts), 5 (7 GPts), 6 (12 GPts)
  if optLME.gamma == 1.8
    optGL.orderGL = 6;
  elseif optLME.gamma == 4.8
    optGL.orderGL = 4;
  else
    optGL.orderGL = 2;
  end
  
  disp('###########################################################################')
  fprintf(1,'Aspect shape parameter gamma=%4.2f\n', gammaV(gg));
  for it=1:length(Tol0V)
    optLME.Tol0  = Tol0V(it);
    disp('========================================================')
    fprintf(1,'Tolerance %8.2e\n', Tol0V(it));
    for n=1:m
      disp('------------------------------------------------')
      % Node Points
      Ny = N_(n);
      hy = Ly/(Ny-1);
      
      Nx = Ny*8-3;
      hx = Lx/(Nx-1);
      fprintf(1,'Quasi Uniform grid %dx%d ---- h=[%5.3f %5.3f]\n', Nx, Ny, hx, hy);
      
      fileIN  = sprintf('data_mesh/timo_mesh_N%03d_ring2.mat', Ny);
      data    = load(fileIN);
      x_n     = data.x_n;
      tri     = data.tri;
      ibnd    = data.ibnd;
      
      nPts    = length(x_n);
      
      h       = 0.5*(hx+hy);
      
      Ndof(n) = sqrt(nPts*2);
      
      fprintf(1,'\tn=%d  Ny=%2d   nPts=%4d   nodal spacing=%6.4f\n',n, Ny, nPts, h);
      
      % Sample points and gauss weights
      optGL.tri = tri;
      [x_s,w_s] = MakeGLSamples2D(x_n, optGL);
      sPts  = size(x_s,1);    %total number of Gauss points
      nElem = size(tri,1);    %number of elements
      gPts  = sPts/nElem;     %number of Gauss points in an element
      
      t=cputime;
      fprintf(1,'\tgamma=%4.2f | order=%2d\n',optLME.gamma,optGL.orderGL);
      
      %% Nodes thermalization, nodes and samples adjacency structures are computed
      optLME.spacing= h;
      optLME.h_n    = ones(nPts,1)*h;

      beta_n         = optLME.gamma/(h*h)*ones(1,nPts);
      range_n        = max(1.1,sqrt(-log(optLME.Tol0)/optLME.gamma))*optLME.h_n;
      optLME.beta    = optLME.gamma/(h*h);
      optLME.beta_n  = beta_n;
      optLME.range_n = range_n;
      
      disp([min(range_n) mean(range_n) max(range_n)]/h)
      
      %% The basis functions are computed
      % adjacency structure with the nearest neighbors nodes to each sample point
      % nodal shape parameter
      n_nears = SamplesAdjacency(x_n,x_n,range_n);
      
      % Samples adjacency
      e_nears={[]};
      for e=1:nElem
        e_nears{e} = unique([n_nears{tri(e,:)}]);
      end
      clear n_nears;
      fprintf(1,'\tcputime e_nears  : %5.2f\n', cputime-t);
      
      t=cputime;
      % Local-max entropy basis functions computation
      optLME.e_nears = e_nears;
      outLME = wrapper_lme_e(x_n,x_s,optLME);
      
      p_samp  = outLME.p_samp;
      dp_samp = outLME.dp_samp;
      fprintf(1,'\tcputime LME      : %5.2f\n', cputime-t);
      
      t=cputime;
      % ------------------------------------------------------------------------
      %% The right hand side rhs is computed
      [rhs,ind_Dirichlet] = Beam_RHS_lme(x_n,optLME,parameters);
      
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
    
    disp('---------------------------------------------------');
    
    %% Save data
    gamma = optLME.gamma;
    Tol0  = optLME.Tol0;
    
    if nu==0.3
      fileOUT = sprintf('normL2_nuA/normL2_lme_g%3.1f_GP%02d_Tol1e%d.mat',gamma,gPts,tpow(it));
    else
      fileOUT = sprintf('normL2_nuB/normL2_lme_g%3.1f_GP%02d_Tol1e%d.mat',gamma,gPts,tpow(it));
    end
    %save(fileOUT,'Ndof','N_','normL2','normE','nnzK','gamma','Tol0','gPts','parameters');
  end
end

%return

%% Analysis of Results
if nu==0.3
  data_fem   = load('normL2_nuA/normL2_fem_GP03.mat');
else
  data_fem   = load('normL2_nuB/normL2_fem_GP03.mat');
end
Ndof_fem   = data_fem.Ndof;
normL2_fem = data_fem.normL2;
normE_fem  = data_fem.normE;

fit_num  = polyfit(log10(1./Ndof(m_id)),log10(normL2(m_id)),1);
disp('- - - - - - - - - - - - - - -');
fprintf ( 1, 'LME  : m=%7.4g  b=%7.4g\n', fit_num(1),  fit_num(2));

% Plots
h_line=[0.01 0.1]*2;
err_line=[0.0001 0.01]*2;
fit_line=polyfit(log10(h_line),log10(err_line),1);
fprintf ( 1, 'line : m=%7.4g  b=%7.4g\n', fit_line(1),fit_line(2));

% Error in L2 norm
figure(1);clf
loglog(1./Ndof_fem,normL2_fem,'bs-','markersize',10,'linewidth',2)
hold on
loglog(1./Ndof,normL2,'ro-',h_line, err_line,'k-', ...
  'LineWidth',2,'MarkerSize',6)
hold off
xlabel('h','fontsize',14,'fontweight','b')
ylabel('Relative Error','fontsize',14,'fontweight','b')
legend('FEM','LME','Location','Best','Orientation','horizontal')
title(strcat('Relative Error in Norm L_2 -- \gamma_{LME}=',num2str(optLME.gamma)),'fontsize',16,'fontweight','b')
set(gca,'XTick',0.02:0.02:0.1)
%xlim([0.001 0.11])


% Node and Sample Points
figure(2)
plot(x_n(:,1),x_n(:,2),'or',...
  x_s(:,1),x_s(:,2),'.b');
axis equal
legend('Node points','Sample points');


% Deformed and undeformed (numerical and analytical) configurations
figure(3)
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
