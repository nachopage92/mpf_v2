%% Timoshenko cantilever beam - cell-based max-entropy

clear all
close all

format short
cd ../..
addpath(userpath_matlab_linux)
cd examples/ex_elasticbeam_e
addpath('ex_toolbox_local')

fprintf(1,'====================================================================\n');

%% CME options
%  Input parameters: CME options
TolNR = 1.e-14;  % Newton-Raphson tolerance fo the nonlinear max-ent problem

optCME.verb  = 0;       % 0:off    1:on
optCME.grad  = 1;       % Computation of the Gradient 0:OFF 1:ON
optCME.isConv= 1;       % flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)

optCME.nring = 0;       % number ring of neighbors to consider for each node: 1, 2, 3.
optCME.spow  = 0;       % w^spow, where w is the approximation to the distance function (R-function)
optCME.mpow  = 0;       % distance function derivative maximum degree of the approximation

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

for spow=3:4
  optCME.spow  = spow;
  for mpow=2:2
    optCME.mpow  = mpow;
    for rr=3:3
      optCME.nring = rr;
      
      disp('###############################################################################')
      fprintf(1,'\tCME  NRing=%d  spow=%d  mpow=%d\n', rr, spow, mpow);
      
      %% Convervence curves
      if optCME.nring == 2
        N_  = [3 5 9 17 33 65];
      else
        N_  = [5 9 17 33 65];
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
      orderGL= [15 10];
      gPtsV  = [54 25];
      optGL.quality = 0.01;
      
      for gp=1:length(orderGL)
        disp('=============================================================')
        %order 1 (1 gPts), 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
        optGL.orderGL = orderGL(gp);
        
        fprintf(1,'\tNumerical integration Gauss-Legendre  order=%2d  gPts=%2d\n',orderGL(gp), gPtsV(gp));
        
        for n=1:m
          disp('------------------------------------------------')
          
          if n>3
            optCME.TolNR = TolNR;
          else
            optCME.TolNR = 1.e-12;
          end
          % Node Points
          Ny = N_(n);
          Nx = N_(n)*8-3;
          hx = Lx/(Nx-1);
          hy = Ly/(Ny-1);
          fprintf(1,'Quasi Uniform grid %dx%d ---- h=[%5.3f %5.3f]\n', Nx, Ny, hx, hy);
          
          
          fileIN  = sprintf('../ex_elasticbeam/data_mesh/timo_mesh_N%03d_ring%d.mat', Ny, optCME.nring);
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
        end
        
        disp('---------------------------------------------------');
        
        %% Save data                 -------------------------------------------------------------
        nring = optCME.nring;
        
        if nu == 0.3
          fileOUT = sprintf('normL2_nuA/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts);
        else
          fileOUT = sprintf('normL2_nuB/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts);
        end
        fprintf(1,'Save output in file: %s\n',fileOUT);
        save(fileOUT,'Ndof','N_','normL2','nnzK','normE','spow','mpow','nring','gPts','parameters');
      end
    end
  end
end


return

%% Analysis of Results       -------------------------------------------------------------
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
fprintf ( 1, 'CME  : m=%7.4g  b=%7.4g\n', fit_num(1),  fit_num(2));

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
legend('FEM','CME','Location','Best','Orientation','horizontal')
title(strcat('Relative Error in Norm L_2 -- nring=',int2str(optCME.nring)),'fontsize',16,'fontweight','b')
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
