%% Timoshenko cantilever beam - Finite Element Method (P1: linear interpolants)
% plane strain state

clear all
close all

format short
addpath('ex_toolbox_local')

fprintf(1,'====================================================================\n');

%save data for this run -- check end lines after "Convervence curves"
saveData = 0; % 1: Yes     0: No

%% FEM options
optFEM.dim   = 2;
optFEM.verb  = 0; %0:off
optFEM.grad  = 1;           % Computation of the Gradient 0:OFF 1:ON
optFEM.hess  = 0;           % Computation of the Hessian  0:OFF 1:ON
optFEM.knn   = 0;

% Nodal Integration options
%cubature order: 1 (1 GPts), 2 (3 GPts), 4 (6 GPts), 5 (7 GPts), 6 (12 GPts)
optGL.orderGL = 6;
fprintf(1,'orderGL=%2d\n', optGL.orderGL);
  
%% Convervence curves
N_ = [81 120 283 581 1121 2132];% 4319 8510 16834];

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
  fprintf(1,'Unstructured set of points N=%d\n', N_(n));

  
  fileIN  = sprintf('data_mesh/mesh_N%05d_ring2.mat', N_(n));
  data    = load(fileIN);
  x_n     = data.x_n;
  tri     = data.tri;
  ibnd    = data.ibnd;
  
  nPts    = size(x_n,1);

  Ndof(n) = sqrt(nPts);
  
  % Sample points and gauss weights
  optGL.tri = tri;
  optGL.quality = 0.001;
  [x_s,w_s] = MakeGLSamples2D(x_n, optGL);
  sPts  = size(x_s,1);  %total number of Gauss points
  nElem = size(tri,1);        %number of elements
  gPts  = sPts/nElem;         %number of Gauss points in an element
  
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

  %% =====================================================================================
  t=cputime;

  optAssem.nPts    = nPts;
  optAssem.e_nears = e_nears;
  optAssem.p_samp  = p_samp;
  optAssem.dp_samp = dp_samp;
  optAssem.w_samp  = w_s;
  optAssem.x_samp  = x_s;

  % ------------------------------------------------------------------------
  %% The right hand side rhs is computed
  rhs = Unstr_RightHandSide(optAssem);

  % ------------------------------------------------------------------------
  %%  The stiffness matrix is assembled
  K = Unstr_StiffnessMatrix_Assembly(optAssem);
  nnzK(n) = nnz(K);

  % ------------------------------------------------------------------------
  %% Dirichlet BCs are applied
  %(x=0,y=0)  u=0
%   ids        = 1:nPts;
%   ind        = ids(abs(x_n(:,1))<1e-04 & abs(x_n(:,2))<1e-04);
%   K(ind,:)   = 0;
%   K(:,ind)   = 0;
%   K(ind,ind) = 1;
%   
%   rhs(ind)   = 0;
  
  % values on the boundary nodes is prescribed
  ids   = 1:nPts;
  id_bnd{1} = ids(abs(x_n(:,1)+0)<1.e-4); % (x=0)
  id_bnd{2} = ids(abs(x_n(:,1)-1)<1.e-4); % (x=1)
  id_bnd{3} = ids(abs(x_n(:,2)+0)<1.e-4); % (y=0)
  id_bnd{4} = ids(abs(x_n(:,2)-1)<1.e-4); % (y=1)
  for i=1:4
    id_segm= id_bnd{i};
    for j=1:length(id_segm)
      jj      = id_segm(j);
      x_bd    = x_n(jj,:);
      u_val   = Unstr_AnalyticalSolution(x_bd);
      rhs     = rhs - K(:,jj)*u_val;
      K(:,jj) = 0;
      K(jj,:) = 0;
      K(jj,jj)= 1;
      rhs(jj) = u_val;
    end
  end
  
  % ------------------------------------------------------------------------
  %% The displacement field is computed: system is solved
  u_num = K\rhs;

  fprintf(1,'\tcputime system: %4.3f\n', cputime-t);
  
  % ------------------------------------------------------------------------
  %% L2 norm and energy norm
  [normL2(n),normE(n)] = ErrorNorms(u_num,optAssem);
end

%% Save data
if saveData == 1
    fileOUT = sprintf('normL2/normL2_fem_GP%02d.mat',gPts);
    save(fileOUT,'Ndof','N_','normL2','normE','nnzK','gPts');
end

%% Analysis of Results
fit_num1 = polyfit(log10(1./Ndof(m_id)),log10(normL2(m_id)),1);
fit_num2 = polyfit(log10(1./Ndof(m_id)),log10(normE (m_id)),1);

disp('- - - - - - - - - - - - - - -');
fprintf ( 1, 'FEM  : normL2 m=%7.4g\n', fit_num1(1));
fprintf ( 1, '     : normE  m=%7.4g\n', fit_num2(1));

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
xlabel('DOF^{-1/2}','fontsize',14,'fontweight','b')
ylabel('L2 Error','fontsize',14,'fontweight','b')
legend('FEM','Location','Best','Orientation','horizontal')
title('Error in Norm L_2','fontsize',16,'fontweight','b')
set(gca,'XTick',0.02:0.02:0.1)
%xlim([0.001 0.11])

figure(2);clf
loglog(1./Ndof,normE,'ro-',h_line, err_line1,'k-', ...
  'LineWidth',2,'MarkerSize',6)
xlabel('DOF^{-1/2}','fontsize',14,'fontweight','b')
ylabel('Energy semi-norm','fontsize',14,'fontweight','b')
legend('FEM','Location','Best','Orientation','horizontal')
title('Error in Energy semi-norm','fontsize',16,'fontweight','b')
set(gca,'XTick',0.02:0.02:0.1)


%% Node and Sample Points
figure(3)
plot(x_n(:,1),x_n(:,2),'or', x_s(:,1),x_s(:,2),'.b');
axis equal
legend('Node points','Sample points');


% Numerical and analytical solutions
u_sol = Unstr_AnalyticalSolution(x_n);

figure(4) 
plot3(x_n(:,1),x_n(:,2),u_sol,'.r')
hold on
trimesh(tri,x_n(:,1),x_n(:,2),u_num)
hold off
xlabel('X')
ylabel('Y')
