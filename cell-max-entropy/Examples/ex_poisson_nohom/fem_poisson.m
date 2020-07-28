%% EXAMPLE: Heat Equation - Finite Element Method with P1

clear all
close all

format short
addpath('ex_toolbox_local')

%% LME options (Machine precision)
optFEM.dim   = 2;
optFEM.verb  = 0; % 0:off    1:on
optFEM.knn   = 0;
optFEM.grad  = 1;           % Computation of the Gradient 0:OFF 1:ON

%% Integration options
%Numerical integration, cubature order
optGL.orderGL = 4;  %order 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
optGL.quality = 0.01;

%% Geometrical Parameters
L     = 1.0;    % length
Ndof  = [5 9 17 33 65 129];

normL2= zeros(length(Ndof),1);
centre= L*[0.5,0.5];

for n=1:length(Ndof)
  disp('------------------------------------------------')
  
  % Node Points
  N     = Ndof(n);
  x_nodes = UniformGrid2D(N, N, L, L, centre);
  nPts  = length(x_nodes);
  h     = L/(N-1);
  
  fprintf(1,'\tGrid size: %2dx%2d\n',N,N);
  
  ids   = 1:nPts;
  id_bnd{1} = ids(abs(x_nodes(:,1)+0)<1.e-6); % (x=-L/2)
  id_bnd{2} = ids(abs(x_nodes(:,1)-L)<1.e-6); % (x= L/2)
  id_bnd{3} = ids(abs(x_nodes(:,2)+0)<1.e-6); % (y=-L/2)
  id_bnd{4} = ids(abs(x_nodes(:,2)-L)<1.e-6); % (y= L/2)
  
  % GAUSS_LEGENDRE
  [x_samples,w_samples,tri] = MakeGLSamples2D(x_nodes, optGL);
  sPts  = size(x_samples,1);  %total number of Gauss points
  nElem = size(tri,1);        %number of elements
  gPts  = sPts/nElem;         %number of Gauss points in an element
  
  fprintf(1,'\torderGL=%2d\n', optGL.orderGL);
  
  
  %% The basis functions are computed
  % adjacency structure with the nearest neighbors nodes to each sample point
  % nodal shape parameter
  t=cputime;
  s_near = SamplesAdjacencyFEM(tri,gPts);
  
  % Local-max entropy basis functions computation
  optFEM.s_near = s_near;
  outFEM = wrapper_fem(x_nodes,x_samples,optFEM);
  
  p_samp  = outFEM.p_samp;
  dp_samp = outFEM.dp_samp;
  fprintf(1,'\tcputime FEM    : %4.3f\n', cputime-t);
  
  
  %% ------------------------------------------------------------------------
  % System ASSEMBLY
  t = cputime;
  optSystem.nPts    = nPts;
  optSystem.s_near  = s_near;
  optSystem.p_samp  = p_samp;
  optSystem.dp_samp = dp_samp;
  optSystem.w_samp  = w_samples;
  
  % RHS computation
  f = zeros(nPts,1); 
  
  % ASSEMBLY: The stiffness matrix K is filled
  K = NoHom_StiffnessMatrix_Assembly(optSystem);
  
  %%  -----------------------------------------------------------------------
  % The temperature value on a line segment is enforced

  % values on the boundary nodes is prescribed
  for i=1:4
    id_segm= id_bnd{i};
    for j=1:length(id_segm)
      jj      = id_segm(j);
      x_bd    = x_nodes(jj,:);
      u_val   = NoHom_AnalyticalSolution(1,x_bd);
      f       = f - K(:,jj)*u_val;
      K(:,jj) = 0;
      K(jj,:) = 0;
      K(jj,jj)= 1;
      f(jj)   = u_val;
    end
  end
  
  %% Compute System Solution
  u_fem = K\f;
  
  fprintf(1,'\tcputime System : %4.3f\n', cputime-t);
  
  %% Error in L2 norm
  t = cputime;
  errL2 = 0;
  for k =1:sPts
    k_near = s_near{k};
    u_k    = u_fem(k_near);
    p_k    = p_samp{k};
    
    u_h   = p_k'*u_k;
    u     = NoHom_AnalyticalSolution(1,x_samples(k,:));
    errL2 = errL2 + (u-u_h)^2 * w_samples(k);
  end
  normL2(n) = sqrt(errL2);
  fprintf(1,'\tcputime L2 norm: %4.3f\n', cputime-t);
end
disp('---------------------------------------------------');

%% Plot results
t = cputime;

%save(sprintf('normL2/normL2_fem_GP%02d.mat',gPts),'Ndof','normL2','gPts');

n_id  = (n-2):n;
mrate = polyfit(log(1./(Ndof(n_id)-1)),log(normL2(n_id)'),1);

mref  = polyfit(log([0.1,0.05]),log([0.0005,0.000125]),1);
fprintf(1,'\tRate of convergence  m=%4.2f  mref=%4.2f\n', mrate(1), mref(1));

% Error in L2 norm
figure(1);clf
loglog(1./(Ndof-1),normL2,'ro-',[0.1,0.05],[0.0005,0.000125],'k-','linewidth',2)
title('FEM  Error in L2 norm  |u-u_h|_2','fontsize',18,'fontweight','b')
xlabel('DOF', 'fontsize',16,'fontweight','b')
ylabel('|u-u_h|_2', 'fontsize',16,'fontweight','b')

% Analytical solution
u_sol = NoHom_AnalyticalSolution(1,x_nodes);

% Numerical solution
figure(2);clf
%trisurf(tri,x_nodes(:,1),x_nodes(:,2),u_fem)
plot3(x_nodes(:,1),x_nodes(:,2),u_fem,'ro')
hold on
plot3(x_nodes(:,1),x_nodes(:,2),u_sol,'b*')
hold off
title('FEM  P_1','fontsize',18,'fontweight','b')
xlabel('X', 'fontsize',16,'fontweight','b')
ylabel('Y', 'fontsize',16,'fontweight','b')
zlabel('Temperature','fontsize',16,'fontweight','b')
fprintf(1,'\tcputime Plot  : %4.3f\n', cputime-t);