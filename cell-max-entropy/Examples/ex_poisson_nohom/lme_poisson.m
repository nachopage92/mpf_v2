%% EXAMPLE: Heat Equation - Local max-ent basis functions

clear all
close all

format short

addpath('ex_toolbox_local')

disp('------------------------------------------------')

%save data for this run -- check end lines of "Convervence curves"
saveData = 0; % 1: Yes     0: No

%% =======================================================================================
%  Input parameters: LME options

optLME.gamma = 0;
optLME.dim   = 2;
optLME.verb  = 0; % 0:off    1:on
optLME.knn   = 0;
optLME.grad  = 1;           % Computation of the Gradient 0:OFF 1:ON
optLME.hess  = 0;
optLME.Tol0  = 1.e-06;
optLME.TolNR = 1.e-14;      %Newton-Raphson tolerance fo the nonlinear max-ent problem

gammaV = [4.8, 1.8, 0.8];

% -------------------------------------------------------------------
% Integration options
%Numerical integration, cubature order
%order 1 (1 gPts), 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
orderGL= [2 4  6 10]; %[1 2 3 4 5  6 10];
gPtsV  = [3 6 12 25]; %[1 3 4 6 7 12 25];
optGL.quality = 0.01;

% -------------------------------------------------------------------
% Geometrical Parameters
L     = 1.0;    % length
Ndof  = [5 9 17 33 65]; %129]

normL2= zeros(length(Ndof),1);
centre= L*[0.5,0.5];

for ig=1:length(gammaV)
  optLME.gamma = gammaV(ig);
  for gp=1:length(orderGL)
    disp('=============================================================')
    %order 1 (1 gPts), 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
    optGL.orderGL = orderGL(gp);
    
    fprintf(1,'\tNumerical integration Gauss-Legendre  order=%2d  gPts=%2d\n',orderGL(gp), gPtsV(gp));
    
    for n=1:length(Ndof)
      disp('------------------------------------------------')
      t=cputime;
      
      % Node Points
      N       = Ndof(n);
      x_nodes = UniformGrid2D(N, N, L, L, centre);
      nPts    = length(x_nodes);
      h       = L/(N-1);
      h_n     = ones(nPts,1)*h;
      
      optLME.spacing= h;
      optLME.h_n    = h_n;
      
      fprintf(1,'\tGrid size: %2dx%2d\n',N,N);
      
      ids   = 1:nPts;
      id_bnd{1} = ids(abs(x_nodes(:,1)+0)<1.e-6); % (x=-L/2)
      id_bnd{2} = ids(abs(x_nodes(:,1)-L)<1.e-6); % (x= L/2)
      id_bnd{3} = ids(abs(x_nodes(:,2)+0)<1.e-6); % (y=-L/2)
      id_bnd{4} = ids(abs(x_nodes(:,2)-L)<1.e-6); % (y= L/2)
      
      % -------------------------------------------------------------------
      % Gauss Legendre quadrature points for numerical integration
      [x_samples,w_samples,tri] = MakeGLSamples2D(x_nodes, optGL);
      sPts  = size(x_samples,1);  %total number of Gauss points
      nElem = size(tri,1);        %number of ellments
      
      % -------------------------------------------------------------------
      % Boundary nodes identifiers
      ebnd   = freeBoundary(triangulation(tri,x_nodes));
      ibnd   = unique(ebnd);
      bPts   = length(ibnd);
      x_bd   = x_nodes(ibnd,:);
      
      %plot triangulation constrained by the boundary edges
      figure(1);clf
      hold on
      plot(x_nodes(:,1),x_nodes(:,2), ...
        'ko','MarkerFaceColor','b','Markersize',10 )
      hold on
      for k = 1:size(ebnd,1)
        plot(x_nodes(ebnd(k,:),1),x_nodes(ebnd(k,:),2),'r-','LineWidth',4)
      end
      plot(x_bd(:,1),x_bd(:,2), ...
        'ko','MarkerFaceColor','c','Markersize',8 )
      triplot(tri,x_nodes(:,1),x_nodes(:,2),'m-','LineWidth',1);
      axis equal
      
      fprintf(1,'\tcputime LME Connectivity Structures: %5.3f seconds\n', cputime-t);
      
      
      %% =======================================================================================
      %  LME: basis functions computation
      t=cputime;
      
      % Nodes thermalization, nodes and samples adjacency structures are computed
      [beta_n,range_n] = NodalThermalization(x_nodes, optLME);
      optLME.beta   = optLME.gamma/(h*h);
      optLME.beta_n = beta_n;
      optLME.range_n= range_n;
      
      disp([min(range_n) mean(range_n) max(range_n)])
      
      % The basis functions are computed
      % adjacency structure with the nearest neighbors nodes to each sample point
      % nodal shape parameter
      s_near = SamplesAdjacency(x_nodes,x_samples,range_n);
      
      % Local-max entropy basis functions computation
      optLME.s_near = s_near;
      
      outLME  = wrapper_lme(x_nodes,x_samples,optLME);
      
      p_samp  = outLME.p_samp;
      dp_samp = outLME.dp_samp;
      
      fprintf(1,'\tcputime LME basis functions    BULK: %5.3f seconds\n', cputime-t);
      
      
      %% =======================================================================================
      
      % -------------------------------------------------------------------
      % System matrix ASSEMBLY
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
      
      %  -----------------------------------------------------------------------
      % values on the boundary nodes is prescribed
      for i=1:3
        id_segm= id_bnd{i};
        for j=1:length(id_segm)
          jj      = id_segm(j);
          x_bd    = x_nodes(jj,:);
          u_val   = 0;
          f       = f - K(:,jj)*u_val;
          K(:,jj) = 0;
          K(jj,:) = 0;
          K(jj,jj)= 1;
          f(jj)   = u_val;
        end
      end
      
      %Least squares fitting to impose the values of U in the control points at the top
      %boundary
      id_top = id_bnd{4};
      u_top  = BndFitting_lme(x_nodes(id_top,1),optLME);
      for j=1:length(id_top)
        jj      = id_top(j);
        x_bd    = x_nodes(jj,:);
        u_val   = u_top(j);
        f       = f - K(:,jj)*u_val;
        K(:,jj) = 0;
        K(jj,:) = 0;
        K(jj,jj)= 1;
        f(jj)   = u_val;
      end
      
      % -------------------------------------------------------------------
      % Compute System Solution
      u_lme = K\f;
      
      fprintf(1,'\tcputime System : %4.3f\n', cputime-t);
      
      %% Error in L2 norm
      t = cputime;
      errL2 = 0;
      for k =1:sPts
        k_near = s_near{k};
        u_k    = u_lme(k_near);
        p_k    = p_samp{k};
        
        u_h   = p_k'*u_k;
        u     = NoHom_AnalyticalSolution(1,x_samples(k,:));
        errL2 = errL2 + (u-u_h)^2 * w_samples(k);
      end
      normL2(n) = sqrt(errL2);
      fprintf(1,'\tcputime L2 norm: %4.3f\n', cputime-t);
    end
    
    disp('---------------------------------------------------');
    
    %% =======================================================================================
    gamma = gammaV(ig);
    Tol0  = optLME.Tol0;
    gPts  = gPtsV(gp);
    
    %% Save convergence results
    if saveData == 1
        fileOUT = sprintf('normL2/normL2_lme_g%3.1f_GP%02d.mat',gamma,gPts);
        save(fileOUT,'Ndof','normL2','gamma','Tol0','gPts');
    end
  end
end

return

%% Plot results
t = cputime;

n_id  = (n-1):n;
mrate = polyfit(log(1./(Ndof(n_id)-1)),log(normL2(n_id)'),1);

mref  = polyfit(log([0.1,0.01]),log([0.0005,0.000005]),1);
fprintf(1,'\tRate of convergence  m=%4.2f  mref=%4.2f\n', mrate(1), mref(1));

data_fem   = load('normL2/normL2_fem.mat');
Ndof_fem   = data_fem.Ndof;
normL2_fem = data_fem.normL2;

% Error in L2 norm
figure(2);clf
loglog(1./(Ndof-1),normL2,'ro-','markersize',10,'linewidth',2)
hold on
loglog(1./(Ndof_fem-1),normL2_fem,'bs-','markersize',10,'linewidth',2)
loglog([0.1,0.01],[0.0005,0.000005],'k-','linewidth',2)
hold off
title('LME vs FEM  Error in L2 norm  |u-u_h|_2','fontsize',18,'fontweight','b')
xlabel('DOF', 'fontsize',16,'fontweight','b')
ylabel('|u-u_h|_2', 'fontsize',16,'fontweight','b')
legend('L_2 LME','L_2 FEM')


% Analytical solution
u_sol = NoHom_AnalyticalSolution(1,x_nodes);

% Numerical solution
figure(3);clf
plot3(x_nodes(:,1),x_nodes(:,2),u_lme,'ro')
hold on
plot3(x_nodes(:,1),x_nodes(:,2),u_sol,'b*')
hold off

return

nSX = 101;
nSY = 101;
figure(4);clf
plotLME2D(x_nodes,v_lme(1:nPts), nSX, nSY, Lx, Ly, [0 0], optLME);
view([-17 4])
title(strcat('LME  \gamma = ',num2str(gamma)),'fontsize',18,'fontweight','b')
xlabel('X', 'fontsize',16,'fontweight','b')
ylabel('Y', 'fontsize',16,'fontweight','b')
zlabel('Temperature','fontsize',16,'fontweight','b')
fprintf(1,'cputime Plot: %4.3f\n', cputime-t);