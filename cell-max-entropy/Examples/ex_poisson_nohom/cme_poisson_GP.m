%% EXAMPLE: Heat Equation with a distributed source - Cell-based max-ent basis functions

clear all
close all

format short

addpath('ex_toolbox_local')

disp('------------------------------------------------')

%% =======================================================================================
%  Input parameters: CME options

optCME.dim   = 2;
optCME.verb  = 0;           % 0:off    1:on
optCME.grad  = 1;           % Computation of the Gradient 0:OFF 1:ON
optCME.spow  = 0;
optCME.mpow  = 0;
optCME.nring = 0;
optCME.TolNR = 1.e-12;      %Newton-Raphson tolerance fo the nonlinear max-ent problem
optCME.isConv= 1;       % flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)

% -------------------------------------------------------------------
% Integration options
%Numerical integration, cubature order
orderGL= [15]; %[1 2 3 4 5  6 10 15];
gPtsV  = [54]; %[1 3 4 6 7 12 25 54];
optGL.quality = 0.01;

for nring=2:3
  optCME.nring = nring;
  for spow=3:4
    optCME.spow = spow;
    
    for mpow=2:3
      optCME.mpow = mpow;
      for gp=1:length(orderGL)
        disp('=============================================================')
        %order 1 (1 gPts), 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
        optGL.orderGL = orderGL(gp);
        
        fprintf(1,'\tNumerical integration Gauss-Legendre  order=%2d  gPts=%2d\n',orderGL(gp), gPtsV(gp));
        
        % -------------------------------------------------------------------
        % Geometrical Parameters
        L     = 1.0;    % length
        
        if optCME.nring == 2
          Ndof  = [5 9 17 33 65 129];
        else
          Ndof  = [9 17 33 65 129];
        end
        normL2= zeros(length(Ndof),1);
        centre= L*[0.5,0.5];
        
        for n=1:length(Ndof)
          t=cputime;
          
          % Node Points
          N       = Ndof(n);
          h       = L/(N-1);
          
          fprintf(1,'\tGrid size: %2dx%2d    nodal spacing=%6.3g\n',N,N,h);
          dataDir= 'data_mesh';
          fileIN = sprintf('%s/square_mesh_N%03d_ring%d.mat', dataDir, N, optCME.nring);
          data   = load(fileIN);
          
          nPts    = length(data.x_n);
          
          x_nodes = data.x_n + ones(nPts,1)*centre;
          
          
          ids   = 1:nPts;
          id_bnd{1} = ids(abs(x_nodes(:,1)+0)<1.e-6); % (x=0)
          id_bnd{2} = ids(abs(x_nodes(:,1)-L)<1.e-6); % (x=L)
          id_bnd{3} = ids(abs(x_nodes(:,2)+0)<1.e-6); % (y=0)
          id_bnd{4} = ids(abs(x_nodes(:,2)-L)<1.e-6); % (y=L)
          
          tri       = data.tri;
          connecStr = data.connecStr;
          
          % -------------------------------------------------------------------
          % Gauss Legendre quadrature points for numerical integration
          [x_samples,w_samples] = MakeGLSamples2D(x_nodes, optGL);
          sPts  = size(x_samples,1);  %total number of Gauss points
          nElem = size(tri,1);        %number of elcments
          
          inside = ones(nElem,1);
          
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
          triplot(tri(inside==1, :),x_nodes(:,1),x_nodes(:,2),'m-','LineWidth',1);
          axis equal
          
          fprintf(1,'\tcputime CME Connectivity Structures: %5.3f seconds\n', cputime-t);
          
          
          %% =======================================================================================
          %  CME: basis functions computation
          
          t=cputime;
          
          optCME.tri           = tri;
          optCME.bd_segm       = connecStr.bd_segm;
          optCME.bd_segm_out   = connecStr.bd_segm_out;
          optCME.bd_nodes      = connecStr.bd_nodes;
          optCME.nodes_in_elem = connecStr.nodes_in_elem;
          
          outCME  = wrapper_cme(x_nodes,x_samples,optCME);
          
          p_samp  = outCME.p_samp;
          dp_samp = outCME.dp_samp;
          s_near  = outCME.s_near;
          
          fprintf(1,'\tcputime CME basis functions    BULK: %5.3f seconds\n', cputime-t);
          
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
              u_val   = NoHom_AnalyticalSolution(1,x_bd);
              f       = f - K(:,jj)*u_val;
              K(:,jj) = 0;
              K(jj,:) = 0;
              K(jj,jj)= 1;
              f(jj)   = u_val;
            end
          end
          
          % -----------------------------------------------------------------------
          %Least squares fitting to impose the values of U in the control points at
          %the top boundary
          id_top = id_bnd{4};
          u_top  = BndFitting_cme(x_nodes,id_top,optCME);
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
          u_cme = K\f;
          
          fprintf(1,'\tcputime System : %4.3f\n', cputime-t);
          
          %% Error in L2 norm
          t = cputime;
          errL2 = 0;
          for k =1:sPts
            k_near = s_near{k};
            u_k    = u_cme(k_near);
            p_k    = p_samp{k};
            
            u_h   = p_k'*u_k;
            u     = NoHom_AnalyticalSolution(1,x_samples(k,:));
            errL2 = errL2 + (u-u_h)^2 * w_samples(k);
          end
          normL2(n) = sqrt(errL2);
          fprintf(1,'\tcputime L2 norm: %4.3f\n', cputime-t);
          disp('-------------------');
        end
        
        disp('---------------------------------------------------');
        
        %% =======================================================================================
        gPts  = gPtsV(gp);
        
        fileOUT = sprintf('normL2/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts);
        %save(fileOUT,'Ndof','normL2','spow','mpow','nring','gPts');
        
      end
    end
  end
end


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
title('FEM  Error in L2 norm  |u-u_h|_2','fontsize',18,'fontweight','b')
xlabel('DOF', 'fontsize',16,'fontweight','b')
ylabel('|u-u_h|_2', 'fontsize',16,'fontweight','b')
legend('L_2 CME','L_2 FEM')


% Analytical solution
u_sol = NoHom_AnalyticalSolution(1,x_nodes);


% Numerical solution
figure(3);clf
plot3(x_nodes(:,1),x_nodes(:,2),u_cme,'ro')
hold on
plot3(x_nodes(:,1),x_nodes(:,2),u_sol,'b*')
hold off