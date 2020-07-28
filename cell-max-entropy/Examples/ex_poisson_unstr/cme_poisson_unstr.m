%% Timoshenko cantilever beam - cell-based max-entropy

clear all
close all

format short
addpath('ex_toolbox_local')

fprintf(1,'====================================================================\n');

%save data for this run -- check end lines of "Convervence curves"
saveData = 0; % 1: Yes     0: No

%% CME options
%  Input parameters: CME options
TolNR = 1.e-12;  % Newton-Raphson tolerance fo the nonlinear max-ent problem

optCME.verb  = 0;       % 0:off    1:on
optCME.grad  = 1;       % Computation of the Gradient 0:OFF 1:ON
optCME.isConv= 1;       % flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)

optCME.nring = 3;       % number ring of neighbors to consider for each node: 1, 2, 3.
optCME.spow  = 4;       % w^spow, where w is the approximation to the distance function (R-function)
optCME.mpow  = 2;       % distance function derivative maximum degree of the approximation

      
disp('###############################################################################')
fprintf(1,'\tCME  NRing=%d  spow=%d  mpow=%d\n', ...
    optCME.nring, optCME.spow, optCME.mpow);

%% Convervence curves
N_     = [81 120 283 581 1121 2132]; %4319 8510 16834];
m      = length(N_);

m_id   = (m-2:m);      %index of the L2 error norm to perform the fit
Ndof   = zeros(m,1);
normL2 = zeros(m,1);
normE  = zeros(m,1);
nnzK   = zeros(m,1);

% -------------------------------------------------------------------
% Integration options
% Numerical integration, cubature order
% cubature order: 1 (1 GPts), 2 (3 GPts), 4 (6 GPts), 5 (7 GPts), 6 (12 GPts)
orderGL= [6];
gPtsV  = [12];

for gp=1:length(orderGL)
    disp('=============================================================')
    %order 1 (1 gPts), 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
    optGL.orderGL = orderGL(gp);
    
    fprintf(1,'\tNumerical integration Gauss-Legendre  order=%2d  gPts=%2d\n',orderGL(gp), gPtsV(gp));
    
    for n=1:length(N_)
        disp('------------------------------------------------')
        
        if n>3
            optCME.TolNR = TolNR;
        else
            optCME.TolNR = 1.e-11;
        end
        % Node Points
        fileIN  = sprintf('data_mesh/mesh_N%05d_ring%d.mat', N_(n), optCME.nring);
        data    = load(fileIN);
        x_n     = data.x_n;
        tri     = data.tri;
        ibnd    = data.ibnd;
        connecStr = data.connecStr;
        
        nPts    = size(x_n,1);
        Ndof(n) = sqrt(nPts);
        
        % Sample points and gauss weights
        optGL.tri = tri;
        optGL.quality = 0.001;
        [x_s,w_s] = MakeGLSamples2D(x_n, optGL);
        sPts  = size(x_s,1);  %total number of Gauss points
        nElem = size(tri,1);        %number of elements
        gPts  = sPts/nElem;         %number of Gauss points in an element
        
        fprintf(1,'\tN=%5d  nring=%d | s=%d  m=%d | order=%2d\n',...
            N_(n), optCME.nring, optCME.spow, optCME.mpow, optGL.orderGL);
        
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
        
        
        %% =============================================================================
        % -------------------------------------------------------------------
        % System matrix ASSEMBLY
        t = cputime;
        
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
        %           ids        = 1:nPts;
        %           ind        = ids(abs(x_n(:,1))<1e-04 & abs(x_n(:,2))<1e-04);
        %           K(ind,:)   = 0;
        %           K(:,ind)   = 0;
        %           K(ind,ind) = 1;
        %
        %           rhs(ind)   = 0;
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
        
        % ------------------------------------------------------------------------
        %% L2 norm and energy norm
        [normL2(n),normE(n)] = ErrorNorms(u_num,optAssem);
        
        fprintf(1,'\tcputime system: %4.3f\n', cputime-t);
    end
    
    disp('---------------------------------------------------');
    
    %% Save data                 -------------------------------------------------------------
    if saveData == 1
        fileOUT = sprintf('normL2/normL2_cme_r%d_s%02d_m%02d_GP%02d.mat',nring,spow,mpow,gPts);
        fprintf(1,'Save output in file: %s\n',fileOUT);
        save(fileOUT,'Ndof','N_','normL2','nnzK','normE','spow','mpow','nring','gPts');
    end
end


%% Analysis of Results       -------------------------------------------------------------
data_fem   = load('normL2/normL2_fem_GP03.mat');
Ndof_fem   = data_fem.Ndof;
normL2_fem = data_fem.normL2;
normE_fem  = data_fem.normE;

fit_num1 = polyfit(log10(1./Ndof(m_id)),log10(normL2(m_id)),1);
fit_num2 = polyfit(log10(1./Ndof(m_id)),log10(normL2(m_id)),1);
disp('- - - - - - - - - - - - - - -');
fprintf ( 1, 'CME  : m=%7.4g  b=%7.4g\n', fit_num1(1),  fit_num1(2));
fprintf ( 1, '     : m=%7.4g  b=%7.4g\n', fit_num2(1),  fit_num2(2));


% Plots
h_line=[0.01 0.1]*2;
err_line=[0.0001 0.01]*2;
fit_line=polyfit(log10(h_line),log10(err_line),1);
fprintf ( 1, 'line : m=%7.4g  b=%7.4g\n', fit_line(1),fit_line(2));

h_line=[0.01 0.1]*2;
err_line1=[0.0050 0.05]*2;
err_line2=[0.0001 0.01]*2;
fit_line1=polyfit(log10(h_line),log10(err_line1),1);
fit_line2=polyfit(log10(h_line),log10(err_line2),1);
fprintf ( 1, 'line : mref1=%7.4g   mref2=%7.4g\n', fit_line1(1),fit_line2(1));


%Error in L2 norm
figure(1);clf
loglog(1./Ndof,normL2,'ro-','LineWidth',2,'MarkerSize',6)
hold on
loglog(1./Ndof_fem,normL2_fem,'bs-',h_line, err_line2,'k-','LineWidth',2,'MarkerSize',6)
hold off
xlabel('DOF^{-1/2}','fontsize',14,'fontweight','b')
ylabel('L2 Error','fontsize',14,'fontweight','b')
legend('CME','FEM','Location','Best','Orientation','horizontal')
title('Error in Norm L_2','fontsize',16,'fontweight','b')
set(gca,'XTick',0.02:0.02:0.1)
%xlim([0.001 0.11])

%Error in energy norm
figure(2);clf
loglog(1./Ndof,normE,'ro-','LineWidth',2,'MarkerSize',6)
hold on
loglog(1./Ndof_fem,normE_fem,'bs-',h_line, err_line1,'k-','LineWidth',2,'MarkerSize',6)
hold off
xlabel('DOF^{-1/2}','fontsize',14,'fontweight','b')
ylabel('Energy semi-norm','fontsize',14,'fontweight','b')
legend('CME','FEM','Location','Best','Orientation','horizontal')
title('Error in Energy semi-norm','fontsize',16,'fontweight','b')
set(gca,'XTick',0.02:0.02:0.1)

% Node and Sample Points
figure(3)
plot(x_n(:,1),x_n(:,2),'or',...
  x_s(:,1),x_s(:,2),'.b');
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