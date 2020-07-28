%% Timoshenko cantilever beam - local max-entropy

clear all
close all

format short
addpath('ex_toolbox_local')

fprintf(1,'====================================================================\n');

%save data for this run -- check end lines after "Convervence curves"
saveData = 0; % 1: Yes     0: No

%% LME options
optLME.dim   = 2;
optLME.verb  = 0; %0:off
optLME.grad  = 1;           % Computation of the Gradient 0:OFF 1:ON
optLME.hess  = 0;           % Computation of the Hessian  0:OFF 1:ON
optLME.TolNR = 1.e-12;
optLME.knn   = 0;

%% Convervence curves
N_     = [81 120 283 581 1121 2132];% 4319 8510 16834];
m      = length(N_);

m_id   = (m-2:m);      %index of the L2 error norm to perform the fit
h_n    = zeros(m,1);
Ndof   = zeros(m,1);

normL2 = zeros(m,1);
normE  = zeros(m,1);
nnzK   = zeros(m,1);

tpow  = 8;%[4,6,8];
Tol0V = 1e-8;%[1e-4,1e-6,1e-8];
gammaV= [1.8];
for gg=1:length(gammaV)
  optLME.gamma = gammaV(gg);      % value of gamma to compute LME
  
  % Nodal Integration options
  %cubature order: 1 (1 GPts), 2 (3 GPts), 4 (6 GPts), 5 (7 GPts), 6 (12 GPts)
  if optLME.gamma == 0.8
    optGL.orderGL = 6;
    gPtsV         = 12;
  elseif optLME.gamma == 1.8
    optGL.orderGL = 5;
    gPtsV         = 7;
  else
    optGL.orderGL = 2;
    gPtsV         = 3;
  end
  
  
  disp('###########################################################################')
  fprintf(1,'Aspect shape parameter gamma=%4.2f\n', gammaV(gg));
  for it=1:length(Tol0V)
    optLME.Tol0  = Tol0V(it);
    disp('========================================================')
    fprintf(1,'Tolerance %8.2e\n', Tol0V(it));
    
    for n=1:length(N_)
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
      optGL.tri     = tri;
      optGL.quality = 0.001;
      [x_s,w_s] = MakeGLSamples2D(x_n, optGL);
      sPts  = size(x_s,1);  %total number of Gauss points
      nElem = size(tri,1);        %number of elements
      gPts  = sPts/nElem;         %number of Gauss points in an element
      
      %% Nodes thermalization, nodes and samples adjacency structures are computed
      optNS.dim   = 2;
      optNS.gamma = 1.5;
      optNS.knn   = 8;
      h_n         = NodalSpacing(x_n,optNS,ibnd);
      
      optLME.h_n    = h_n;
      optLME.knn    = 1;
      [beta_n,range_n] = NodalThermalization(x_n, optLME);
      optLME.beta_n  = beta_n;
      optLME.range_n = range_n;
      
      disp([min(range_n)/min(h_n) mean(range_n)/mean(h_n) max(range_n)/max(h_n)])
      
      %% The basis functions are computed
      t=cputime;
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
      
      %       for e=1:nElem
      %         ik = ((e-1)*gPtsV+1):e*gPtsV;
      %         figure(3);clf
      %         plot(x_n(:,1),x_n(:,2),'or',x_s(ik,1),x_s(ik,2),'.b');
      %         hold on
      %         plot(x_n(e_nears{e},1),x_n(e_nears{e},2),'ks','markersize',14,'markerfacecolor','g');
      %         plot(x_n(tri(e,:),1),x_n(tri(e,:),2),'md','markerfacecolor','c','markersize',14,'linewidth',2)
      %         hold off
      %         axis equal
      %         pause
      %       end
      
      for e=1:nElem
        for g=1:gPtsV
          k = (e-1)*gPtsV+g;
          if isempty(p_samp{k})==1
            figure(3);clf
            plot(x_n(:,1),x_n(:,2),'or',x_s(:,1),x_s(:,2),'.b');
            hold on
            plot(x_n(e_nears{e},1),x_n(e_nears{e},2),'ks','markersize',14,'markerfacecolor','g');
            plot(x_s(k,1),x_s(k,2),'md','markerfacecolor','c','markersize',14)
            hold off
            axis equal
            pause
          end
        end
      end
      
      
      %% =================================================================================
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
%       %(x=0,y=0)  u=0
%       ids        = 1:nPts;
%       ind        = ids(abs(x_n(:,1))<1e-04 & abs(x_n(:,2))<1e-04);
%       K(ind,:)   = 0;
%       K(:,ind)   = 0;
%       K(ind,ind) = 1;
%       
%       rhs(ind)   = 0;

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
    
    disp('---------------------------------------------------');
    
    %% Save data
    if saveData == 1
        gamma = optLME.gamma;
        Tol0  = optLME.Tol0;
        
        fileOUT = sprintf('normL2/normL2_lme_g%3.1f_GP%02d_Tol1e%d.mat',gamma,gPts,tpow(it));
        save(fileOUT,'Ndof','N_','normL2','normE','nnzK','gamma','Tol0','gPts');
    end
  end
end

%% Analysis of Results
data_fem   = load('normL2/normL2_fem_GP03.mat');
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
legend('LME','FEM','Location','Best','Orientation','horizontal')
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
legend('LME','FEM','Location','Best','Orientation','horizontal')
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
