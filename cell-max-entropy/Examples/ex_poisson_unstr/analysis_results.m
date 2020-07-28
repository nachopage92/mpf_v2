%% EXAMPLE: Cantilever beam - cell-based max-ent basis functions

clear all
close all

format long

%number of Gauss points for CME
GP = 12;

disp('------------------------------------------------')

data_ex  = load('normL2/normL2_exact.mat');
normL2_ex= data_ex.normL2(end);
normE_ex = data_ex.normE(end);

%% Plot results -----------------------------------------------------

% FEM -----------------------------------------------------
data_fem1  = load('normL2/normL2_fem_GP01.mat');
Ndof_fem1  = data_fem1.Ndof;
normL2_fem1= data_fem1.normL2/normL2_ex;
normE_fem1 = data_fem1.normE/normE_ex;

data_fem2  = load('normL2/normL2_fem_GP03.mat');
Ndof_fem2  = data_fem2.Ndof;
normL2_fem2= data_fem2.normL2/normL2_ex;
normE_fem2 = data_fem2.normE/normE_ex;

data_fem3  = load('normL2/normL2_fem_GP06.mat');
Ndof_fem3  = data_fem3.Ndof;
normL2_fem3= data_fem3.normL2/normL2_ex;
normE_fem3 = data_fem3.normE/normE_ex;

data_fem4  = load('normL2/normL2_fem_GP12.mat');
Ndof_fem4  = data_fem4.Ndof;
normL2_fem4= data_fem4.normL2/normL2_ex;
normE_fem4 = data_fem4.normE/normE_ex;


% LME -----------------------------------------------------
data_lme   = load('normL2/normL2_lme_g0.8_GP12_Tol1e6.mat');
Ndof_lme1  = data_lme.Ndof;
normL2_lme1= data_lme.normL2/normL2_ex;
normE_lme1 = data_lme.normE/normE_ex;

data_lme   = load('normL2/normL2_lme_g1.8_GP07_Tol1e6.mat');
Ndof_lme2  = data_lme.Ndof;
normL2_lme2= data_lme.normL2/normL2_ex;
normE_lme2 = data_lme.normE/normE_ex;

data_lme   = load('normL2/normL2_lme_g4.8_GP06_Tol1e6.mat');
Ndof_lme3  = data_lme.Ndof;
normL2_lme3= data_lme.normL2/normL2_ex;
normE_lme3 = data_lme.normE/normE_ex;

data_lme   = load('normL2/normL2_lme_g1.8_GP12_Tol1e6.mat');
Ndof_lme4  = data_lme.Ndof;
normL2_lme4= data_lme.normL2/normL2_ex;
normE_lme4 = data_lme.normE/normE_ex;



% % CME -----------------------------------------------------
% --------------------- NRING 2 --------------------
data_cme    = load(sprintf('normL2/normL2_cme_r2_s03_m02_GP%02d.mat',GP));
Ndof_cme21  = data_cme.Ndof;
normL2_cme21= data_cme.normL2/normL2_ex;
normE_cme21 = data_cme.normE/normE_ex;

data_cme    = load(sprintf('normL2/normL2_cme_r2_s04_m02_GP%02d.mat',GP));
Ndof_cme22  = data_cme.Ndof;
normL2_cme22= data_cme.normL2/normL2_ex;
normE_cme22 = data_cme.normE/normE_ex;

data_cme    = load(sprintf('normL2/normL2_cme_r2_s03_m03_GP%02d.mat',GP));
Ndof_cme23  = data_cme.Ndof;
normL2_cme23= data_cme.normL2/normL2_ex;
normE_cme23 = data_cme.normE/normE_ex;

data_cme    = load(sprintf('normL2/normL2_cme_r2_s04_m03_GP%02d.mat',GP));
Ndof_cme24  = data_cme.Ndof;
normL2_cme24= data_cme.normL2/normL2_ex;
normE_cme24 = data_cme.normE/normE_ex;

% --------------------- NRING 3 --------------------
data_cme   = load(sprintf('normL2/normL2_cme_r3_s03_m02_GP%02d.mat',GP));
Ndof_cme31  = data_cme.Ndof;
normL2_cme31= data_cme.normL2/normL2_ex;
normE_cme31 = data_cme.normE/normE_ex;

data_cme   = load(sprintf('normL2/normL2_cme_r3_s04_m02_GP%02d.mat',GP));
Ndof_cme32  = data_cme.Ndof;
normL2_cme32= data_cme.normL2/normL2_ex;
normE_cme32 = data_cme.normE/normE_ex;

data_cme    = load(sprintf('normL2/normL2_cme_r3_s03_m03_GP%02d.mat',GP));
Ndof_cme33  = data_cme.Ndof;
normL2_cme33= data_cme.normL2/normL2_ex;
normE_cme33 = data_cme.normE/normE_ex;

data_cme    = load(sprintf('normL2/normL2_cme_r3_s04_m03_GP%02d.mat',GP));
Ndof_cme34  = data_cme.Ndof;
normL2_cme34= data_cme.normL2/normL2_ex;
normE_cme34 = data_cme.normE/normE_ex;

mref  = polyfit(log([0.05,0.005]),log([0.001,0.00001]),1);
fprintf(1,'\tRate of convergence  mref=%4.2f\n', mref(1));
mref  = polyfit(log([0.05,0.005]),log([0.001,0.0001]),1);
fprintf(1,'\tRate of convergence  mref=%4.2f\n', mref(1));


%% ##### Error in L2 norm
figure(1);clf
loglog(1./Ndof_fem3,normL2_fem3,'go-','markersize',12,'linewidth',2)
hold on
loglog(1./Ndof_cme21,normL2_cme21,'rs-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme22,normL2_cme22,'ro-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme23,normL2_cme23,'r^-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme24,normL2_cme24,'r<-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme31,normL2_cme31,'ms-.','markersize',12,'linewidth',2)
loglog(1./Ndof_cme32,normL2_cme32,'m*-.','markersize',12,'linewidth',2)
loglog(1./Ndof_cme33,normL2_cme33,'m^-.','markersize',12,'linewidth',2)
loglog(1./Ndof_cme34,normL2_cme34,'m<-.','markersize',12,'linewidth',2)
loglog(1./Ndof_lme1,normL2_lme1,'bo--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog(1./Ndof_lme2,normL2_lme2,'bs--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog(1./Ndof_lme3,normL2_lme3,'bv--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog([0.05,0.005],[0.4,0.004],'k-','linewidth',2)
loglog(1./Ndof_fem1,normL2_fem1,'g*-','markersize',12,'linewidth',2)
hold off
title('FEM - CME - LME - Error in L2 norm  |u-u_h|_2','fontsize',20,'fontweight','b')
xlabel('DOF^{-1/2}', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
legend(...
  ' FEM',...
  ' CME (N_R = 2, s = 3, m = 2)',...
  ' CME (N_R = 2, s = 4, m = 2)',...
  ' CME (N_R = 2, s = 3, m = 3)',...
  ' CME (N_R = 2, s = 4, m = 3)',...
  ' CME (N_R = 3, s = 3, m = 2)',...
  ' CME (N_R = 3, s = 4, m = 2)',...
  ' CME (N_R = 3, s = 3, m = 3)',...
  ' CME (N_R = 3, s = 4, m = 3)',...
  ' LME (\gamma = 0.8)',...
  ' LME (\gamma = 1.8)',...
  ' LME (\gamma = 4.8)',...
  ' slope=2')

figure(2);clf
loglog(1./Ndof_fem3,normE_fem3,'go-','markersize',12,'linewidth',2)
hold on
loglog(1./Ndof_cme21,normE_cme21,'rs-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme22,normE_cme22,'ro-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme23,normE_cme23,'r^-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme24,normE_cme24,'r<-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme31,normE_cme31,'ms-.','markersize',12,'linewidth',2)
loglog(1./Ndof_cme32,normE_cme32,'m*-.','markersize',12,'linewidth',2)
loglog(1./Ndof_cme33,normE_cme33,'m^-.','markersize',12,'linewidth',2)
loglog(1./Ndof_cme34,normE_cme34,'m<-.','markersize',12,'linewidth',2)
loglog(1./Ndof_lme1,normE_lme1,'bo--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog(1./Ndof_lme2,normE_lme2,'bs--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog(1./Ndof_lme3,normE_lme3,'bv--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog([0.05,0.005],[0.4,0.04],'k-','linewidth',2)
loglog(1./Ndof_fem1,normE_fem1,'g*-','markersize',12,'linewidth',2)
hold off
title('FEM - CME - LME - Error in Energy norm','fontsize',20,'fontweight','b')
xlabel('DOF^{-1/2}', 'fontsize',24,'fontweight','b')
ylabel('Energy error', 'fontsize',24,'fontweight','b')
legend(...
  ' FEM',...
  ' CME (N_R = 2, s = 3, m = 2)',...
  ' CME (N_R = 2, s = 4, m = 2)',...
  ' CME (N_R = 2, s = 3, m = 3)',...
  ' CME (N_R = 2, s = 4, m = 3)',...
  ' CME (N_R = 3, s = 3, m = 2)',...
  ' CME (N_R = 3, s = 4, m = 2)',...
  ' CME (N_R = 3, s = 3, m = 3)',...
  ' CME (N_R = 3, s = 4, m = 3)',...
  ' LME (\gamma = 0.8)',...
  ' LME (\gamma = 1.8)',...
  ' LME (\gamma = 4.8)',...
  ' slope=1')


%%
disp('=============')
m_id  = (length(Ndof_fem3)-1):length(Ndof_fem3);
m_id2 = (length(Ndof_cme22)-1):length(Ndof_cme22);
m_id3 = (length(Ndof_cme31)-1):length(Ndof_cme31);
fit_fem_L  = polyfit(log10(1./Ndof_fem3 (m_id )),log10(normL2_fem3 (m_id )),1);
fit_fem_E  = polyfit(log10(1./Ndof_fem3 (m_id )),log10(normE_fem3  (m_id )),1);
fit_lme_L  = polyfit(log10(1./Ndof_lme1 (m_id )),log10(normL2_lme1 (m_id )),1);
fit_lme_E  = polyfit(log10(1./Ndof_lme1 (m_id )),log10(normE_lme1  (m_id )),1);
fit_cme_L2 = polyfit(log10(1./Ndof_cme22(m_id2)),log10(normL2_cme22(m_id2)),1);
fit_cme_E2 = polyfit(log10(1./Ndof_cme22(m_id2)),log10(normE_cme22 (m_id2)),1);
fit_cme_L3 = polyfit(log10(1./Ndof_cme31(m_id3)),log10(normL2_cme31(m_id3)),1);
fit_cme_E3 = polyfit(log10(1./Ndof_cme31(m_id3)),log10(normE_cme31 (m_id3)),1);


disp('- - - - - - - - - - - - - - -');
h_line1=[0.012 0.08]; err_line1=[0.0000225 0.001];
h_line2=[0.006 0.04]; err_line2=[0.1 0.67];
fit_line1=polyfit(log10(h_line1),log10(err_line1),1);
fit_line2=polyfit(log10(h_line2),log10(err_line2),1);


fprintf ( 1, 'FEM  normL2: m=%7.3g\n', fit_fem_L (1));
fprintf ( 1, 'LME  normL2: m=%7.3g  gamma=0.8\n', fit_lme_L (1));
fprintf ( 1, 'CME  normL2: m=%7.3g  NRing=2\n', fit_cme_L2(1));
fprintf ( 1, 'CME  normL2: m=%7.3g  NRing=3\n', fit_cme_L3(1));
fprintf ( 1, 'ref slope1 : m=%7.3g\n', fit_line1(1));

fprintf ( 1, 'FEM  normE : m=%7.3g\n', fit_fem_E(1));
fprintf ( 1, 'LME  normE : m=%7.3g  gamma=0.8\n', fit_lme_E(1));
fprintf ( 1, 'CME  normE : m=%7.3g  NRing=2\n', fit_cme_E2(1));
fprintf ( 1, 'CME  normE : m=%7.3g  NRing=3\n', fit_cme_E3(1));
fprintf ( 1, 'ref slope2 : m=%7.3g\n', fit_line2(1));


figure(10);clf
subplot(1,2,1)
loglog(1./Ndof_fem3,normL2_fem3,'rd:','markersize',12,'linewidth',2)
hold on
loglog(1./Ndof_lme1,normL2_lme1,'bo--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_lme2,normL2_lme2,'bs--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_lme3,normL2_lme3,'bv--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_cme22,normL2_cme22,'ko-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme31,normL2_cme31,'ks-','markersize',12,'linewidth',2)
loglog(h_line1,err_line1,'k-','linewidth',2)
hold off
xlabel('DOF^{-1/2}', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
xlim([0.005,0.12])
ylim([1e-05,10])

subplot(1,2,2)
loglog(1./Ndof_fem3,normE_fem3,'rd:','markersize',12,'linewidth',2)
hold on
loglog(1./Ndof_lme1,normE_lme1,'bo--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_lme2,normE_lme2,'bs--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_lme3,normE_lme3,'bv--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_cme22,normE_cme22,'ko-','markersize',12,'linewidth',2)
loglog(1./Ndof_cme31,normE_cme31,'ks-','markersize',12,'linewidth',2)
loglog(h_line2,err_line2,'k-','linewidth',2)
hold off
xlabel('DOF^{-1/2}', 'fontsize',24,'fontweight','b')
ylabel('Energy error', 'fontsize',24,'fontweight','b')
legend(...
  ' FEM ',...
  ' LME (\gamma = 0.8)',...
  ' LME (\gamma = 1.8)',...
  ' LME (\gamma = 4.8)',...
  ' CME (N_R = 2)',...
  ' CME (N_R = 3)',...
  ' slope reference','Location','SouthEast')
xlim([0.005,0.12])
ylim([1e-05,10])
% Create textbox
annotation('textbox',...
  [0.269444444444445 0.724978647686831 0.0322916666666667 0.0640569395017794],'String',{'2'},...
  'FontWeight','bold',...
  'FontSize',20,...
  'LineStyle','none');
annotation('textbox',...
  [0.710972222222223 0.79750533807829 0.0322916666666667 0.0640569395017794],'String',{'1'},...
  'FontWeight','bold',...
  'FontSize',20,...
  'LineStyle','none');



%%

figure(20);clf
subplot(1,2,1)
loglog(1./Ndof_fem1,normL2_fem1,'r*:','markersize',12,'linewidth',2)
hold on
loglog(1./Ndof_fem2,normL2_fem2,'rd:','markersize',12,'linewidth',2)
loglog(1./Ndof_fem3,normL2_fem3,'rd:','markersize',12,'linewidth',2)
loglog(1./Ndof_fem4,normL2_fem4,'rd:','markersize',12,'linewidth',2)
loglog(1./Ndof_lme1,normL2_lme2,'bo--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_lme2,normL2_lme4,'bs--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(h_line1,err_line1,'k-','linewidth',2)
hold off
xlabel('DOF^{-1/2}', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
xlim([0.005,0.12])
ylim([1e-05,1000])

subplot(1,2,2)
loglog(1./Ndof_fem1,normE_fem1,'r*:','markersize',12,'linewidth',2)
hold on
loglog(1./Ndof_fem2,normE_fem2,'rd:','markersize',12,'linewidth',2)
loglog(1./Ndof_fem3,normE_fem3,'rd:','markersize',12,'linewidth',2)
loglog(1./Ndof_fem4,normE_fem4,'rd:','markersize',12,'linewidth',2)
loglog(1./Ndof_lme1,normE_lme2,'bo--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(1./Ndof_lme2,normE_lme4,'bs--','markersize',12,'linewidth',2,'markerfacecolor',[0.0431372560560703 0.517647087574005 0.780392169952393])
loglog(h_line2,err_line2,'k-','linewidth',2)
hold off
xlabel('DOF^{-1/2}', 'fontsize',24,'fontweight','b')
ylabel('Energy error', 'fontsize',24,'fontweight','b')
legend(...
  ' FEM (GP =  1)',...
  ' FEM (GP =  3)',...
  ' FEM (GP =  6)',...
  ' FEM (GP = 12)',...
  ' LME (\gamma = 1.8, GP =  6)',...
  ' LME (\gamma = 1.8, GP = 12)',...
  ' slope reference','Location','SouthEast')
xlim([0.005,0.12])
ylim([1e-05,1000])
% Create textbox
annotation('textbox',...
  [0.269444444444445 0.724978647686831 0.0322916666666667 0.0640569395017794],'String',{'2'},...
  'FontWeight','bold',...
  'FontSize',20,...
  'LineStyle','none');
annotation('textbox',...
  [0.710972222222223 0.79750533807829 0.0322916666666667 0.0640569395017794],'String',{'1'},...
  'FontWeight','bold',...
  'FontSize',20,...
  'LineStyle','none');

