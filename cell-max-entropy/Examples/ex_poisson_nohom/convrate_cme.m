%% EXAMPLE: Heat Equation with a distributed source - cell-based max-ent basis functions

clear all
close all

format short

disp('------------------------------------------------')
data_ex  = load('normL2/normL2_exact.mat');
normL2_ex= data_ex.normL2(end);


%% =======================================================================================
% Plot results

mref  = polyfit(log([0.1,0.01]),log([0.0001,0.000001]),1);

%% FME ---------------------------------------------------------------
disp('FME -------------------')
data_fem  = load('normL2/normL2_fem_GP03.mat');
Ndof_fem  = data_fem.Ndof;
errL2_fem = data_fem.normL2/normL2_ex;
n         = length(Ndof_fem);
n_id      = (n-1):n;
mrate     = polyfit(log(1./(Ndof_fem(n_id)-1)),log(errL2_fem(n_id)'),1);
fprintf(1,'\tRate of convergence  FEM  m=%4.2f  mref=%4.2f\n', mrate(1), mref(1));


%% CME ---------------------------------------------------------------
disp('CME -- ring 2 ---------')
% ring 2 --------------------
data_cme   = load('normL2/normL2_cme_r2_s03_m02_GP03.mat');
Ndof_cme21 = data_cme.Ndof;
errL2_cme21= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme21);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme21(n_id)-1)),log(errL2_cme21(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts= 3\n', mrate(1));

data_cme   = load('normL2/normL2_cme_r2_s03_m02_GP06.mat');
Ndof_cme22 = data_cme.Ndof;
errL2_cme22= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme22);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme22(n_id)-1)),log(errL2_cme22(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts= 6\n', mrate(1));
 
data_cme   = load('normL2/normL2_cme_r2_s03_m02_GP12.mat');
Ndof_cme23 = data_cme.Ndof;
errL2_cme23= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme23);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme23(n_id)-1)),log(errL2_cme23(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts=12\n', mrate(1));

data_cme   = load('normL2/normL2_cme_r2_s03_m02_GP25.mat');
Ndof_cme24 = data_cme.Ndof;
errL2_cme24= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme24);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme24(n_id)-1)),log(errL2_cme24(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts=25\n', mrate(1));

disp('CME -- ring 3 ---------')
% ring 3 --------------------
data_cme   = load('normL2/normL2_cme_r3_s03_m02_GP03.mat');
Ndof_cme31 = data_cme.Ndof;
errL2_cme31= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme31);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme31(n_id)-1)),log(errL2_cme31(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts= 3\n', mrate(1));

data_cme   = load('normL2/normL2_cme_r3_s03_m02_GP06.mat');
Ndof_cme32 = data_cme.Ndof;
errL2_cme32= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme32);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme32(n_id)-1)),log(errL2_cme32(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts= 6\n', mrate(1));
 
data_cme   = load('normL2/normL2_cme_r3_s03_m02_GP12.mat');
Ndof_cme33 = data_cme.Ndof;
errL2_cme33= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme33);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme33(n_id)-1)),log(errL2_cme33(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts=12\n', mrate(1));

data_cme   = load('normL2/normL2_cme_r3_s03_m02_GP25.mat');
Ndof_cme34 = data_cme.Ndof;
errL2_cme34= data_cme.normL2/normL2_ex;
n          = length(Ndof_cme34);
n_id       = (n-1):n;
mrate      = polyfit(log(1./(Ndof_cme34(n_id)-1)),log(errL2_cme34(n_id)'),1);
fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts=25\n', mrate(1));

%% LME ---------------------------------------------------------------
disp('LME -- \gamma=0.8,1.8,4.8 -----')
data_lme  = load('normL2/normL2_lme_g0.8_GP12.mat');
Ndof_lme1 = data_lme.Ndof(2:end);
errL2_lme1= data_lme.normL2(2:end)/normL2_ex;
n         = length(Ndof_lme1);
n_id      = (n-1):n;
mrate     = polyfit(log(1./(Ndof_lme1(n_id)-1)),log(errL2_lme1(n_id)'),1);
fprintf(1,'\tRate of convergence  LME  m=%4.2f  gPts=12\n', mrate(1));

data_lme  = load('normL2/normL2_lme_g1.8_GP06.mat');
Ndof_lme2 = data_lme.Ndof(2:end);
errL2_lme2= data_lme.normL2(2:end)/normL2_ex;
n         = length(Ndof_lme2);
n_id      = (n-1):n;
mrate     = polyfit(log(1./(Ndof_lme2(n_id)-1)),log(errL2_lme2(n_id)'),1);
fprintf(1,'\tRate of convergence  LME  m=%4.2f  gPts= 6\n', mrate(1));


data_lme  = load('normL2/normL2_lme_g4.8_GP03.mat');
Ndof_lme3 = data_lme.Ndof(2:end);
errL2_lme3= data_lme.normL2(2:end)/normL2_ex;
n         = length(Ndof_lme3);
n_id      = (n-1):n;
mrate     = polyfit(log(1./(Ndof_lme3(n_id)-1)),log(errL2_lme3(n_id)'),1);
fprintf(1,'\tRate of convergence  LME  m=%4.2f  gPts= 3\n', mrate(1));

% Error in L2 norm
figure(1);clf
loglog(1./(Ndof_fem-1),errL2_fem,'go-','markersize',12,'linewidth',2)
hold on
loglog(1./(Ndof_cme21-1),errL2_cme21,'rs-','markersize',12,'linewidth',2)
loglog(1./(Ndof_cme22-1),errL2_cme22,'r*-','markersize',12,'linewidth',2)
loglog(1./(Ndof_cme23-1),errL2_cme23,'r^-','markersize',12,'linewidth',2)
loglog(1./(Ndof_cme24-1),errL2_cme24,'r<-','markersize',12,'linewidth',2)
loglog(1./(Ndof_cme31-1),errL2_cme31,'ms-.','markersize',12,'linewidth',2)
loglog(1./(Ndof_cme32-1),errL2_cme32,'m*-.','markersize',12,'linewidth',2)
loglog(1./(Ndof_cme33-1),errL2_cme33,'m^-.','markersize',12,'linewidth',2)
loglog(1./(Ndof_cme34-1),errL2_cme34,'m<-.','markersize',12,'linewidth',2)
loglog(1./(Ndof_lme1-1),errL2_lme1,'bo--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog(1./(Ndof_lme2-1),errL2_lme2,'bs--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog(1./(Ndof_lme3-1),errL2_lme3,'bv--','markersize',12,'linewidth',2,'markerfacecolor','c')
loglog([0.1,0.01],[0.005,0.00005],'k-','linewidth',2)
hold off
title('FEM - CME - LME - Error in L2 norm  |u-u_h|_2','fontsize',20,'fontweight','b')
xlabel('nodal spacing', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
legend('FEM',...
  'CME  GP03  NRing=2  s=3  m=2',...
  'CME  GP06  NRing=2  s=3  m=2',...
  'CME  GP12  NRing=2  s=3  m=2',...
  'CME  GP25  NRing=2  s=3  m=2',...
  'CME  GP03  NRing=3  s=3  m=2',...
  'CME  GP06  NRing=3  s=3  m=2',...
  'CME  GP12  NRing=3  s=3  m=2',...
  'CME  GP25  NRing=3  s=3  m=2',...
  'LME  GP12  \gamma=0.8',...
  'LME  GP06  \gamma=1.8',...
  'LME  GP03  \gamma=4.8')


%% =======================================================================================
%% CME -- ring 2 -------------------------------------------------------------
disp('Full Gauss Points checking  - ring 2 --------------------')
GP = [1 3 4 6 7 12 25];
Ndof_cme2  = {[]};
errL2_cme2 = {[]};
L2_cme2    = zeros(1,length(GP));
for gp=1:length(GP)
  data_cme      = load(sprintf('normL2/normL2_cme_r2_s03_m02_GP%02d.mat',GP(gp)));
  Ndof_cme2{gp} = data_cme.Ndof;
  errL2_cme2{gp}= data_cme.normL2/normL2_ex;
  
  L2_cme2(gp) = data_cme.normL2(end)/normL2_ex;
  
  n     = length(Ndof_cme2{gp});
  n_id  = (n-1):n;
  mrate = polyfit(log(1./(Ndof_cme2{gp}(n_id)-1)),log(errL2_cme2{gp}(n_id)'),1);
  fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts=%2d\n', mrate(1), GP(gp));
end

% Error in L2 norm
figure(2);clf
loglog(1./(Ndof_fem-1),errL2_fem,'go-','markersize',10,'linewidth',2)
hold on
loglog(1./(Ndof_cme2{1}-1),errL2_cme2{1},'rs-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme2{2}-1),errL2_cme2{2},'r*-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme2{3}-1),errL2_cme2{3},'r^-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme2{4}-1),errL2_cme2{4},'r<-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme2{5}-1),errL2_cme2{5},'r>-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme2{6}-1),errL2_cme2{6},'rx-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme2{7}-1),errL2_cme2{7},'r+-','markersize',10,'linewidth',2)
loglog([0.1,0.01],[0.0001,0.000001],'k-','linewidth',2)
hold off
title('CME  Error in L2 norm  |u-u_h|_2','fontsize',20,'fontweight','b')
xlabel('nodal spacing', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
legend('L_2 FEM  GP03',...
  'L_2  CME  GP01 | ring=2  s=3  m=2',...
  'L_2  CME  GP03 | ring=2  s=3  m=2',...
  'L_2  CME  GP04 | ring=2  s=3  m=2',...
  'L_2  CME  GP06 | ring=2  s=3  m=2',...
  'L_2  CME  GP07 | ring=2  s=3  m=2',...
  'L_2  CME  GP12 | ring=2  s=3  m=2',...
  'L_2  CME  GP25 | ring=2  s=3  m=2')





%% CME -- ring 3 -------------------------------------------------------------
disp('Full Gauss Points checking  - ring 3 --------------------')
Ndof_cme3  = {[]};
errL2_cme3 = {[]};
L2_cme3    = zeros(1,length(GP));
for gp=1:length(GP)
  data_cme      = load(sprintf('normL2/normL2_cme_r3_s03_m02_GP%02d.mat',GP(gp)));
  Ndof_cme3{gp} = data_cme.Ndof;
  errL2_cme3{gp}= data_cme.normL2/normL2_ex;
  
  L2_cme3(gp) = data_cme.normL2(end)/normL2_ex;
  
  n     = length(Ndof_cme3{gp});
  n_id  = (n-1):n;
  mrate = polyfit(log(1./(Ndof_cme3{gp}(n_id)-1)),log(errL2_cme3{gp}(n_id)'),1);
  fprintf(1,'\tRate of convergence  CME  m=%4.2f  gPts=%2d\n', mrate(1), GP(gp));
end

% Error in L2 norm
figure(3);clf
loglog(1./(Ndof_fem-1),errL2_fem,'go-','markersize',10,'linewidth',2)
hold on
loglog(1./(Ndof_cme3{1}-1),errL2_cme3{1},'bs-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme3{2}-1),errL2_cme3{2},'b*-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme3{3}-1),errL2_cme3{3},'b^-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme3{4}-1),errL2_cme3{4},'b<-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme3{5}-1),errL2_cme3{5},'b>-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme3{6}-1),errL2_cme3{6},'bx-','markersize',10,'linewidth',2)
loglog(1./(Ndof_cme3{7}-1),errL2_cme3{7},'b+-','markersize',10,'linewidth',2)
loglog([0.1,0.01],[0.0001,0.000001],'k-','linewidth',2)
hold off
title('CME  Error in L2 norm  |u-u_h|_2','fontsize',20,'fontweight','b')
xlabel('nodal spacing', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
legend('L_2 FEM  GP03',...
  'L_2  CME  GP01 | ring=3  s=3  m=2',...
  'L_2  CME  GP03 | ring=3  s=3  m=2',...
  'L_2  CME  GP04 | ring=3  s=3  m=2',...
  'L_2  CME  GP06 | ring=3  s=3  m=2',...
  'L_2  CME  GP07 | ring=3  s=3  m=2',...
  'L_2  CME  GP12 | ring=3  s=3  m=2',...
  'L_2  CME  GP25 | ring=3  s=3  m=2')