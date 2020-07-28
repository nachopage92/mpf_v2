%% EXAMPLE: Heat Equation with a distributed source - Elcment max-ent basis functions

clear all
close all

format short

disp('------------------------------------------------')



%% =======================================================================================
%% CME -------------------------------------------------------------
disp('Full Gauss Points checking  - ring 2 --------------------')
%GP = [1 3 4 6 7 12 16 25];
GP = [3 4 6 7 12 13 16 25 54];
errL2_cme21 = zeros(length(GP),6);
errL2_cme22 = zeros(length(GP),6);
errL2_cme23 = zeros(length(GP),6);
errL2_cme24 = zeros(length(GP),6);

errL2_cme31 = zeros(length(GP),5);
errL2_cme32 = zeros(length(GP),5);
errL2_cme33 = zeros(length(GP),5);
errL2_cme34 = zeros(length(GP),5);


data_ex  = load('normL2/normL2_exact.mat');
normL2_ex= data_ex.normL2(end);

id2 = 5;
id3 = 4;

for gp=1:length(GP)
  % CME -- ring 2 -------------------------------------------------------------
  data_cme          = load(sprintf('normL2/normL2_cme_r2_s03_m02_GP%02d.mat',GP(gp)));
  errL2_cme21(gp,:) = data_cme.normL2/normL2_ex;

  data_cme          = load(sprintf('normL2/normL2_cme_r2_s04_m02_GP%02d.mat',GP(gp)));
  errL2_cme22(gp,:) = data_cme.normL2/normL2_ex;
% 
  data_cme          = load(sprintf('normL2/normL2_cme_r2_s03_m03_GP%02d.mat',GP(gp)));
  errL2_cme23(gp,:) = data_cme.normL2/normL2_ex;
% 
  data_cme          = load(sprintf('normL2/normL2_cme_r2_s04_m03_GP%02d.mat',GP(gp)));
  errL2_cme24(gp,:) = data_cme.normL2/normL2_ex;
end
N2 = data_cme.Ndof(id2);

for gp=1:length(GP)
  % CME -- ring 3 -------------------------------------------------------------
  data_cme          = load(sprintf('normL2/normL2_cme_r3_s03_m02_GP%02d.mat',GP(gp)));
  errL2_cme31(gp,:) = data_cme.normL2/normL2_ex;

  data_cme          = load(sprintf('normL2/normL2_cme_r3_s04_m02_GP%02d.mat',GP(gp)));
  errL2_cme32(gp,:) = data_cme.normL2/normL2_ex;

  data_cme          = load(sprintf('normL2/normL2_cme_r3_s03_m03_GP%02d.mat',GP(gp)));
  errL2_cme33(gp,:) = data_cme.normL2/normL2_ex;

  data_cme          = load(sprintf('normL2/normL2_cme_r3_s04_m03_GP%02d.mat',GP(gp)));
  errL2_cme34(gp,:) = data_cme.normL2/normL2_ex;
end
N3 = data_cme.Ndof(id3);

% Error in L2 norm
figure(1);clf
semilogy(GP,errL2_cme21(:,id2),'rs-','markersize',12,'linewidth',2)
hold on
semilogy(GP,errL2_cme22(:,id2),'r*-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme23(:,id2),'r^-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme24(:,id2),'r<-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme31(:,id3),'ms-.','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme32(:,id3),'m*-.','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme33(:,id3),'m^-.','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme34(:,id3),'m<-.','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme22(:,end),'r*-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme22(:,end),'r*-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme23(:,end),'r^-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme24(:,end),'r<-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme31(:,end),'ms-.','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme32(:,end),'m*-.','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme33(:,end),'m^-.','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme34(:,end),'m<-.','markersize',12,'linewidth',2)
hold off
title(sprintf('CME - Error in L2 norm  |u-u_h|_2 - %dx%d',N2,N3),'fontsize',20,'fontweight','b')
xlabel('number of Gauss points', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
legend(...
  'CME NRing=2  s=3  m=2',...
  'CME NRing=2  s=4  m=2',...
  'CME NRing=2  s=3  m=3',...
  'CME NRing=2  s=4  m=3',...
  'CME NRing=3  s=3  m=2',...
  'CME NRing=3  s=4  m=2',...
  'CME NRing=3  s=3  m=3',...
  'CME NRing=3  s=4  m=3')



%% Error in L2 norm
figure(2);clf
semilogy(GP,errL2_cme22(:,id2),'ko--','markersize',12,'linewidth',2)
hold on
semilogy(GP,errL2_cme31(:,id3),'ks--','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme22(:,end),'ko-','markersize',12,'linewidth',2)
semilogy(GP,errL2_cme32(:,end),'ks-','markersize',12,'linewidth',2)
hold off
xlabel('number of Gauss points', 'fontsize',24,'fontweight','b')
ylabel('L_2 error', 'fontsize',24,'fontweight','b')
legend(...
  'CME NRing=2 (s=4,m=2) h=L/64',...
  'CME NRing=3 (s=3,m=2) h=L/64',...
  'CME NRing=2 (s=4,m=2) h=L/128',...
  'CME NRing=3 (s=3,m=2) h=L/128')