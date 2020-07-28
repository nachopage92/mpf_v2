%% EXAMPLE: Heat Equation with a distributed source - Finite Element Method with P1

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
optGL.orderGL = 6;  %order 2(3 gPts), 3(4), 4(6), 5(7), 6(12), 10(25)
optGL.quality = 0.01;


%% Convervence curves
N_ = [81 289 581 1121 2132 4319 8510 16834];

m      = length(N_);
m_id   = (m-2:m);      %index of the L2 error norm to perform the fit
h_n    = zeros(m,1);
Ndof   = zeros(m,1);
normL2 = zeros(m,1);
normE  = zeros(m,1);

for n=1:length(N_)
  disp('------------------------------------------------')
  % Node Points
  fileIN  = sprintf('data_mesh/mesh_N%05d_ring2.mat', N_(n));
  data    = load(fileIN);
  x_n     = data.x_n;
  tri     = data.tri;
  ibnd    = data.ibnd;
  
  nPts    = size(x_n,1);
  
  Ndof(n) = sqrt(nPts);
  
  % GAUSS_LEGENDRE
  [x_s,w_s] = MakeGLSamples2D(x_n, optGL);
  sPts = length(w_s);
  fprintf(1,'\torderGL=%2d\n', optGL.orderGL);

  
  %% L2 norm
  t = cputime;
  
  %Displacement and stress tensor exact solutions
  [u_s,du_s] = AnalyticalSolution(x_s);
  
  errL2 = 0;
  errE  = 0;
  for k =1:sPts
    %L2 nor of the error
    errL2 = errL2 + sum(u_s(k).^2) * w_s(k);
    
    %Energy semi-norm
    errE  = errE + (du_s(k,:)*du_s(k,:)') * w_s(k);
  end
  normL2(n) = sqrt(errL2);
  normE(n)  = sqrt(0.5*errE);
  
  fprintf(1,'\tL2 norm: %15.13g    Energy Error:%15.13g  computed  in  %6.3f\n', normL2(n), normE(n), cputime-t);
end
save('normL2/normL2_exact.mat','Ndof','N_','normL2','normE');

%% plot results
figure(1);clf
loglog(1./(Ndof),abs(normL2-normL2(end)),'ro-')
title('L2 norm  |u|_2','fontsize',20,'fontweight','b')
xlabel('Nodal spacing', 'fontsize',24,'fontweight','b')
ylabel('|u|_2', 'fontsize',24,'fontweight','b')

figure(2);clf
loglog(1./(Ndof),abs(normE-normE(end)),'ro-')
xlabel('Nodal spacing', 'fontsize',24,'fontweight','b')
ylabel('Energy norm', 'fontsize',24,'fontweight','b')
disp('---------------------------------------------------');
