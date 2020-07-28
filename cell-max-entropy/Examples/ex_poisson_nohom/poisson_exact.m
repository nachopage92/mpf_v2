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

%% Geometrical Parameters
L     = 1.0;    % length
Ndof  = [3 5 9 17 33 65 129 257 513 1025 2049];

normL2= zeros(length(Ndof),1);
centre= L*[0.5,0.5];

for n=1:length(Ndof)
  disp('------------------------------------------------')
  
  % Node Points
  N     = Ndof(n);
  x_nodes = UniformGrid2D(N, N, L, L, centre);
  
  fprintf(1,'\tGrid size: %2dx%2d\n',N,N);
  
  % GAUSS_LEGENDRE
  [x_samples,w_samples] = MakeGLSamples2D(x_nodes, optGL);

  fprintf(1,'\torderGL=%2d\n', optGL.orderGL);

  
  %% L2 norm
  t = cputime;
  u     = AnalyticalSolution(1,x_samples);
  errL2 = sum(u.^2 .* w_samples);
  normL2(n) = sqrt(errL2);
  fprintf(1,'\tL2 norm: %15.13g  computed in  %6.3f\n', normL2(n), cputime-t);
end

save('normL2/normL2_exact.mat','Ndof','normL2');

%% plot results
figure(1);clf
loglog(1./(Ndof-1),abs(normL2-normL2(end)),'ro-')
title('L2 norm  |u|_2','fontsize',20,'fontweight','b')
xlabel('Nodal spacing', 'fontsize',24,'fontweight','b')
ylabel('|u|_2', 'fontsize',24,'fontweight','b')

disp('---------------------------------------------------');
