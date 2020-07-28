%% Timoshenko cantilever beam - Finite Element Method (P1: linear interpolants)
% plane strain state

clear all
close all

format short
addpath('ex_toolbox_local')


fprintf(1,'====================================================================\n');

%% Convervence curves
N_ = [81 289 581 1121 2132]; % 4319 8510 16834];


for n=length(N_)-2:length(N_)
  disp('------------------------------------------------')
  % Node Points
  fprintf(1,'Unstructured set of points N=%d\n', N_(n));
  
  
  fileIN  = sprintf('data_mesh/mesh_N%05d_ring2.mat', N_(n));
  data    = load(fileIN);
  x_n     = data.x_n;
  tri     = data.tri;
  ibnd    = data.ibnd;
  
  % Deformed and undeformed (numerical and analytical) configurations
  [~,~,u_sol] = AnalyticalSolution(x_n);
  
  figure(1)
  plot3(x_n(:,1),x_n(:,2),u_sol,'.r')
  hold on
  trimesh(tri,x_n(:,1),x_n(:,2),u_sol)
  hold off
  xlabel('X')
  ylabel('Y')
%   axis equal
  pause
end