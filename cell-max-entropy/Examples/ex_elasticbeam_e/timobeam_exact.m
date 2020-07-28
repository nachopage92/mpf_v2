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

%% Matrial Parameters
D  = 1;         % height
L  = 4*D;       % length

P  = -1000;     % external load
E  =  1e07;     % Young modulus
nu =   0.499;   % Poisson coefficient: 0.3, 0.499 (quasi incomp)

% sigma = C_stiff*epsilon
C_stiff = E/(1+nu)/(1-2*nu)*[1-nu,   nu,          0 ;...
                               nu, 1-nu,          0 ;...
                                0,    0, (1-2*nu)/2];
                              
%epsilon = S_stiff*sigma
S_stiff = (1+nu)/E*[1-nu,  -nu, 0;
                     -nu, 1-nu, 0;
                       0,    0, 2];

parameters.L  = L;
parameters.D  = D;
parameters.E  = E;
parameters.nu = nu;
parameters.P  = P;
Lx = L;
Ly = 0.5 * D;

%% Geometrical Parameters
%% Convervence curves
N_     = [3 5 9 17 33 65];

m      = length(N_);
m_id   = (m-2:m);      %index of the L2 error norm to perform the fit
h_n    = zeros(m,1);
Ndof   = zeros(m,1);
normL2 = zeros(m,1);
normE  = zeros(m,1);

centre= 0.5*[Lx,Ly];

for n=1:length(N_)
  disp('------------------------------------------------')
  % Node Points
  Ny = N_(n);
  Nx = N_(n)*8-3;
  hx = Lx/(Nx-1);
  hy = Ly/(Ny-1);
  fprintf(1,'Quasi Uniform grid %dx%d ---- h=[%5.3f %5.3f]\n', Nx, Ny, hx, hy);

  % Node Points
  x_nodes = UniformGrid2D(Nx, Ny, Lx, Ly, centre);
  nPts    = size(x_nodes,1);
  Ndof(n) = sqrt(nPts*2);
  
  
  % GAUSS_LEGENDRE
  [x_samples,w_samples] = MakeGLSamples2D(x_nodes, optGL);
  sPts = length(w_samples);
  fprintf(1,'\torderGL=%2d\n', optGL.orderGL);

  
  %% L2 norm
  t = cputime;
  
  %Displacement and stress tensor exact solutions for the Timoshenko's cantilever beam
  [u_s,sgm_s] = Beam_AnalyticalSolution(x_samples,P,D,L,E,nu);
  
  errL2 = 0;
  errE  = 0;
  for k =1:sPts
    %L2 nor of the error
    errL2 = errL2 + sum(u_s(k,:).^2) * w_samples(k);
    
    %Energy semi-norm   
    sgm_  = sgm_s(k,:)';
    
    eps_  = S_stiff*sgm_;      %epsilon = S_stiff*sigma
    
    errE  = errE  + eps_'*sgm_ * w_samples(k);
  end
  normL2(n) = sqrt(errL2);
  normE(n)  = sqrt(0.5*errE);
  
  fprintf(1,'\tL2 norm: %15.13g    Energy Error:%15.13g  computed  in  %6.3f\n', normL2(n), normE(n), cputime-t);
end

if nu==0.3
  save('normL2_nuA/normL2_exact.mat','Ndof','N_','normL2','normE');
else
  save('normL2_nuB/normL2_exact.mat','Ndof','N_','normL2','normE');
end

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
