%%
% This function plots the prior functions based in R-functions
% for a given node set, as well as the gradients and the approximation
% for discontinuous nodal data
%

clear all
close all
format long

% ========================================================================================
%% Input parameters

%  Input parameters: CME options
optCME.verb  = 0;      % 0:off    1:on
optCME.grad  = 1;      % Computation of the Gradient 0:OFF 1:ON
optCME.spow  = 1;      % w^spow, where w is the approximation to the distance function (R-function)
optCME.mpow  = 2;      % distance function derivative maximum degree of the approximation

zmax = 1.;

x_a   = [1,0];

% ========================================================================================
%% Preprocessing: nodes and mesh
fact= 1;
Nx  = 100;
Ny  = 100;
Lx  = fact*2;
Ly  = fact*2;

%x_n = [-1 0;1 0];       edges = [1 2];
%x_n = [-1 0;0 0;1 0];     edges = [1 2; 2 3];
% x_n = [-2  1;2  1;
%         2 -1;-2 -1];  edges = [2 1;3 2;4 3;1 4];
% x_n = fact*[-2  1;0  1; 2  1;2 0;
%             2 -1;0 -1;-2 -1;-0 0];  edges = [2 1;3 2;4 3;5 4;6 5;7 6;8 7;1 8];

delta=0.1;
x_n = fact*[-1  1;0  1; 1  1;1 0;
             1 -1;0 -1;-1 -1;-1 0];  edges = [2 1;3 2;4 3;5 4;6 5;7 6;8 7;1 8];
x_n = x_n + delta*[0 0;0  1;0 0; 1 0;
                   0 0;0 -1;0 0;-1 0];

% a0  = 0.0;
% x_n = fact*[1.0  0.0; 0.5  0.0; a0   a0; 
%             0.0  1.0;-1.0  1.0;-1.0  0.0; 
%            -1.0 -1.0;-0.5 -1.0; 0.0 -1.0; 
%             0.5 -1.0; 1.0 -1.0; 1.0 -0.5; 1.0 0.0];
% edges = [1 2;2 3; 3 4;4 5;5 6;6 7;7 8;8 9;9 10;10 11; 11 12;12 13];

Lx = max(Lx,Lx+2*delta);
Ly = max(Ly,Ly+2*delta);
x_s = UniformGrid2D(Nx,Ny,Lx,Ly,[0 0]);
% x_s = (-1.5:1.5:1.5);
% x_s = [x_s' zeros(length(x_s),1)];

figure(1)
plot(x_n(:,1),x_n(:,2),'ro','markersize',12,'linewidth',2)
hold on
plot(x_s(:,1),x_s(:,2),'bx')
% plot(x_a(1),x_a(2),'ks','markersize',15,'MarkerFaceColor','y','linewidth',2)
plot(x_s(2250,1),x_s(2250,2),'ks','markersize',15,'MarkerFaceColor','y','linewidth',2)
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'r*-','markersize',12,'linewidth',2)
end
hold off
axis equal

%% RFunction - ------------------------------------------------------
[ddist,gdist]= Rfunction_equiv_(edges,x_n,x_s,optCME.spow,optCME.mpow);
% ddist  = ddist/max(ddist);

in=inpolygon(x_s(:,1),x_s(:,2),x_n(:,1),x_n(:,2));
ddist(in==0)=0;
gdist(in==0,:)=0;
ddist_ = zeros(Ny,Nx);
gdist_ = zeros(Ny,Nx,2);
for i=1:Ny
  for j=1:Nx
    ddist_(i,j) = ddist((i-1)*Nx+j);
    gdist_(i,j,:) = gdist((i-1)*Nx+j,:);
  end
end

%% Plot basis functions at a given node

figure(2);clf
plot(x_n(:,1),x_n(:,2),'ro-','markersize',12,'linewidth',2)
hold on
%plot3(x_s(:,1),x_s(:,2),ddist,'bx')
surf(x_s(1:Nx,1),x_s(1:Nx:end,2),ddist_)
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'r*-','markersize',12,'linewidth',2)
end
hold off
axis normal
camlight
shading interp
set(gca,'FontName','times')
set(gca,'FontSize',20)
title('R-Function Equivalence')
xlabel('X')
ylabel('Y')
zlabel('Z')
zlim([0 zmax])
caxis([0.0 zmax])
colorbar('YLimMode','manual','Ylim',[0 zmax],'YTickMode','manual',...
  'YTick',[0.0,0.2,0.4,0.6,0.8,1.0],...
  'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'},...
  'FontSize',20,'FontName','times');


figure(3);clf
plot(x_n(:,1),x_n(:,2),'ro-','markersize',12,'linewidth',2)
hold on
%plot3(x_s(:,1),x_s(:,2),ddist,'bx')
surf(x_s(1:Nx,1),x_s(1:Nx:end,2),gdist_(:,:,1))
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'r*-','markersize',12,'linewidth',2)
end
hold off
axis normal
camlight
shading interp
set(gca,'FontName','times')
set(gca,'FontSize',20)
xlabel('X')
ylabel('Y')
zlabel('Z')
title('R-Function Equivalence - Grad X')

figure(4);clf
plot(x_n(:,1),x_n(:,2),'ro-','markersize',12,'linewidth',2)
hold on
%plot3(x_s(:,1),x_s(:,2),ddist,'bx')
surf(x_s(1:Nx,1),x_s(1:Nx:end,2),gdist_(:,:,2))
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'r*-','markersize',12,'linewidth',2)
end
hold off
axis normal
camlight
shading interp
set(gca,'FontName','times')
set(gca,'FontSize',20)
xlabel('X')
ylabel('Y')
zlabel('Z')
title('R-Function Equivalence - Grad Y')