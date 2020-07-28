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
spow  = 1;      % w^spow, where w is the approximation to the distance function (R-function)
mpow  = 6;       % distance function derivative maximum degree of the approximation

zmax = 2;

x_a   = [1,0];

% ========================================================================================
%% Preprocessing: nodes and mesh

% Example 2 -----------------------------------------
Nx  = 300;
Ny  = 200;
Lx  = 3;
Ly  = 2;
x_n = [-1 -0.5;0 0.5;1 -0.5];     edges = [1 2; 2 3];


x_s = UniformGrid2D(Nx,Ny,Lx,Ly,[0 0]);
% x_s = (-1.5:1.5:1.5);
% x_s = [x_s' zeros(length(x_s),1)];

figure(1)
plot(x_n(:,1),x_n(:,2),'ro','markersize',12,'linewidth',2)
hold on
plot(x_s(:,1),x_s(:,2),'bx')
% plot(x_a(1),x_a(2),'ks','markersize',15,'MarkerFaceColor','y','linewidth',2)
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'r*-','markersize',12,'linewidth',4)
end
hold off
axis equal

%% RFunction - ------------------------------------------------------
ddist = Rfunction_equiv_(edges,x_n,x_s,spow,mpow);

ddist_ = zeros(Ny,Nx);
for i=1:Ny
  for j=1:Nx
    ddist_(i,j) = ddist((i-1)*Nx+j);
  end
end

%% Plot basis functions at a given node

figure1 = figure('Renderer','OpenGL','InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[4.5 4.5 1],'Layer','top','FontSize',30,...
  'FontName','times',...
  'XTick',[-1.5,-1,-0.5,0,0.5,1,1.5],...
  'XTickLabel',{'-1.5','-1','-0.5','0','0.5','1','1.5'},...
  'YTick',[-1,-0.5,0,0.5,1],...
  'YTickLabel',{'-1','-0.5','0','0.5','1'},...
  'DataAspectRatio',[1 1 1]);

box(axes1,'on');
hold(axes1,'all');

isolins = (0:0.1:zmax);

xlim([-Lx/2 Lx/2])
ylim([-Ly/2 Ly/2])

% Create contour
contour(x_s(1:Nx,1),x_s(1:Nx:end,2),ddist_,'LineWidth',1,'LevelStep',0.001,...
  'LevelList',isolins,...
  'Fill','on',...
  'Parent',axes1);
hold on
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'k-',...
      'markersize',12,'linewidth',6,'MarkerFaceColor','k')
end
hold off

% Create light
light('Parent',axes1,'Style','local',...
  'Position',[37.2617031399255 42.1752557155139 55.0128835732709]);

set(gca,'FontName','times')
set(gca,'FontSize',30)
%title('R-Function Equivalence')
%xlabel('X','FontSize',30);
%ylabel('Y','FontSize',30);


% Create colorbar
%colorbar('peer',axes1,'FontSize',20,'FontName','times');
caxis([0.0 zmax])
colorbar('YLimMode','manual','Ylim',[0 zmax],'YTickMode','manual',...
  'YTick',[0.0,0.5,1.0,1.5,2.0],...
  'YTickLabel',{'0.0','0.5','1.0','1.5','2.0'},...
  'FontSize',30,'FontName','times');
