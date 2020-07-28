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

% Preprocessing: nodes and mesh

% %Example 1 -----------------------------------------
Nx  = 300;
Ny  = 150;
Lx  = 2;
Ly  = 1;
x_n = [-0.5 0;0.5 0];       edges = [1 2];

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

%% TRIMMING Functions - ------------------------------------------------------
[f,df,t,dt] = trimming(x_n(1,:),x_n(2,:),x_s);

fdist_ = zeros(Ny,Nx);
tdist_ = zeros(Ny,Nx);
rdist_ = zeros(Ny,Nx);
for i=1:Ny
  for j=1:Nx
    ff = f((i-1)*Nx+j);
    tt = t((i-1)*Nx+j);
    fdist_(i,j) = ff;
    tdist_(i,j) = tt;
    ff2 = ff*ff;
    val = sqrt(tt*tt+ff2*ff2)-tt;
    rdist_(i,j) = sqrt(ff2 + 0.25*val*val);
  end
end


fmin = -0.5;  fmax = 0.50;
tmin = -1.0;  tmax = 0.25;
rmin =  0.0;  rmax = 1.00;

fstep = (fmax-fmin)/20; fisolins = (fmin:fstep:fmax);
tstep = (tmax-tmin)/20; tisolins = (tmin:tstep:tmax);
rstep = (rmax-rmin)/20; risolins = (rmin:rstep:rmax);

%% Plot f-distance functions
figure1 = figure('Renderer','OpenGL','InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[4.5 4.5 1],'Layer','top','FontSize',24,...
  'FontName','times',...
    'XTick',[-1.5,-1,-0.5,0,0.5,1,1.5],...
  'XTickLabel',{'-1.5','-1','-0.5','0','0.5','1','1.5'},...
  'YTick',[-1,-0.5,0,0.5,1],...
  'YTickLabel',{'-1','-0.5','0','0.5','1'},...
  'DataAspectRatio',[1 1 1]);

box(axes1,'on');
hold(axes1,'all');

xlim([-Lx/2 Lx/2])
ylim([-Ly/2 Ly/2])

% Create contour
contour(x_s(1:Nx,1),x_s(1:Nx:end,2),fdist_,'LineWidth',1,'LevelStep',0.001,...
  'Parent',axes1,'Fill','on','LevelList',fisolins);
hold on
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'k-','markersize',12,'linewidth',8)
end
hold off

% Create light
light('Parent',axes1,'Style','local',...
  'Position',[37.2617031399255 42.1752557155139 55.0128835732709]);

set(gca,'FontName','times')
set(gca,'FontSize',24)
%title('R-Function Equivalence')
%xlabel('X','FontSize',24);
%ylabel('Y','FontSize',24);


% Create colorbar
%colorbar('peer',axes1,'FontSize',20,'FontName','times');
caxis([fmin fmax])
colorbar('YLimMode','manual','Ylim',[fmin fmax],...
  'FontSize',24,'FontName','times');
%colorbar('YLimMode','manual','Ylim',[tmin tmax],'YTickMode','manual',...
%  'YTick',[0.0,0.2,0.4,0.6,0.8,1.0],...
%  'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'},...
%  'FontSize',24,'FontName','times');


%% Plot t-distance functions
figure2 = figure('Renderer','OpenGL','InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes2 = axes('Parent',figure2,'PlotBoxAspectRatio',[4.5 4.5 1],'Layer','top','FontSize',24,...
  'FontName','times',...
    'XTick',[-1.5,-1,-0.5,0,0.5,1,1.5],...
  'XTickLabel',{'-1.5','-1','-0.5','0','0.5','1','1.5'},...
  'YTick',[-1,-0.5,0,0.5,1],...
  'YTickLabel',{'-1','-0.5','0','0.5','1'},...
  'DataAspectRatio',[1 1 1]);

box(axes2,'on');
hold(axes2,'all');

xlim([-Lx/2 Lx/2])
ylim([-Ly/2 Ly/2])

% Create contour
contour(x_s(1:Nx,1),x_s(1:Nx:end,2),tdist_,'LineWidth',1,'LevelStep',0.001,...
  'Parent',axes2,'Fill','on','LevelList',tisolins);
hold on
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'k-','markersize',12,'linewidth',8)
end
hold off

% Create light
light('Parent',axes2,'Style','local',...
  'Position',[37.2617031399255 42.1752557155139 55.0128835732709]);

set(gca,'FontName','times')
set(gca,'FontSize',24)
%title('R-Function Equivalence')
%xlabel('X','FontSize',24);
%ylabel('Y','FontSize',24);


% Create colorbar
%colorbar('peer',axes1,'FontSize',20,'FontName','times');
caxis([tmin tmax])
colorbar('YLimMode','manual','Ylim',[tmin tmax],...
  'FontSize',24,'FontName','times');
% colorbar('YLimMode','manual','Ylim',[tmin tmax],'YTickMode','manual',...
%   'YTick',[0.0,0.2,0.4,0.6,0.8,1.0],...
%   'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'},...
%   'FontSize',24,'FontName','times');


%% Plot rho-distance functions
figure3 = figure('Renderer','OpenGL','InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes3 = axes('Parent',figure3,'PlotBoxAspectRatio',[4.5 4.5 1],'Layer','top','FontSize',24,...
  'FontName','times',...
    'XTick',[-1.5,-1,-0.5,0,0.5,1,1.5],...
  'XTickLabel',{'-1.5','-1','-0.5','0','0.5','1','1.5'},...
  'YTick',[-1,-0.5,0,0.5,1],...
  'YTickLabel',{'-1','-0.5','0','0.5','1'},...
  'DataAspectRatio',[1 1 1]);

box(axes3,'on');
hold(axes3,'all');

xlim([-Lx/2 Lx/2])
ylim([-Ly/2 Ly/2])

% Create contour
contour(x_s(1:Nx,1),x_s(1:Nx:end,2),rdist_,'LineWidth',1,'LevelStep',0.001,...
  'Parent',axes3,'Fill','on','LevelList',risolins);
hold on
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'k-','markersize',12,'linewidth',6)
end
hold off

% Create light
light('Parent',axes3,'Style','local',...
  'Position',[37.2617031399255 42.1752557155139 55.0128835732709]);

set(gca,'FontName','times')
set(gca,'FontSize',24)
%title('R-Function Equivalence')
%xlabel('X','FontSize',24);
%ylabel('Y','FontSize',24);


% Create colorbar
%colorbar('peer',axes1,'FontSize',20,'FontName','times');
caxis([rmin rmax])
colorbar('YLimMode','manual','Ylim',[rmin rmax],...
  'FontSize',24,'FontName','times');
% colorbar('YLimMode','manual','Ylim',[rmin rmax],'YTickMode','manual',...
%   'YTick',[0.0,0.2,0.4,0.6,0.8,1.0],...
%   'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'},...
%   'FontSize',24,'FontName','times');
