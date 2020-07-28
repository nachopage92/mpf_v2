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

domain = 'circle24'; %'circle18'; %'circle24'

%  Input parameters: CME options
optCME.verb  = 0;       % 0:off    1:on
optCME.grad  = 1;       % Computation of the Gradient 0:OFF 1:ON
optCME.spow  = 2;       % w^spow, where w is the approximation to the distance function (R-function)
optCME.mpow  = 2;       % distance function derivative maximum degree of the approximation
optCME.nring = 2;       % number ring of neighbors to consider for each node: 1, 2, 3.
optCME.isConv= 1;       % flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)
optCME.polyCase='C';
optCME.pUnity= 1;       %partition of unity is enforced if is defined as 1
optCME.cleanColinear = 1; %if 1 then the colinear edges will be merged

% Set auxiliary variables for visualization
reflevel     = 6;         %level of refinement to plot basis functions and gradients
shrink_fct   = 1.e-3;     %boundary points are shrinked to the interior: shrink_bnd(x_n,ebnd,shrink_fct)
x_pfun       = [0 0]; %seed point for searching the nearest node whose shape function is computed

savePlot     = 0;         %save figures of shape functions and derivatives in png and fig format

fig_polyring = sprintf('%s_polyring_r%d',domain,optCME.nring);
fig_shpfun   = sprintf('%s_shpfunct_r%d_s%02d_m%02d_pu%d',domain,optCME.nring,optCME.spow,optCME.mpow,optCME.pUnity);
fig_DshpfunX = sprintf('%s_shpfunDX_r%d_s%02d_m%02d_pu%d',domain,optCME.nring,optCME.spow,optCME.mpow,optCME.pUnity);

limDP        = 1;        %maximum value of the gradient of p_nod in X and Y 
limP         = 0.5;      %maximum value of the basis functions of p_nod

% ========================================================================================
%% Preprocessing: nodes and mesh

% Node set, and some pre-processing
data = load(strcat('data/',domain,'_mesh.mat'));

x_n  = data.x_n;
tri  = data.tri;
ibnd = data.ibnd';

[tri,x_n] = repair_trimesh(tri,x_n);

nPts   = size(x_n,1);
nElem  = size(tri,1);
inside = ones(nElem,1);
C1     = ones(nPts, 1);

bPts   = length(ibnd);
x_bd   = x_n(ibnd,:);

figure(1);clf
trimesh(tri,x_n(:,1),x_n(:,2),C1,...
  'Facecolor','none','FaceLighting','none','EdgeColor','k','EdgeLighting','flat')
hold on
plot(x_n(:,1),x_n(:,2), ...
  'ko','MarkerFaceColor','b','Markersize',6 )
hold off
view([0,90])
axis equal


% ========================================================================================
%% Create the connectivity structures
t=cputime;
optCME.tri  = tri;
optCME.ibnd = ibnd;
connecStr = ConnectivityStructures_elem(x_n,optCME);

nodes_in_elem = connecStr.nodes_in_elem;
bd_segm       = connecStr.bd_segm;
bd_segm_out   = connecStr.bd_segm_out;
bd_nodes      = connecStr.bd_nodes;
first_ring    = connecStr.first_ring;
last_ring     = connecStr.last_ring;

optCME.connec_struct = connecStr;
optCME.bd_segm     = bd_segm;
optCME.bd_segm_out = bd_segm_out;
optCME.bd_nodes    = bd_nodes;

fprintf(1,'cputime CME Connectivity Structures: %5.3f seconds\n', cputime-t);

% ========================================================================================
%% Display graphically the neighborhoods

%We define a non-convex domain
figure(2);clf
hold on
plot(x_n(:,1),x_n(:,2), ...
  'ko','MarkerFaceColor','b','Markersize',10 )
hold on

%We define a non-convex domain
ebnd   = freeBoundary(triangulation(tri,x_n));
for k = 1:size(ebnd,1)
  plot(x_n(ebnd(k,:),1),x_n(ebnd(k,:),2),'r-','LineWidth',4)
end
plot(x_bd(:,1),x_bd(:,2), ...
  'ko','MarkerFaceColor','c','Markersize',8 )

%Triangulation constrained by the boundary edges
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'m-','LineWidth',1);
axis equal

if 1 == 0
  % Plot the elements influenced by a node and the nodes participating in prior
  for i=1:length(x_n)
    h11=triplot(tri(last_ring{i}, :),x_n(:,1),x_n(:,2),'k-','LineWidth',2);
    neighbornodes = unique(tri(last_ring{i},:));
    h2=plot(x_n(neighbornodes,1),x_n(neighbornodes,2),...
      'go','Markersize',10 ,'LineWidth',2);
    h3= [];
    for ii = 1:size(bd_segm{i},1)
      h3(ii)=plot(x_n(bd_segm{i}(ii,:),1),x_n(bd_segm{i}(ii,:),2),...
        'rx-','Markersize',10 ,'LineWidth',4);
    end
    %h33= [];
    %for ii = 1:size(bd_segm_out{i},1)
    %    h33(ii)=plot(x_n(bd_segm_out{i}(ii,:),1),x_n(bd_segm_out{i}(ii,:),2),...
    %        'mo-','Markersize',10 ,'LineWidth',4);
    %end
    
    h1=plot(x_n(i,1),x_n(i,2), 'k*','Markersize',20 ,'LineWidth',3);
    axis equal
    pause
    delete(h1)
    delete(h11)
    delete(h2)
    delete(h3)
    %delete(h33)
  end
  
  % Plot the nodes affecting an element
  for j = 1:numel_
    if inside(j)==1
      h1=triplot(tri(j, :),x_n(:,1),x_n(:,2),'c-','LineWidth',4);
      neighbornodes = nodes_in_elem{j};
      h2=plot(x_n(neighbornodes,1),x_n(neighbornodes,2),...
        'ko','Markersize',10 ,'LineWidth',4);
      pause
      delete(h1)
      delete(h2)
    end
  end
  return
end

% ========================================================================================
%% Compute basis functions in a given element
[~,inod]=min(dist(x_pfun,x_n'));

ie_near = zeros(length(first_ring{inod}),1);
for i=1:length(ibnd)
  tri_n   = tri(first_ring{inod},:);
  ie_near = ie_near + any(tri_n(:,:)==ibnd(i),2);
end
[nn_max,iel]=max(ie_near);
iel     = first_ring{inod}(iel);

enear= nodes_in_elem{iel}; % Nodes affecting the element
optCME.enear = enear;

if inside(iel) == 1
  % A set of evaluation points within this element
  x_s = [];
  x_s(1,:) = x_n(tri(iel,1),:)/6 + x_n(tri(iel,2),:)/6 + x_n(tri(iel,3),:)*2/3 ;
  x_s(2,:) = x_n(tri(iel,1),:)/6 + x_n(tri(iel,3),:)/6 + x_n(tri(iel,2),:)*2/3 ;
  x_s(3,:) = x_n(tri(iel,3),:)/6 + x_n(tri(iel,2),:)/6 + x_n(tri(iel,1),:)*2/3 ;
  tic
  
  %------------------------------------------------------------------
  % we compute the prior
  plot(x_n(tri(iel,:),1),x_n(tri(iel,:),2),'ks',...
    x_n(enear,1),x_n(enear,2),'go',...
    x_s(:,1),x_s(:,2),'bx','Markersize',10 ,'LineWidth',2);
  
  [q,dq]=priorf_elem(x_n,x_s,optCME);
  toc
end


% ========================================================================================
%% Check derivatives
if 1 == 1
  h = 1.e-7;
  % A set of evaluation points within this element
  x_s = [];
  x_s(1,:) = x_n(tri(iel,1),:)/6 + x_n(tri(iel,2),:)/6 + x_n(tri(iel,3),:)*2/3 ;
  [q,dq]=priorf_elem(x_n,x_s,optCME);
  
  % X-axis
  x_s_p = x_s + [h 0];
  q_p   = priorf_elem(x_n,x_s_p,optCME);
  
  x_s_m = x_s - [h 0];
  q_m   = priorf_elem(x_n,x_s_m,optCME);
  dq_num = (q_p-q_m)/2/h;
  disp('Error in the derivatives: X-axis')
  disp(norm(dq_num-dq(:,:,1)))
  
  % Y-axis
  x_s_p = x_s + [0 h];
  q_p   = priorf_elem(x_n,x_s_p,optCME);
  
  x_s_m = x_s - [0 h];
  q_m   = priorf_elem(x_n,x_s_m,optCME);
  dq_num = (q_p-q_m)/2/h;
  disp('Error in the derivatives: Y-axis')
  disp(norm(dq_num-dq(:,:,2)))
end


% ========================================================================================
%% Plot a given shape function
[~,inod]=min(dist(x_pfun,x_n'));

segmIN  = bd_segm{inod};
segmOUT = bd_segm_out{inod};
segm1   = segmOUT(:,1);
segm2   = segmOUT(:,2);
ijoin   = intersect(segm1,segm2);


% -------------------------------------------------------------------
% Here is generated a shrinked node set such that only the boundary nodes are
% moved to the interior in the normal direction to the boundary edge
x_b = shrink_bnd(x_n,ebnd,shrink_fct);

figure(10);clf
plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','b','Markersize',6 )
hold on
plot(x_b(:,1),x_b(:,2),'r.',x_b(ibnd,1),x_b(ibnd,2),'m*')
plot(x_n(ijoin,1),x_n(ijoin,2),'bo','MarkerFaceColor','c','Markersize',18,'LineWidth',1)
plot(x_n(inod,1),x_n(inod,2),'ks','MarkerFaceColor','y','Markersize',15,'LineWidth',2)
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'m-','LineWidth',1);
for ek=1:size(segmOUT,1)
  plot(x_n(segmOUT(ek,:),1),x_n(segmOUT(ek,:),2),'r-','LineWidth',4)
end
for ek=1:size(segmIN,1)
  plot(x_n(segmIN(ek,:),1),x_n(segmIN(ek,:),2),'kd-','MarkerFaceColor','g','Markersize',12,'LineWidth',2)
end
hold off
axis equal


%% -------------------------------------------------------------------
clear  q_nod
clear dq_nod
clear tri_r_
clear x_r_
for i = 1: length(last_ring{inod})
  iel  = last_ring{inod}(i);
  enear= nodes_in_elem{iel};
  
  x_t = x_b(tri(iel,:),:);
  
  x_t = [x_t, [0 0 0]'];
  tri_el = [1 2 3];
  [x_r,tri_r] = subdivide_tri( x_t, tri_el);
  for ir=1:reflevel
    [x_r,tri_r] = subdivide_tri( x_r, tri_r);
  end
  x_r(:,3) =[];
  tri_r_(i)={tri_r};
  x_r_(i)  ={x_r};
  tic
  optCME.enear = enear;
  [q,dq]=priorf_elem(x_n,x_r,optCME);
  toc
  q_nod(i)={q(:,enear==inod)};
  dq_nod_(:,1)=dq(:,enear==inod,1);
  dq_nod_(:,2)=dq(:,enear==inod,2);
  dq_nod(i)={dq_nod_};
end


% ========================================================================================
%% Plot a given shape function
% -------------------------------------------------------------------
% Here is generated a shrinked node set such that only the boundary nodes are
% moved to the interior in the normal direction to the boundary edge
segmIN  = bd_segm{inod};
segmOUT = bd_segm_out{inod};
segm1   = segmOUT(:,1);
segm2   = segmOUT(:,2);
ijoin   = intersect(segm1,segm2);

figure(10);clf
plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','k','Markersize',6 )
hold on
%plot(x_b(:,1),x_b(:,2),'r.',x_b(ibnd,1),x_b(ibnd,2),'m*')
% plot(x_n(ijoin,1),x_n(ijoin,2),'bo','MarkerFaceColor','c','Markersize',18,'LineWidth',2)
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'k-','LineWidth',1);
if optCME.isConv == 1
  for ek=1:size(segmIN,1)
    plot(x_n(segmIN(ek,:),1),x_n(segmIN(ek,:),2),'k-','LineWidth',6)
    plot(x_n(segmIN(ek,:),1),x_n(segmIN(ek,:),2),'ko','MarkerFaceColor','y','Markersize',12,'LineWidth',2)
  end
  title(strcat('Domain is considered as Convex'),'fontsize',20)
else
  for ek=1:size(segmOUT,1)
    plot(x_n(segmOUT(ek,:),1),x_n(segmOUT(ek,:),2),'r-','LineWidth',6)
  end
  for ek=1:size(segmIN,1)
    plot(x_n(segmIN(ek,:),1),x_n(segmIN(ek,:),2),'k--','LineWidth',3)
    plot(x_n(segmIN(ek,:),1),x_n(segmIN(ek,:),2),'kd','MarkerFaceColor','g','Markersize',12,'LineWidth',2)
  end
  title(strcat('Domain is considered as NON Convex'),'fontsize',20)
end
plot(x_n(inod,1),x_n(inod,2),'ks','MarkerFaceColor','y','Markersize',16,'LineWidth',2)

%plot(x_r(:,1),x_r(:,2),'bx',x_n(enear(sum(q)==0),1),x_n(enear(sum(q)==0),2),'rv','MarkerFaceColor','g','Markersize',15,'LineWidth',2)
hold off
axis equal
set(gca,'XTickLabel',{},'YTickLabel',{},'FontSize',20,'FontName','times');
xlim([min(x_n(:,1)) max(x_n(:,1))])
ylim([min(x_n(:,2)) max(x_n(:,2))])

if savePlot==1
  print('-dpng','-r300',strcat('figures_png/',fig_polyring,'.png'))
end

%% -------------------------------------------------------------------
%Plot basis functions at a given node
plot_domain(x_bd,x_n,ebnd);
hold on
for i = 1: length(last_ring{inod})
  C = fix(q_nod{i}*70)+10;
  trisurf(tri_r_{i},x_r_{i}(:,1),x_r_{i}(:,2),q_nod{i},C,'CDataMapping','direct')
end
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'k-','LineWidth',1);
view([30 25])
axis normal
camlight
shading interp
lighting phong
%material metal
material shiny
%colormap([0.7 0.7 0.7])
set(gca,'ZTickLabel',{'0','0.25','0.50','0.75'},'ZTick',[0 0.25 0.5 0.75],...
  'XTickLabel',{'-1','0','1'},'XTick',[-1,0,1],...
  'YTickLabel',{'-1','0','1'},'YTick',[-1,0,1],...
  'FontSize',20,...
  'FontName','times',...
  'LineWidth',1,'box','off');
% xlabel('X')
% ylabel('Y')
xlim([min(x_n(:,1)) max(x_n(:,1))])
ylim([min(x_n(:,2)) max(x_n(:,2))])
zlim([0 limP])
if savePlot==1
  print('-dpng','-r300',strcat('figures_png/',fig_shpfun,'.png'))
end


%% -------------------------------------------------------------------
%Plot gradient in X of the basis functions at a given node
plot_domain(x_bd,x_n,ebnd);
hold on
for i = 1: length(last_ring{inod})
  C1 = round((dq_nod{i}(:,1))*32)+32;
  trisurf(tri_r_{i},x_r_{i}(:,1),x_r_{i}(:,2),dq_nod{i}(:,1),C1,'CDataMapping','direct')
end
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'k-','LineWidth',1);
view([30 25])
axis normal
camlight
shading interp
set(gca,'ZTickLabel',{sprintf('-%d',limDP),'0',sprintf('-%d',limDP)},'ZTick',[-limDP 0 limDP],...
  'XTickLabel',{'-1','0','1'},'XTick',[-1,0,1],...
  'YTickLabel',{'-1','0','1'},'YTick',[-1,0,1],...
  'FontSize',20,...
  'FontName','times',...
  'LineWidth',1,'box','off');
xlim([min(x_n(:,1)) max(x_n(:,1))])
ylim([min(x_n(:,2)) max(x_n(:,2))])
zlim(1.1*[-limDP limDP])
%xlabel('X')
%ylabel('Y')
%title('Gradient X')
if savePlot==1
  print('-dpng','-r300',strcat('figures_png/',fig_DshpfunX,'.png'))
end
