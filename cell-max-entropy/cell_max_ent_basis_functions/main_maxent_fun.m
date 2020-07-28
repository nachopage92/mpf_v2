%%
% This function plots the element-based maximum entropy (CME) approximants
% for a given node set, as well as the gradients and the approximation
% for discontinuous nodal data
%
% Reference:
% [1] Marino Arroyo and Michael Ortiz, Local maximum-entropy approximation
%     schemes: a seamless bridge between finite elements and meshfree methods,
%     International Journal for Numerical Methods in Engineering, 65:2167-2202 (2006).


clear all
%close all

format long

% ========================================================================================
%% Input parameters

domain = 'circle';%'grid09';%'qhole'; %'hole'; %'slice060';'quoin090';%'ele'; %'grid' 'circle'

%  Input parameters: CME options
optCME.verb  = 0;       % 0:off    1:on
optCME.grad  = 1;       % Computation of the Gradient 0:OFF 1:ON
optCME.TolNR = 1.e-12;  % Newton-Raphson tolerance fo the nonlinear max-ent problem

optCME.nring = 2; %R    % number ring of neighbors to consider for each node: 1, 2, 3.
optCME.spow  = 2; %S    % w^spow, where w is the approximation to the distance function (R-function)
optCME.mpow  = 2; %M    % distance function derivative maximum degree of the approximation
optCME.isConv= 1;       % flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)
optCME.polyCase='C';

% Set auxiliary variables for visualization ----------------------------------------------

%seed point for searching the nearest node whose shape function is computed
x_pfun       = [0.0 0.0];

reflevel     = 3;        %level of refinement to plot basis functions and gradients
shrink_fct   = 1e-3;     %boundary points are shrinked to the interior: shrink_bnd(x_n,ebnd,shrink_fct)

savePlot     = 0;        %save figures of shape functions and derivatives in png and fig format
saveVTK      = 0;
limDP        = 2;        %maximum value of the gradient of p_nod in X and Y 
limP         = 0.5;      %maximum value of the basis functions of p_nod

fig_polyring = sprintf('%s_polyring_r%d_s%02d_m%02d_ptA_%s%d',domain,optCME.nring,optCME.spow,optCME.mpow,optCME.isConv,optCME.polyCase);
fig_shpfun   = sprintf('%s_shpfunct_r%d_s%02d_m%02d_ptA_%s%d',domain,optCME.nring,optCME.spow,optCME.mpow,optCME.isConv,optCME.polyCase);
fig_DshpfunX = sprintf('%s_shpfunDX_r%d_s%02d_m%02d_ptA_%s%d',domain,optCME.nring,optCME.spow,optCME.mpow,optCME.isConv,optCME.polyCase);
fig_DshpfunY = sprintf('%s_shpfunDY_r%d_s%02d_m%02d_ptA_%s%d',domain,optCME.nring,optCME.spow,optCME.mpow,optCME.isConv,optCME.polyCase);

% ========================================================================================
%% Preprocessing: nodes and mesh

% Node set, and some pre-processing
data=load(strcat('data/',domain,'_mesh.mat'));

theta= 0;
Rot  = [cos(theta) sin(theta);
        -sin(theta) cos(theta)];
x_n  = data.p*Rot;
%x_n  = data.x_n*Rot;
x_pfun=x_pfun*Rot;
tri  = data.tri;
[tri,x_n] = repair_trimesh(tri,x_n);

nPts = size(x_n,1);
C1   = ones(nPts, 1);

ids  = 1:nPts;

% -------------------------------------------------------------------
% Boundary nodes identifiers
ebnd   = freeBoundary(triangulation(tri,x_n));
ibnd   = unique(ebnd);
bPts   = length(ibnd);
x_bd   = x_n(ibnd,:);

nElem   = size(tri,1);
badelem = zeros(nElem,1);


%checking wrong elements with zero area and deleting them
for e=1:nElem
  if TriangleArea(x_n(tri(e,:),:)) < 1.0e-06
    badelem(e) = 1;
  end
end
fprintf(1,'Bad elements with area < 1e-06:  %d\n',sum(badelem));
tri(badelem==1,:)  = [];
nElem  = size(tri,1);
inside = ones(nElem,1);

figure(1);clf
trimesh(tri,x_n(:,1),x_n(:,2),C1,...
  'Facecolor','none','FaceLighting','none','EdgeColor','k','EdgeLighting','flat')
hold on
plot(x_n(:,1),x_n(:,2), ...
  'ko','MarkerFaceColor','b','Markersize',6 )
hold off
view([0,90])
axis equal


%We define a non-convex domain with alpha shapes
figure(2);clf
hold on
plot(x_n(:,1),x_n(:,2), ...
  'ko','MarkerFaceColor','b','Markersize',8 )
hold on

%We define a non-convex domain
%ebnd = freeBoundary(triangulation(tri,x_n));
for k = 1:size(ebnd,1)
  plot(x_n(ebnd(k,:),1),x_n(ebnd(k,:),2),'r-','LineWidth',3)
end
plot(x_bd(:,1),x_bd(:,2), ...
  'ko','MarkerFaceColor','c','Markersize',8 )

%Triangulation constrained by the boundary edges
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'m-','LineWidth',1);
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

TolLag = optCME.TolNR;
spow   = optCME.spow;
mpow   = optCME.mpow;

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

optCME.enear = nodes_in_elem{iel};
if inside(iel) == 1
  % A set of evaluation points within this element
  x_s = [];
  x_s(1,:) = x_n(tri(iel,1),:)/6 + x_n(tri(iel,2),:)/6 + x_n(tri(iel,3),:)*2/3 ;
  x_s(2,:) = x_n(tri(iel,1),:)/6 + x_n(tri(iel,3),:)/6 + x_n(tri(iel,2),:)*2/3 ;
  x_s(3,:) = x_n(tri(iel,3),:)/6 + x_n(tri(iel,2),:)/6 + x_n(tri(iel,1),:)*2/3 ;
  tic
  
  %------------------------------------------------------------------
  % we compute the prior
  near= nodes_in_elem{iel}; % Nearest Nodes affecting the element
  id_el = [tri(iel,:) tri(iel,1)];

  plot(x_n(near,1),x_n(near,2),'ks','MarkerFaceColor','g','Markersize',10 ,'LineWidth',2);
  plot(x_n(id_el,1),x_n(id_el,2),'r*-',...
    x_s(:,1),x_s(:,2),'bx','Markersize',10 ,'LineWidth',3);
  
  % prior Rfunctions ------------------------------------------------------
  [q,dq]=priorf_elem(x_n,x_s,optCME);
  optCME.q  = q;
  optCME.dq = dq;
  
  % cell-based max-ent functions ------------------------------------------
  [p,dp]=shapef_elem(x_n,x_s,optCME);
  toc
end


% ========================================================================================
%% Check derivatives
if 1 == 1
  h = 1.e-7;
  TolCHECK=1.e-12;
  % A set of evaluation points within this element
  x_s = [];
  x_s(1,:) = x_n(tri(iel,1),:)/6 + x_n(tri(iel,2),:)/6 + x_n(tri(iel,3),:)*2/3 ;
  x_s(2,:) = x_n(tri(iel,1),:)/6 + x_n(tri(iel,3),:)/6 + x_n(tri(iel,2),:)*2/3 ;
  x_s(3,:) = x_n(tri(iel,3),:)/6 + x_n(tri(iel,2),:)/6 + x_n(tri(iel,1),:)*2/3 ;
  
  %   [q,dq]=priorf_elem(x_n,x_s,optCME);
  %   optCME.q  = q;
  %   optCME.dq = dq;
  %   [p,dp]=shapef_elem(x_n,x_s,optCME);
  [p,dp,q,dq]=shapef_elem_(x_n,x_s,optCME);
  
  % X-axis
  x_s_p = x_s + [1;1;1]*[h 0];
  %   [q_p,dq_p]=priorf_elem(x_n,x_s_p,optCME);
  %   optCME.q  = q_p;
  %   optCME.dq = dq_p;
  %   [p_p,dp_p]=shapef_elem(x_n,x_s_p,optCME);
  [p_p,dp_p,q_p,dq_p]=shapef_elem_(x_n,x_s_p,optCME);
  
  x_s_m = x_s - [1;1;1]*[h 0];
  %   [q_m,dq_m]=priorf_elem(x_n,x_s_m,optCME);
  %   optCME.q  = q_m;
  %   optCME.dq = dq_m;
  %   [p_m,dp_m]=shapef_elem(x_n,x_s_m,optCME);
  [p_m,dp_m,q_m,dq_m]=shapef_elem_(x_n,x_s_m,optCME);
  dq_num = (q_p-q_m)/2/h;
  dp_num = (p_p-p_m)/2/h;
  disp('Error in the derivatives: X-axis')
  disp([norm(dq_num-dq(:,:,1)) norm(dp_num-dp(:,:,1))])
  
  % Y-axis
  x_s_p = x_s + [1;1;1]*[0 h];
  %   [q_p,dq_p]=priorf_elem(x_n,x_s_p,optCME);
  %   optCME.q  = q_p;
  %   optCME.dq = dq_p;
  %   [p_p,dp_p]=shapef_elem(x_n,x_s_p,optCME);
  [p_p,dp_p,q_p,dq_p]=shapef_elem_(x_n,x_s_p,optCME);  
  
  x_s_m = x_s - [1;1;1]*[0 h];
  %   [q_m,dq_m]=priorf_elem(x_n,x_s_m,optCME);
  %   optCME.q  = q_m;
  %   optCME.dq = dq_m;
  %   [p_m,dp_m]=shapef_elem(x_n,x_s_m,optCME);
  [p_m,dp_m,q_m,dq_m]=shapef_elem_(x_n,x_s_m,optCME);
  dq_num = (q_p-q_m)/2/h;
  dp_num = (p_p-p_m)/2/h;
  disp('Error in the derivatives: Y-axis')
  disp([norm(dq_num-dq(:,:,2)) norm(dp_num-dp(:,:,2))])
end


% ========================================================================================
%% Plot a given shape function
% -------------------------------------------------------------------
% Here is generated a shrinked node set such that only the boundary nodes are
% moved to the interior in the normal direction to the boundary edge
if optCME.nring>1
  x_b = shrink_bnd(x_n,ebnd,shrink_fct);
else
  x_b = x_n;
end

segmIN  = bd_segm{inod};
segmOUT = bd_segm_out{inod};
segm1   = segmOUT(:,1);
segm2   = segmOUT(:,2);
ijoin   = intersect(segm1,segm2);

figure(10);clf
plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','b','Markersize',6 )
hold on
%plot(x_b(:,1),x_b(:,2),'r.',x_b(ibnd,1),x_b(ibnd,2),'m*')
plot(x_n(ijoin,1),x_n(ijoin,2),'bo','MarkerFaceColor','c','Markersize',18,'LineWidth',2)
plot(x_n(inod,1),x_n(inod,2),'ks','MarkerFaceColor','y','Markersize',15,'LineWidth',2)
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'k-','LineWidth',1);
if optCME.isConv == 1
  for ek=1:size(segmIN,1)
    plot(x_n(segmIN(ek,:),1),x_n(segmIN(ek,:),2),'m--','LineWidth',4)
    plot(x_n(segmIN(ek,:),1),x_n(segmIN(ek,:),2),'md','MarkerFaceColor','g','Markersize',12,'LineWidth',2)
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
hold off
axis equal
set(gca,'XTickLabel',{},'YTickLabel',{},'FontSize',20,'FontName','times');
xlim([min(x_n(:,1)) max(x_n(:,1))])
ylim([min(x_n(:,2)) max(x_n(:,2))])

if savePlot==1
  print('-dpng','-r300',strcat('figures_png_poly/',fig_polyring,'.png'))
end

%% -------------------------------------------------------------------
tic
clear  p_nod
clear dp_nod
clear tri_r_
clear x_r_

%length of the domain bounding box
dx    = abs(max(x_n(:,1))-min(x_n(:,1)));
dy    = abs(max(x_n(:,2))-min(x_n(:,2)));
bbox  = sqrt(dx*dx + dy*dy);

nfail=0;
xs_f = [];
for i = 1: length(last_ring{inod})
  iel = last_ring{inod}(i);
  x_t = x_b(tri(iel,:),:);
  
  if optCME.nring ==1
    for j = 1:3
      x_c=sum(x_n(tri(iel,:),:),1)/3;
      x_t(j,:)=x_t(j,:) + (x_c-x_t(j,:))*bbox/1000;
    end
  end
  
  x_t = [x_t, [0 0 0]'];
  tri_el = [1 2 3];
  [x_r,tri_r] = subdivide_tri( x_t, tri_el);
  for ir=1:reflevel
    [x_r,tri_r] = subdivide_tri( x_r, tri_r);
  end
  x_r(:,3) =[];
  tri_r_(i)={tri_r};
  x_r_(i)  ={x_r};
  
  optCME.enear = nodes_in_elem{iel};
  
  % prior Rfunctions ------------------------------------------------------
%   [q,dq]=priorf_elem(x_n,x_r,optCME);
%   optCME.q  = q;
%   optCME.dq = dq;
  
  % cell-based max-ent functions ------------------------------------------
  [p,dp,q,dq,ifail]=shapef_elem_(x_n,x_r,optCME);
  
  if mod(i,round(length(last_ring{inod})/10))==0 || i==length(last_ring{inod})
    fprintf(1,'%5.1f%% -- Elapsed time is %6.2f seconds\n',i/length(last_ring{inod})*100, toc);
  end
  
  p_nod(i)={p(:,nodes_in_elem{iel}==inod)};
  dp_nod_(:,1)=dp(:,nodes_in_elem{iel}==inod,1);
  dp_nod_(:,2)=dp(:,nodes_in_elem{iel}==inod,2);
  dp_nod(i)={dp_nod_};
  
  %Check 0th and 1st order consistency
  err0 = max(abs(sum(p,2)-1));
  kk   = p*x_n(nodes_in_elem{iel},:);
  err1 = max(max(abs(kk-x_r)));
  err2 = max(abs(sum(dp(:,:,1),2))) + max(abs(sum(dp(:,:,2),2)));
  if err0 > TolLag
    disp('Error in 0th order consistency')
    disp(err0)
  end
  if err1 > TolLag*1e2
    disp('Error in 1st order consistency')
    disp(err1)
  end
  if err2 > TolLag*1e4
    disp('Error in sum grad p = 0')
    disp(err2)
  end
  nfail = nfail + sum(ifail);
  xs_f  = [xs_f;x_r(ifail==1,:)];
end
fprintf(1,'The number of sample points where the N-R fails is =%d\n',nfail);

if nfail>0
  figure(10)
  hold on
  plot(xs_f(:,1),xs_f(:,2),'rd','MarkerFaceColor','g','Markersize',12,'LineWidth',2)
  hold off
end

%% -------------------------------------------------------------------
%Plot basis functions at a given node
p_max   = 0;
idnears = tri(last_ring{inod},:);
idnears = unique(idnears(:));

plot_domain(x_bd,x_n,ebnd);
hold on
for i = 1: length(last_ring{inod})
  C = fix(p_nod{i}*70)+10;
  trisurf(tri_r_{i},x_r_{i}(:,1),x_r_{i}(:,2),p_nod{i},C,'CDataMapping','direct')
  
  atria = nn_prepare(x_r_{i});
  [~,nearest] = range_search(x_r_{i},atria,x_n(idnears,:),0.001);
  ind   = unique([nearest{:,1}]);
  %plot3(x_r_{i}(ind,1),x_r_{i}(ind,2),p_nod{i}(ind)+0.001,'ko','MarkerFaceColor','k','Markersize',4)
  p_max = max([p_max max(p_nod{i}(ind))]);
end
fprintf(1,'The max value of the basis functions is p_max=%f\n',p_max);
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'k-','LineWidth',1);
view([-20 25])
axis normal
camlight
shading interp
lighting phong
%material metal
material shiny
%colormap([0.7 0.7 0.7])
set(gca,'ZTickLabel',{'0','0.25','0.50','0.75'},'ZTick',[0 0.25 0.5 0.75],...
  'XTickLabel',{},'YTickLabel',{},...
  'FontSize',20,...
  'FontName','times');
xlabel('X')
ylabel('Y')
xlim([min(x_n(:,1)) max(x_n(:,1))])
ylim([min(x_n(:,2)) max(x_n(:,2))])
zlim([0 limP])
if savePlot==1
  print('-dpng','-r300',strcat('figures_png_poly/',fig_shpfun,'.png'))
  %savefig(strcat('figures_fig/',fig_shpfun,'.fig'))
end



%% -------------------------------------------------------------------
%Plot gradient in X of the basis functions at a given node
plot_domain(x_bd,x_n,ebnd);
hold on
for i = 1: length(last_ring{inod})
  C1 = fix((dp_nod{i}(:,1)+limDP)/(2*limDP)*70)+10;
  trisurf(tri_r_{i},x_r_{i}(:,1),x_r_{i}(:,2),dp_nod{i}(:,1),C1,'CDataMapping','direct')
end
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'k-','LineWidth',1);
view([-20 25])
axis normal
camlight
shading interp
set(gca,'ZTickLabel',{sprintf('-%d',limDP),'0',sprintf('-%d',limDP)},'ZTick',[-limDP 0 limDP],...
  'XTickLabel',{},'YTickLabel',{},...
  'FontSize',20,...
  'FontName','times');
xlim([min(x_n(:,1)) max(x_n(:,1))])
ylim([min(x_n(:,2)) max(x_n(:,2))])
zlim([-limDP limDP])
xlabel('X')
ylabel('Y')
title('Gradient X')
if savePlot==1
  print('-dpng','-r300',strcat('figures_png_poly/',fig_DshpfunX,'.png'))
  %savefig(strcat('figures_fig/',fig_DshpfunX,'.fig'))
end

%% -------------------------------------------------------------------
%Plot gradient in Y of the basis functions at a given node
plot_domain(x_bd,x_n,ebnd);
hold on
for i = 1: length(last_ring{inod})
  C2 = fix((dp_nod{i}(:,2)+limDP)/(2*limDP)*70)+10;
  trisurf(tri_r_{i},x_r_{i}(:,1),x_r_{i}(:,2),dp_nod{i}(:,2),C2,'CDataMapping','direct')
end
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'k-','LineWidth',1);
view([-20 25])
axis normal
camlight
shading interp
set(gca,'ZTickLabel',{sprintf('-%d',limDP),'0',sprintf('-%d',limDP)},'ZTick',[-limDP 0 limDP],...
  'XTickLabel',{},'YTickLabel',{},...
  'FontSize',20,...
  'FontName','times');
xlim([min(x_n(:,1)) max(x_n(:,1))])
ylim([min(x_n(:,2)) max(x_n(:,2))])
zlim([-limDP limDP])
xlabel('X')
ylabel('Y')
title('Gradient Y')
if savePlot==1
  print('-dpng','-r300',strcat('figures_png_poly/',fig_DshpfunY,'.png'))
  %savefig(strcat('figures_fig/',fig_DshpfunY,'.fig'))
end

if saveVTK == 1
  %% Generates a VTK file with the Shape function of the inod point
  
  nElem_r    = size(tri_r_{1},1);
  nElemTot_r = length(last_ring{inod})*nElem_r;
  tri_tot_r  = zeros(nElemTot_r,3);
  
  sPts_r     = size(x_r_{1},1);
  sPtsTot_r  = length(last_ring{inod})*sPts_r;
  x_tot_r    = zeros(sPtsTot_r,2);
  p_tot_r    = zeros(sPtsTot_r,1);
  dp_tot_r1  = zeros(sPtsTot_r,1);
  dp_tot_r2  = zeros(sPtsTot_r,1);
  
  for i = 1: length(last_ring{inod})
    epos_r = ((i-1)*nElem_r+1):i*nElem_r;
    tri_tot_r(epos_r,:) = tri_r_{i}+sPts_r*(i-1);
    
    xpos_r = ((i-1)*sPts_r+1):i*sPts_r;
    x_tot_r(xpos_r,:) = x_r_{i};
    
    p_tot_r(xpos_r) = p_nod{i};
    dp_tot_r1(xpos_r) = dp_nod{i}(:,1);
    dp_tot_r2(xpos_r) = dp_nod{i}(:,2);
  end
  
  
  plot_domain(x_bd,x_n);
  hold on
  trisurf(tri_tot_r,x_tot_r(:,1),x_tot_r(:,2),p_tot_r)
  view([-20 25])
  axis normal
  axis([min(x_n(:,1)) max(x_n(:,1)) min(x_n(:,2)) max(x_n(:,2)) ...
    -.1 1 ])
  camlight
  shading interp
  set(gca,'FontName','times')
  set(gca,'FontSize',20)
  
  %surf2vtk('data_vtk/elem_mesh.vtk',[x_n zeros(nPts,1)],tri,zeros(nPts,1))
  surf2vtk(strcat('data_vtk/',fig_shpfunct,'_ptA.vtk'),[x_tot_r  p_tot_r] ,tri_tot_r,p_tot_r)
  surf2vtk(strcat('data_vtk/',fig_shpfunDX,'_ptA.vtk'),[x_tot_r dp_tot_r1],tri_tot_r,dp_tot_r1)
  surf2vtk(strcat('data_vtk/',fig_shpfunDY,'_ptA.vtk'),[x_tot_r dp_tot_r2],tri_tot_r,dp_tot_r2)
end