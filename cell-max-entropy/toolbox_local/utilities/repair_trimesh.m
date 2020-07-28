function [tri,x_n] = repair_trimesh(tri,x_n)
% function [tri,x_n] = repair_trimesh(tri,x_n)

% -------------------------------------------------------------------
% Boundary nodes identifiers
ebnd   = freeBoundary(triangulation(tri,x_n));
ibnd   = unique(ebnd);
bPts   = length(ibnd);
x_bd   = x_n(ibnd,:);

nElem  = size(tri,1);
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

%% --------------------------------------------------------------------
id_elem = (1:nElem)';

el_bad = zeros(nElem,1);
for i=1:bPts
  if sum(any(tri(:,:)==ibnd(i),2))>0
    el_bad = el_bad + any(tri(:,:)==ibnd(i),2);
  end
end
el_bad = id_elem(el_bad==3);
id_bad = unique(tri(el_bad,:));

el_near = zeros(length(el_bad),1);
for e=1:length(el_bad)
  ee_near = zeros(nElem,1);
  id_tri  = tri(el_bad(e),:);
  for i=1:3
    if sum(any(tri(:,:)==id_tri(i),2))>0
      ee_near = ee_near + any(tri(:,:)==id_tri(i),2);
    end
  end
  el_near(e) = setdiff(id_elem(ee_near==2),el_bad(e));
end
id_near = unique(tri(el_near,:));

id_newA=[];
id_newB=[];
for e=1:length(el_bad)
  id_badA  = tri(el_bad(e),:);
  id_badB  = tri(el_near(e),:);
  id_union = union(id_badA,id_badB);
  id_inter = intersect(id_badA,id_badB);
  id_newA  = [setdiff(id_union,id_inter) id_inter(1)];
  id_newB  = [setdiff(id_union,id_inter) id_inter(2)];
  tri(el_bad(e),:) = id_newA;
  tri(el_near(e),:)= id_newB;
end

%% Check the ordering
el_change = [el_bad el_near];
for i=1:length(el_change)
  e   = el_change(i);
  x_e = x_n(tri(e,:),:);
  x21=x_e(2,:)-x_e(1,:);
  x32=x_e(3,:)-x_e(2,:);
  v_e = cross([x21 0],[x32 0]);
  if v_e(3)<0
    tri(e,:) = [tri(e,1) tri(e,3) tri(e,2)];
  end
end

%Check all triangles
for e=1:nElem
  x_e = x_n(tri(e,:),:);
  x21=x_e(2,:)-x_e(1,:);
  x32=x_e(3,:)-x_e(2,:);
  v_e = cross([x21 0],[x32 0]);
  if v_e(3)<0
    tri(e,:) = [tri(e,1) tri(e,3) tri(e,2)];
  end
end

iter = 0;
iend = 0;
while iend == 0
  nbad   = 0;
  el_bad = zeros(nElem,1);
  for e=1:nElem
    x_e = x_n(tri(e,:),:);
    x21=x_e(2,:)-x_e(1,:);
    x32=x_e(3,:)-x_e(2,:);
    v_e = cross([x21 0],[x32 0]);
    if v_e(3)<0
      tri(e,:)  = [tri(e,1) tri(e,3) tri(e,2)];
      el_bad(e) = 1;
      nbad = nbad+1;
    end
  end
  
  if nbad == 0
    iend = 1;
  end
  if iter > 10
    id_el = 1:nElem;
    disp(id_el(el_bad==1))
    error('Wrong ordered elements') 
  end
  iter = iter + 1;
end

return
%% --------------------------------------------------------------------
nPts   = size(x_n,1);
C1     = zeros(nPts, 1);

figure(1);clf
trimesh(tri,x_n(:,1),x_n(:,2),C1,...
  'Facecolor','none','FaceLighting','none','EdgeColor','k','EdgeLighting','flat')
hold on
plot(x_n(:,1),x_n(:,2), ...
  'ko','MarkerFaceColor','k','Markersize',4 )
plot(x_n(id_bad,1),x_n(id_bad,2), 'r^','Markersize',12,'linewidth',2)
plot(x_n(id_near,1),x_n(id_near,2), 'bv','Markersize',12,'linewidth',2)
hold off
view([0,90])
axis equal


figure(2);clf
trimesh(tri,x_n(:,1),x_n(:,2),C1,...
  'Facecolor','none','FaceLighting','none','EdgeColor','k','EdgeLighting','flat')
hold on
plot(x_n(:,1),x_n(:,2), ...
  'ko','MarkerFaceColor','k','Markersize',4 )
plot(x_n(id_newA,1),x_n(id_newA,2), 'g<','Markersize',12,'linewidth',2)
plot(x_n(id_newB,1),x_n(id_newB,2), 'c>','Markersize',12,'linewidth',2)
hold off
view([0,90])
axis equal


%We define a non-convex domain with alpha shapes
figure(3);clf
hold on
plot(x_n(:,1),x_n(:,2), ...
  'ko','MarkerFaceColor','b','Markersize',10 )
hold on

%We define a non-convex domain
%ebnd = freeBoundary(triangulation(tri,x_n));
for k = 1:size(ebnd,1)
  plot(x_n(ebnd(k,:),1),x_n(ebnd(k,:),2),'r-','LineWidth',4)
end
plot(x_bd(:,1),x_bd(:,2), ...
  'ko','MarkerFaceColor','c','Markersize',8 )

%Triangulation constrained by the boundary edges
triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'m-','LineWidth',1);
axis equal

