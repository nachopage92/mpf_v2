%% Purpose: Connectivity Structures Generator - cell-based max-ent basis functions

clear all
close all

format short

addpath('ex_toolbox_local')

disp('------------------------------------------------')

%% =======================================================================================
%  Input parameters: CME options
makePlot     = 0;
optCME.dim   = 2;
optCME.verb  = 0; % 0:off    1:on
optCME.nring = 2;

% -------------------------------------------------------------------
% Geometrical Parameters
L     = 1.0;    % length
Ndof  = [3 5 9 17 33 65 129];

for n=1:length(Ndof)
  disp('------------------------------------------------')
  
  % Node Points
  N    = Ndof(n);
  x_n  = SquareNodes(0, N, L);
  nPts = size(x_n,1);
  fprintf(1,'\tGrid size: %2dx%2d\n',N,N);
  
  % -------------------------------------------------------------------
  % Gauss Legendre quadrature points for numerical integration
  tri = delaunay(x_n(:,1),x_n(:,2));
   
  % -------------------------------------------------------------------
  % Boundary nodes identifiers
  ebnd   = freeBoundary(triangulation(tri,x_n));
  ibnd   = unique(ebnd);
  bPts   = length(ibnd);
  x_bd   = x_n(ibnd,:);
   
  nElem = size(tri,1);        %number of elements
  
  %checking wrong elements with zero area and deleting them
  badelem = zeros(nElem,1);
  for e=1:nElem
    if TriangleArea(x_n(tri(e,:),:)) < 1.0e-06
      badelem(e) = 1;
    end
  end
  if sum(badelem)>0
    fprintf(1,'Bad elements with area < 1e-06:  %d\n',sum(badelem));
  end
  tri(badelem==1,:)  = [];
  nElem  = size(tri,1);
  inside = ones(nElem,1); 
  
  
  % -------------------------------------------------------------------
  % Rapair triangles which have three nodes on the boundary and order all
  [tri,x_n] = repair_trimesh(tri,x_n);
  
  %plot triangulation constrained by the boundary edges
  if makePlot==1
    figure(1);clf
    hold on
    plot(x_n(:,1),x_n(:,2), ...
      'ko','MarkerFaceColor','b','Markersize',10 )
    hold on
    for k = 1:size(ebnd,1)
      plot(x_n(ebnd(k,:),1),x_n(ebnd(k,:),2),'r-','LineWidth',4)
    end
    plot(x_bd(:,1),x_bd(:,2), ...
      'ko','MarkerFaceColor','c','Markersize',8 )
    triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'m-','LineWidth',1);
    axis equal
  end  
  
  %% =======================================================================================
  %  EME: Connectivity structures
  
  t=cputime;
  optCME.tri  = tri;
  optCME.ibnd = ibnd;
  connecStr   = ConnectivityStructures_elem(x_n,optCME);
  
  bd_segm = connecStr.bd_segm;
  for i=1:nPts
    if sum(isempty(bd_segm{i}))>0 && makePlot==1
      figure(2);clf
      hold on
      plot(x_n(:,1),x_n(:,2), ...
        'ko','MarkerFaceColor','b','Markersize',10 )
      plot(x_n(i,1),x_n(i,2),'ko','MarkerFaceColor','c','Markersize',8 )
      triplot(tri(inside==1, :),x_n(:,1),x_n(:,2),'m-','LineWidth',1);
      axis equal
      pause
    end
  end
  
  fileOUT = sprintf('data_mesh/square_mesh_N%03d_ring%d.mat', N, optCME.nring);
  %save(fileOUT,'tri','x_n','ibnd','connecStr')
  
  fprintf(1,'\tcputime CME Connectivity Structures: %5.3f seconds\n', cputime-t);
end
