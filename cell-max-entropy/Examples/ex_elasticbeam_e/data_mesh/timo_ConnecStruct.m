clear all
close all

fprintf(1,'====================================================================\n');
  
%% Matrial Parameters
L  = 4;         % length
D  = 1;         % height
Lx = L;
Ly = 0.5 * D;

makePlot     = 1;
optCME.dim   = 2;
optCME.verb  = 0; % 0:off    1:on
optCME.nring = 2;

%% Convervence curves
NChen  = [3 5 9 17 33 65 129];

for n=1:length(NChen)
  
  % Node Points
  Ny = NChen(n);
  Nx = NChen(n)*8-3;
  hx = Lx/(Nx-1);
  hy = Ly/(Ny-1);
  fprintf(1,'Uniform grid %dx%d ---- h=[%5.3f %5.3f]===============\n', Nx, Ny, hx, hy);

  % This line is important
  center= 0.5*[Lx Ly];
  x_n   = UniformGrid2D(Nx, Ny, Lx, Ly, center);
  nPts  = size(x_n,1);
  
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

  fileOUT = sprintf('timo_mesh_N%03d_ring%d.mat', Ny, optCME.nring);
  %save(fileOUT,'tri','x_n','ibnd','connecStr')
  
  fprintf(1,'\tcputime CME Connectivity Structures: %5.3f seconds\n', cputime-t);
  pause
end
