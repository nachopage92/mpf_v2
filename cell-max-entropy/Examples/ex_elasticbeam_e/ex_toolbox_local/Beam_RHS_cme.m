function [rhs,ind_bcD] = Beam_RHS_cme(x_n,options,parameters)
% function [rhs,ind_Dirichlet] = Beam_RHS_cme(x_n,options,parameters);
%
% This function computes the right hand side rhs
%
% Input
%    x_n   : node points
%    parameters: L (length), D (diameter), nu (Poisson coefficient), E
%                (Young modulus)
%    options   : lme options
%
% Output
%    rhs    : right hand side
%    ind_bcD: position of Dirichlet's nodes


%% =======================================================================================
%  Input parameters: CME options
if isfield(options,'nodes_in_elem')
  nodes_in_elem = options.nodes_in_elem;
else
  error('Auxiliary arrays at each node indicating "nodes_in_elem" is not defined')
end
if isfield(options, 'tri')
  tri = options.tri;
else
  error('Input triangulation is needed')
end

%% =======================================================================================
%  CME: basis functions computation
nPts = size(x_n,1);    % Number of nodes
nDim = size(x_n,2);

% -------------------------------------------------------------------
% Searching Non-Delaunay Triangulations
tr     = triangulation(tri,x_n);
Cedges = edges(tr);
DTri   = delaunayTriangulation(x_n, Cedges);
dt2tr  = nan(size(DTri,1),1);
ic     = incenter(tr); %coordinates of the center of each triangle

dtid   = pointLocation(DTri,ic);
dt2tr(dtid) = 1:size(tr,1);


%% =======================================================================================

% Material parameters
P  = parameters.P;
D  = parameters.D;
L  = parameters.L;

I = D^3/12;

%% The BCs nodes are classified
ids   = (1:nPts);

%Nodes with imposed Dirichlet BCs -------------------------
%(x=0,y=0)
ind_bcD(1) = ids(abs(x_n(:,1))<1.e-6 & abs(x_n(:,2))<1.e-6);
%(x=0,y=D/2)
ind_bcD(2) = ids(abs(x_n(:,1))<1.e-6 & abs(x_n(:,2)-0.5*D)<1.e-6);

%Boundary with imposed loads-------------------------------

%Boundary normals in CCW (counter clock wise)
bd_n      = zeros(3,2);

% (x=0)
id_bnd{1} = ids(abs(x_n(:,1))<1.e-6);
bd_n(1,:) = [-1  0];  % (x=0)

% (x=L)
id_bnd{2} = ids(abs(x_n(:,1)-L)<1.e-6);
bd_n(2,:) = [ 1  0];  % (x=L)

% (y=0)
id_bnd{3} = ids(abs(x_n(:,2))<1.e-6);
bd_n(3,:) = [ 0 -1];  % (y=0)



%% =======================================================================================
% The rhs is initialized
rhs    = zeros(2*nPts,1);

% options for 1D integration and 1D cme
opt_int1D.orderGL = 10;
opt_int1D.pruning = 0;
options1D      = options;
options1D.dim  = 1;
options1D.grad = 0;

for ib=1:3
  n_  =   bd_n(ib,:)';
  id_ = id_bnd{ib};
  
  if abs(n_(1)) > 0
    i1 = 2;
    i2 = 1;
  else
    i1 = 1;
    i2 = 2;
  end
  
  [eta_n,id_sort] = sort(x_n(id_,i1));
  
  id_glob = id_(id_sort);
  x_bd    = x_n(id_glob,:);
  
  [eta_s,wth_s] = MakeGLSamples1D(eta_n, opt_int1D);
  
  sPts= length(eta_s);
  x_s = zeros(sPts,2);
  x_s(:,i1) = eta_s;
  
  %x_s(:,i2) = (mean(x_bd(:,i2)) - n_(i2)*1e-10)*ones(sPts,1);
  x_s(:,i2) = mean(x_bd(:,i2))*ones(sPts,1);
    
  figure(20);clf
  plot(x_n(:,1),x_n(:,2),'ro')
  hold on
  plot(x_n(id_,1),x_n(id_,2),'r*','markersize',12,'linewidth',2)
  plot(x_s(:,1),x_s(:,2),'bx','markersize',14)
  axis equal
  
  % -------------------------------------------------------------------
  % returns the index of the enclosing Delaunay triangle for each point in x_s
  ja_tri = pointLocation(DTri,x_s);
  
  ja_tri = dt2tr(ja_tri);
  
  [ja_tri,id_loc] = sort(ja_tri);
  ja_bin = unique(ja_tri);
  ia_tri = histc(ja_tri,ja_bin);
  
  ia_tri = cumsum(ia_tri);   %cummulative number of sample points in elements
  ja_tri = ja_tri(ia_tri);   %global element index containing the sample points
  
  sElem  = length(ja_tri);   %number of elements containing sample points
  
  
  %x_s(:,i2) = mean(x_bd(:,i2))*ones(sPts,1);
  
  % -------------------------------------------------------------------
  % Compute basis functions in each element containing sample points
  
  sInit = 0;
  for e = 1: sElem
    iel     = ja_tri(e);              %global element index
    
    enear   = nodes_in_elem{iel};     %nearest nodes to the e-th element
    
    options.enear = enear;
    
    % rescaling ---------------------------------------------
    xe_n  = x_n(enear,:);
    
    %length of the bounding box containing the nearest nodes
    bbox = 0;
    for i=1:nDim
      dx   = abs(max(xe_n(:,i))-min(xe_n(:,i)));
      bbox = bbox + dx*dx;
    end
    bbox    = sqrt(bbox);
    scale   = 1/bbox/4;
    
    x_n_    = x_n*scale;
    options.h = options.node_spacing*scale;
    
    kPts    = ia_tri(e)-sInit;        %number of sample points in the e-th triangle
    id_samp = (sInit+1):(sInit+kPts); %IDs of the sample points belonging to the e-th triangle
    id_samp = id_loc(id_samp);
    xe_s_   = x_s(id_samp,:)*scale;   %sample points coordinates
    
%     figure(21);clf
%     plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','k','Markersize',10 )
%     hold on
%     plot(x_bd(:,1),x_bd(:,2),'ko','MarkerFaceColor','c','Markersize',8 )
%     triplot(tri,x_n(:,1),x_n(:,2),'m-','LineWidth',1);
%     plot(x_s(:,1),x_s(:,2),'bx','MarkerFaceColor','b','Markersize',10 )
%     plot(x_n(tri(iel,:),1),x_n(tri(iel,:),2),'ks','MarkerFaceColor','g','Markersize',14,'linewidth',1)
%     plot(x_s(id_samp,1),x_s(id_samp,2),'mx','MarkerFaceColor','b','Markersize',14,'linewidth',2)
%     hold off
%     axis equal
%     pause
    
    
    % prior Rfunctions --------------------------------------
    [q,dq]=priorf_elem(x_n_,xe_s_,options);
    
    %[C,IA,IB] = intersect(A,B,'stable')
    %C = A(IA) and C = B(IB);
    [enear_,Inear] = intersect(enear,id_glob,'stable');
    
    % cell-based max-ent functions --------------------------
    options.q  = q(:,Inear);
    options.dq = dq(:,Inear,i1);
    options.enear = enear_;
    p = shapef_elem(x_n_(:,i1),xe_s_(:,i1),options);
    
    
    % Right hand side
    rhs_x = zeros(length(enear_),1);
    rhs_y = zeros(length(enear_),1);
    for k=1:kPts
      s    = id_samp(k);
      
      p_k  = p(k,:)';
      w_k  = wth_s(s);
      
      % Boundary nodes corresponding to x=0  #-------------
      if ib==1
        y_k  = x_s(s,2);
        rhs_x = rhs_x + P*L/I*y_k*p_k*w_k;
        rhs_y = rhs_y - P/(2*I)*(D^2/4-y_k^2)*p_k*w_k;
        
      % Boundary nodes corresponding to x=L  #-------------
      elseif ib==2
        y_k  = x_s(s,2);
        rhs_y = rhs_y + P/(2*I)*(D^2/4-y_k^2)*p_k*w_k;
        
      % Boundary nodes corresponding to y=0  #-------------
      else
        rhs_x = rhs_x - P/(2*I)*D^2/4*p_k*w_k;
      end
      
    end
    rhs(2*enear_-1) = rhs(2*enear_-1) + rhs_x;
    rhs(2*enear_)   = rhs(2*enear_)   + rhs_y;
    
    sInit = ia_tri(e);
  end
end
