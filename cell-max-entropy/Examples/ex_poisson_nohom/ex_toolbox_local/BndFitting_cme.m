function u_n = BndFitting_cme(x_n,id_top,options)
% function u_n = BndFitting_cme(x_n,id_top,options)
%
% INPUT:
%     x_nodes: points in 1D
%
% OUTPUT:
%     u_nodes: value at the nodes/control points

%% =======================================================================================
%  Input parameters: LME options
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

nPts = length(id_top);
L0   = min(x_n(id_top,1));
L1   = max(x_n(id_top,1));

nDim = size(x_n,2);
% -------------------------------------------------------------------
x_s  = zeros(10*nPts,nDim);
sPts = size(x_s,1);  %total number of sample points

x_s(:,1) = linspace(L0,L1,sPts)';
x_s(:,2) = ones(sPts,1)*mean(x_n(id_top,nDim));



%% =======================================================================================
%  CME: basis functions computation

% -------------------------------------------------------------------
% Searching Non-Delaunay Triangulations
tr     = triangulation(tri,x_n);
Cedges = edges(tr);
DTri   = delaunayTriangulation(x_n, Cedges);
dt2tr  = nan(size(DTri,1),1);
ic     = incenter(tr);

dtid   = pointLocation(DTri,ic);
dt2tr(dtid) = 1:size(tr,1);


% -------------------------------------------------------------------
% returns the index of the enclosing Delaunay triangle for each point in x_s
ja_tri = pointLocation(DTri,x_s);

ja_tri = dt2tr(ja_tri);

[ja_tri,id_loc] = sort(ja_tri);
ja_bin = unique(ja_tri);
ia_tri = histc(ja_tri,ja_bin);

ia_tri = cumsum(ia_tri);   %cummulative number of sample points in elcments
ja_tri = ja_tri(ia_tri);   %global elcment index containing the sample points

sElem  = length(ja_tri);   %number of elcments containing sample points

% -------------------------------------------------------------------
% Compute basis functions in each elcment containing sample points
p_cme  = {zeros(sPts,1)};

s_near = {zeros(sPts,1)};

sInit = 0;
for e = 1: sElem
  iel     = ja_tri(e);              %global elcment index
  enear   = nodes_in_elem{iel};     %nearest nodes to the e-th elcment
   
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
  
  kPts    = ia_tri(e)-sInit;        %number of sample points in the e-th triangle
  id_samp = (sInit+1):(sInit+kPts); %IDs of the sample points belonging to the e-th triangle
  id_samp = id_loc(id_samp);
  xe_s_   = x_s(id_samp,:)*scale;   %sample points coordinates
  
  % prior Rfunctions --------------------------------------
  [q,dq]=priorf_elem(x_n_,xe_s_,options);

  %[C,IA,IB] = intersect(A,B,'stable')
  %C = A(IA) and C = B(IB);
  [enear_top,Inear,Itop] = intersect(enear,id_top,'stable');
  
  % cell-based max-ent functions --------------------------
  options.q  = q(:,Inear);
  options.dq = dq(:,Inear,1);
  options.enear = enear_top;
  p = shapef_elem(x_n_(:,1),xe_s_(:,1),options);
 
  for j=1:kPts
    s = id_samp(j);  
       
    p_cme{s} = p(j,:)';
    s_near{s}= Itop;
  end
  sInit = ia_tri(e);
end




%% =======================================================================================
% Least squares system assembly 
u_s  = x_s(:,1).*(1-x_s(:,1));

nn=0;
for k=1:sPts
  nn = max(nn, length(p_cme{k}));
end
nn = min(nn,nPts);

M = spalloc(nPts,nPts,nn*nPts);
f = zeros(nPts,1);

for k=1:sPts
  k_near = s_near{k};
  p_k    = p_cme{k};
  
  M(k_near,k_near) = M(k_near,k_near) + (p_k*p_k');
  f(k_near) = f(k_near) + p_k*u_s(k);
end

u_n = M\f;

u_sh = zeros(sPts,1);
for k=1:sPts
  k_near = s_near{k};
  p_k    = p_cme{k};
  
  u_sh(k) = p_k'*u_n(k_near);
end

figure(10)
plot(x_s(:,1),u_s,'bx-','linewidth',2)
hold on
plot(x_s(:,1),u_sh,'k+-','linewidth',2)
plot(x_n(id_top,1),u_n,'ro','markersize',12)
hold off
fprintf(1,'\t|u - u_h| in L2 is %8.2e\n',norm(u_s-u_sh));
pause(0.1)