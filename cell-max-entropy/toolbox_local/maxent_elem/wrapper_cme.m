function output = wrapper_cme(x_n,x_s,options)
% function output = wrapper_cme(x_n,x_s,options)
% INPUT:
% x_n:      nodal set (nPts, nDim)
% x_s  :    sample point where the shape functions are evaluated (in an elcment)
%
% options: parameters setting
%   options.bd_segm    : segments defining the polygon of the last ring, which are not on the boundary 
%   options.bd_segm_out: segments defining the polygon of the last ring including segments on the boundary 
%   options.bd_nodes   : boundary identifier (0:false, 1:true)
%   options.spow       : w^spow, where w is the approximation to the distance function (R-function)
%   options.mpow       : distance function derivative maximum degree of the approximation
%   options.isConv     : flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)
%   options.TolNR      : tolerance for Newton-Raphson's iterations
%   nodes_in_elem      : nearest nodes to the e-th element
%
% OUTPUT:
%   p_cme  : are the shape functions values for each sample point
%   dp_cme : gradient for the shape functions at the sample points
%   s_near : list of nearest neighbors for each sample point

%% Basic input checkings and variable definitions ========================================
if nargout > 1
  error('The number of output is too big');
end
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


sPts   = size(x_s,1);
nDim   = size(x_s,2);

% ========================================================================================
%% Searching Non-Delaunay Triangulations
tr     = triangulation(tri,x_n);
Cedges = edges(tr);
DTri   = delaunayTriangulation(x_n, Cedges);
dt2tr  = nan(size(DTri,1),1);
ic     = incenter(tr);

dtid = pointLocation(DTri,ic);
dt2tr(dtid) = 1:size(tr,1);


% ========================================================================================
% returns the index of the enclosing Delaunay triangle for each point in x_s
ja_tri = pointLocation(DTri,x_s);

ja_tri = dt2tr(ja_tri);

[ja_tri,id_loc] = sort(ja_tri);
ja_bin = unique(ja_tri);
ia_tri = histc(ja_tri,ja_bin);

ia_tri = cumsum(ia_tri);   %cummulative number of sample points in elcments
ja_tri = ja_tri(ia_tri);   %global elcment index containing the sample points

sElem  = length(ja_tri);   %number of elcments containing sample points

% ========================================================================================
%% Compute basis functions in each elcment containing sample points
p_cme  = {zeros(sPts,1)};
dp_cme = {zeros(sPts,1)};

s_near = {zeros(sPts,1)};

sInit = 0;
for e = 1: sElem
  iel     = ja_tri(e);              %global element index
  enear   = nodes_in_elem{iel};     %nearest nodes to the e-th element
   
  options.enear = enear;
  
  % rescaling -------------------------------------------------------------
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
  
  % prior Rfunctions ------------------------------------------------------
  [q,dq]=priorf_elem(x_n_,xe_s_,options);
  options.q  = q;
  options.dq = dq;
  
  % cell-based max-ent functions ------------------------------------------
  [p,dp] = shapef_elem(x_n_,xe_s_,options);
 
  for j=1:kPts
    s = id_samp(j);  
       
    p_cme{s} = p(j,:)';
    dp_cme{s}= ([dp(j,:,1);dp(j,:,2)]')*scale;  %recover rescale
    s_near{s}= enear;
  end
  sInit = ia_tri(e);
end

output.p_samp  = p_cme;
output.dp_samp = dp_cme;
output.s_near  = s_near;
