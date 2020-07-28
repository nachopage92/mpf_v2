function output = ConnectivityStructures_elem(x_n,options)
% function output = ConnectivityStructures_elem(x_n,options)
%
% INPUT:
%   x_n    : node coordinates
%
%   options: parameters setting
%     options.nring : the last ring of neighbors to consider for each node: 1, 2, 3.
%     options.tri   : Triangulation
%     options.ibnd  : nodes's ID indicating the nodes on the boundary of the domain/mesh
%
% OUTPUT:
%   output.nodes_in_elem : map to each element, the real nodes affecting it
%   output.bd_segm       : segments defining the polygon of the last ring, which are not on the boundary
%   output.bd_segm_out   : segments defining the polygon of the last ring including segments on the boundary
%   output.bd_nodes      : boundary identifier (0:false, 1:true)
%   output.first_ring    : nodes's ID up to the first ring of neighbors
%   output.last_ring     : nodes's ID up to the last ring of neighbors

% ========================================================================================
%% Input parameters
if isfield(options, 'tri')
  tri = options.tri;
else
  error('Input DelaunayTriangulation is needed')
end

if isfield(options, 'nring')
  nring = options.nring;
else
  nring = 2;
end


if isfield(options, 'polyCase')
  polyCase = options.polyCase;
else
  polyCase = 'C';
end

if nring>3
  error('The maximum of rings is 3')
end

if polyCase == 'C'
  ibnd = options.ibnd;
  ebnd = [];
else
  ebnd = freeBoundary(triangulation(tri,x_n));
  ibnd = unique(ebnd);
end


if isfield(options, 'cleanColinear')
  cleanColinear = options.cleanColinear;
else
  cleanColinear = 1;
end


nPts   = size(x_n,1);

% ========================================================================================
%% Create the connectivity structures
nElem  = size(tri,1);
% inside = ones(nElem,1);
% count  = 1:nElem;

%% ############################ 1-RING ##################################
%first ring of elements around one node
first_ring = {[]};
ncount = zeros(nPts,1);
for e=1:nElem
  ncount(tri(e,:)) = ncount(tri(e,:)) + 1;
end
for i=1:nPts
  first_ring{i} = zeros(1,ncount(i));
  ncount(i) = 1;
end
for e=1:nElem
  ering = tri(e,:); %elements in the second ring of the node j
  for k=1:length(ering)
    i    = ering(k);
    ipos = ncount(i);
    first_ring{i}(ipos) = e;
    ncount(i) = ipos + 1;
  end
end

if nring > 1
  %% ############################ 2-RING ##################################
  %second ring of elements around one node
  second_ring = {[]};
  for i = 1:nPts
    neighbornodes = unique(tri(first_ring{i},:));
    second_ring{i}= unique([first_ring{neighbornodes}]);
  end
  
  if nring == 2
    last_ring=second_ring;
  else
    %% ########################## 3-RING ####################################
    %third ring of elements around one node
    third_ring = {[]};
    for i = 1:nPts
      neighbornodes = unique(tri(second_ring{i},:));
      third_ring{i} = unique([first_ring{neighbornodes}]);
    end
    last_ring = third_ring;
  end
else
  last_ring = first_ring;
end

% %--------------------------------------------------------------------
% %first ring of elements around one node
% for i = 1:nPts
%   first_ring(i) = {count(any(tri(:,:)==i,2))};
%   first_ring(i) = {first_ring{i}(inside(first_ring{i})==1)};
% end
%
% if nring > 1
%   %--------------------------------------------------------------------
%   %second ring of elements around one node
%   bd_nodes = zeros(nPts,1);
%   for i = 1:nPts
%     neighbornodes = unique(tri(first_ring{i},:));
%     list = [];
%     for j = 1:length(neighbornodes)
%       list = [list first_ring{neighbornodes(j)}];
%     end
%     second_ring(i)= {unique(list)};
%   end
%   if nring == 2
%     last_ring=second_ring;
%   else
%     %--------------------------------------------------------------------
%     %third ring of elements around one node
%     for i = 1:nPts
%       neighbornodes = unique(tri(second_ring{i},:));
%       list = [];
%       for j = 1:length(neighbornodes)
%         list = [list first_ring{neighbornodes(j)}];
%       end
%       third_ring(i)= {unique(list)};
%     end
%     last_ring=third_ring;
%   end
% else
%   last_ring=first_ring;
% end
%

%--------------------------------------------------------------------
% map to each element, the real nodes affecting it
nodes_in_elem = {[]};

ecount = zeros(nElem,1);
for j=1:nPts
  ecount(last_ring{j}) = ecount(last_ring{j}) + 1;
end
for e=1:nElem
  nodes_in_elem{e} = zeros(1,ecount(e));
  ecount(e) = 1;
end
for j=1:nPts
  jring = last_ring{j}; %elements in the third ring of the node j
  for k=1:length(jring)
    e    = jring(k);
    jpos = ecount(e);
    nodes_in_elem{e}(jpos) = j;
    ecount(e) = jpos + 1;
  end
end

% %--------------------------------------------------------------------
% % map to each element, the real nodes affecting it
% for i=1:nElem
%   list = [];
%   for j=1:nPts
%     if any(last_ring{j}==i)
%       list = [list j];
%     end
%   end
%   if inside(i)==1
%     nodes_in_elem(i) = {list};
%   else
%     nodes_in_elem(i) = {[]};
%   end
% end


%% --------------------------------------------------------------------
% %% Vectors at the boundary
% nEdg = size(ebnd,1);
% vecEdg = zeros(nEdg,2);
% for ed=1:nEdg
%   vecEdg(ed,:) = x_n(ebnd(ed,2),:) - x_n(ebnd(ed,1),:);
% end

node_first = zeros(nPts,1);
for i = 1:length(ibnd)
  el_list = first_ring{ibnd(i)};
  id_list = unique(tri(el_list,:));
  node_first(id_list) = 1;
end
% node_first(ibnd) = 0;


bd_nodes = zeros(nPts,1);
bd_segm_out={[]};
bd_segm    ={[]};
for i = 1:nPts
  if nring > 1
    %put togheter the nodes in the "last" ring of the actual node I
    nodes2 = unique(tri(last_ring{i},:));
    
    if nring == 2
      %put togheter the nodes in the first ring of the actual node I
      nodes1 = unique(tri(first_ring{i},:));
    else
      %put togheter the nodes in the second ring of the actual node I
      nodes1 = unique(tri(second_ring{i},:));
    end
    %remove the first/second ring nodes to the "last" one
    tmp    = setdiff(nodes2,nodes1);
    
    %remove also the actual node I
    tmp    = setdiff(tmp,i);
    %find the intersection of the "last" ring nodes with nodes in the boundary
    %note: this could be empty or will contain the actual node I and some nodes in the first ring
    tmp_   = intersect(nodes2,ibnd);
    %put all togheter: the nodes in the "last" ring without the first ring of nodes plus
    %the boundary nodes of the first ring (can contain the actual node I)
    bd_sec_ring = union(tmp,tmp_);
  else
    %put togheter the nodes in the first ring of the actual node I
    nodes1 = unique(tri(first_ring{i},:)); 
    %remove the actual node I
    tmp    = setdiff(nodes1,i);
    %find the intersection of the first ring with the nodes in the boundary
    %note: this could be empty or will contain the actual node I and some nodes in the first ring
    tmp_   = intersect(nodes1,ibnd);
    %put all togheter: the nodes in the first ring without the node I plus the boundary
    %nodes (can contain the actual node I)
    bd_sec_ring = union(tmp,tmp_);
  end
  
  %%
  % sort polygon; ONLY WORKS FOR TRIANGULAR ELEMENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  poly = bd_sec_ring(1);
  for ii = 1:length(bd_sec_ring)-1
    node1= poly(end);
    aux  = any(tri(last_ring{i},:)==node1,2); % elements containing node1
    aux2 = intersect(unique(tri(last_ring{i}(aux),:)),bd_sec_ring);
    aux2 = setdiff(aux2,poly); %potential nodes of the boundary connecting to node1
    
    %HERE IS THE PROBLEM! IT SHOULD BE MODIFIED TO CONSIDER CONNECTED EDGES/FACES
    for j = 1:length(aux2)
      aux_ = any(tri(last_ring{i},:)==aux2(j),2); % elements containing aux2(j)
      if sum(aux&aux_) == 1 % boundary edge
        poly=[poly aux2(j)];
        break
      end
    end
%     figure(2);clf
%     plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','b','Markersize',6 )
%     hold on
%     trimesh(tri,x_n(:,1),x_n(:,2),zeros(nPts,1),'Facecolor','none','FaceLighting','none','EdgeColor','b','EdgeLighting','flat')
%     triplot(tri(last_ring{i}, :),x_n(:,1),x_n(:,2),'k-','LineWidth',3);
% %    quadmesh(tri,x_n(:,1),x_n(:,2),C1,'Facecolor','none','FaceLighting','none','EdgeColor','k','EdgeLighting','flat')
% %    quadplot(tri(last_ring{i}, :),x_n(:,1),x_n(:,2),'k-','LineWidth',2);
%     
%     plot(x_n(neighbornodes,1),x_n(neighbornodes,2),'go','Markersize',12 ,'LineWidth',2);
%     plot(x_n(bd_sec_ring,1),x_n(bd_sec_ring,2),'bs','Markersize',14 ,'LineWidth',2);
%     plot(x_n(node1,1),x_n(node1,2),'rs','Markersize',14 ,'LineWidth',2);
%     plot(x_n(i,1),x_n(i,2), 'k*','Markersize',20 ,'LineWidth',3);
%     plot(x_n(poly,1),x_n(poly,2), 'mv','Markersize',20 ,'LineWidth',3);
%     axis equal
%     pause
  end
  
  
  %%
  if norm(union(bd_sec_ring, poly')-bd_sec_ring)~=0
    disp('Problem 0')
    return
  end
  aux= mod(1:length(poly),length(poly))+1;
  
  % Set the segments at the boundary of the last ring
  segm = [poly;poly(aux)]';
  
  %------------------------------------------------------------------
  % merging colinear edges
  if cleanColinear == 1
    segm_out = clean_colinearOUT(x_n,segm);
  else
    segm_out = segm;
  end
  bd_segm_out{i} = segm_out;
  
  %------------------------------------------------------------------
  % Remove elements in boundary of domain: only first ring nodes
  id_list = unique(tri(first_ring{i},:));
  segm = clean_colinearIN(polyCase,i,x_n,ibnd,ebnd,id_list,segm_out,segm);
  
  bd_segm{i} = segm;  
  
  %------------------------------------------------------------------
  if any(i==ibnd)
    bd_nodes(i) = 1;
  end
end


output.nodes_in_elem = nodes_in_elem;
output.bd_segm       = bd_segm;
output.bd_segm_out   = bd_segm_out;
output.bd_nodes      = bd_nodes;
output.first_ring    = first_ring;
output.last_ring     = last_ring;







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge colinear segments
function segm = clean_colinearOUT(x_n,segm)

iter = 0;
oend = 0;
while oend==0
  nSgm = size(segm,1);
  
  segm_new = segm;
  ek  = 1;
  ed  = 1;
  iend= 0;
  while iend==0
    v1 = x_n(segm(ed,2),:) - x_n(segm(ed,1),:);
    v2 = x_n(segm(ed+1,2),:) - x_n(segm(ed+1,1),:);
    
    if(norm(v1)) < eps
      error('v1 have zero length')
    end
    if(norm(v2)) < eps
      error('v2 have zero length')
    end
    
    ie1 = segm(ed,1);
    ie2 = segm(ed,2);
    
    if acos(v1*v2'/(norm(v1)*norm(v2))) < 0.05
      ie2 = segm(ed+1,2);
      ed = ed + 1;
    end
    segm_new(ek,1) = ie1;
    segm_new(ek,2) = ie2;
    ed = ed + 1;
    ek = ek + 1;
    
    if ed > nSgm-1
      if ie2 ~= segm_new(1,1)
        segm_new(ek,1) = ie2;
        segm_new(ek,2) = segm_new(1,1);
      else
        ek = ek - 1;
      end
      segm_new(ek+1:end,:) = [];
      iend=1;
    end
  end
  
  if size(segm,1) == size(segm_new,1) && iter > 5
    oend = 1;
  else
    nSgm_new = size(segm_new,1);
    id_ord   = [2:nSgm_new 1];
    segm_new = segm_new(id_ord,:);
  end
  iter = iter + 1;
  segm = segm_new;
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segm = clean_colinearIN(polyCase,i,x_n,ibnd,ebnd,id_list,segm_out,segm)
% [A] ---------------------------------------------------------------
id = zeros(size(ebnd,1),1);
for kk=1:length(id_list)
  id = id + any(ebnd(:,:)==id_list(kk),2);
end
ebnd_i = ebnd(id>0,:);

if polyCase == 'A';
  segm = segmIN_caseA(i,x_n,ebnd_i,ibnd,id_list,segm);
elseif polyCase == 'B';
  segm = segmIN_caseB(i,x_n,ebnd_i,ibnd,id_list,segm);
elseif polyCase == 'C';
  segm = segmIN_caseC(ibnd,segm);
else
  error('Polygon case is not defined')
end

% [B] Performed to avoid many colinear edges ------------------------
to_glue = zeros(size(segm_out,1),1);
for ek=1:size(segm_out,1)
  %distance from the node i to the segment out (infinite line)
  x1 = x_n(segm_out(ek,1),:);
  x2 = x_n(segm_out(ek,2),:);
  s  = x2 - x1;
  ds = norm(s);
  
  for ii=1:size(segm,1)
    s1 = x_n(segm(ii,1),:) - x1;
    da = abs(s1(1)*s(2) - s1(2)*s(1))/ds;
    
    s1 = x_n(segm(ii,2),:) - x1;
    db = abs(s1(1)*s(2) - s1(2)*s(1))/ds;
    
    if da/ds < 0.05 && db/ds < 0.05
      to_glue(ek)=1;
      break;
    end
  end
end
segm = segm_out(to_glue==1,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segm = segmIN_caseA(i,x_n,ebnd_i,ibnd,id_list,segm)
to_remove = zeros(size(segm,1),1);

for ii=1:size(segm,1)
  
  if any(segm(ii,1)==ibnd) && any(segm(ii,2)==ibnd)
    
    %distance from the node i to the segment (infinite line)
    x1 = x_n(segm(ii,1),:);
    x2 = x_n(segm(ii,2),:);
    s1 = x_n(i,:) - x1;
    s  = x2 - x1;
    da = abs(s1(1)*s(2) - s1(2)*s(1))/norm(s);
    
    v1 = x_n(segm(ii,2),:) - x_n(segm(ii,1),:);
    
    for ee = 1:size(ebnd_i,1)
      if sum(any(ebnd_i(ee,1)==id_list)) && sum(any(ebnd_i(ee,2)==id_list))
        
        %distance from the node i to the segment (infinite line)
        x1 = x_n(ebnd_i(ee,1),:);
        x2 = x_n(ebnd_i(ee,2),:);
        s1 = x_n(i,:) - x1;
        s  = x2 - x1;
        db = abs(s1(1)*s(2) - s1(2)*s(1))/norm(s);
        
        v2 = x_n(ebnd_i(ee,2),:) - x_n(ebnd_i(ee,1),:);
        if acos(abs(v1*v2')/(norm(v1)*norm(v2))) < 0.05 && abs(da-db)/norm(s)<0.05
          to_remove(ii) = 1;
          break;
        end
      end
    end
    %to_remove(ii) = 1;
  end
end
segm(to_remove==1,:)=[];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segm = segmIN_caseB(i,x_n,ebnd_i,ibnd,id_list,segm)
to_remove = zeros(size(segm,1),1);

for ii=1:size(segm,1)
  
  if any(segm(ii,1)==ibnd) && any(segm(ii,2)==ibnd)
    
    %distance from the node i to the segment (infinite line)
    x1 = x_n(segm(ii,1),:);
    x2 = x_n(segm(ii,2),:);
    s1 = x_n(i,:) - x1;
    s  = x2 - x1;
    da = abs(s1(1)*s(2) - s1(2)*s(1))/norm(s);
    
    v1 = x_n(segm(ii,2),:) - x_n(segm(ii,1),:);
    
    for ee = 1:size(ebnd_i,1)
      if sum(any(ebnd_i(ee,1)==id_list)) || sum(any(ebnd_i(ee,2)==id_list))
        
        %distance from the node i to the segment (infinite line)
        x1 = x_n(ebnd_i(ee,1),:);
        x2 = x_n(ebnd_i(ee,2),:);
        s1 = x_n(i,:) - x1;
        s  = x2 - x1;
        db = abs(s1(1)*s(2) - s1(2)*s(1))/norm(s);
        
        v2 = x_n(ebnd_i(ee,2),:) - x_n(ebnd_i(ee,1),:);
        if acos(abs(v1*v2')/(norm(v1)*norm(v2))) < 0.05 && abs(da-db)/norm(s)<0.05
          to_remove(ii) = 1;
          break;
        end
      end
    end
    %to_remove(ii) = 1;
  end
end
segm(to_remove==1,:)=[];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segm = segmIN_caseC(ibnd,segm)
to_remove = zeros(size(segm,1),1);
for ii=1:size(segm,1)
  if any(segm(ii,1)==ibnd) && any(segm(ii,2)==ibnd)
    to_remove(ii) = 1;
  end
end
segm(to_remove==1,:)=[];
