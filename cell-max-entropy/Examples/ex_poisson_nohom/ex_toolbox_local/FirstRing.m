function [n_nears,e_nears]= FirstRing(X,T,selfexcl)
% function [n_nears,e_nears] = FirstRing(X,T,selfexcl)
%
% this function receives a set of points X of size [nPts x nDim] and a
% connectivity T of size [nElem x nVert]. For each node returns a list of
% its connected neighbor nodes (also itself is contained but can be
% removed).
% Aditionally a list of the elements which share each node can be given.
%
% INPUT:
%  X        : point coordinates
%  T        : connectivity
%  selfexcl : Self exclusion from the lisnt. Using selfexcl = 1 means: exclude self-matches
%             Default value is elfexcl = 0.
%
% OUTPUT:
%   n_nears: nodes connected by a edge
%   e_nears: elements wich share a given node/vertex 

if nargin<2
  error('bad setting')
end
if nargin==2
  selfexcl=0;
end

% nElem: number of elements    nVert: number of vertices
[nElem,nVert]=size(T);

% nPts: number of nodes, which belong to nDim
nPts=size(X,1);

%% ELEMENTS
e_counter=zeros(nPts,1);
for e=1:nElem
  j = T(e,:);
  e_counter(j) = e_counter(j)+1;
end

e_nears={[]};
for i=1:nPts
  e_nears{i}=zeros(e_counter(i),1);
end

e_counter=zeros(nPts,1);
for e=1:nElem
  for i=1:nVert
    j = T(e,i);
    jj= e_counter(j)+1;
    e_nears{j}(jj) = e;
    e_counter(j)   = jj;
  end
end


%% NODES
i_nears=zeros(nVert*max(e_counter),1);
n_nears={zeros(nPts,1)};
for i=1:nPts
  jj=1;
  for ie=1:length(e_nears{i})
    elem=T(e_nears{i}(ie),:);
    for ii=1:nVert
      i_nears(jj)=elem(ii);
      jj=jj+1;
    end
  end
  i_nears_    = i_nears(1:(jj-1));
  n_nears{i}  = unique(sort(i_nears_));
  if selfexcl==1
    n_nears{i}= setdiff(n_nears{i},i);
  end
end