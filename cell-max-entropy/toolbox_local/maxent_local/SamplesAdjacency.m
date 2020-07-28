function [is_near,ds_near]= SamplesAdjacency(x_n,x_s,range_n)
% function [is_near,ds_near]= SamplesAdjacency(x_n,x_s,range_n)
% This funtion computes to each sample point its neihbors in x_n
%
% is_nears: index identifiers for the nearest nodes to each sample
% ds_nears: distance between the nearest nodes and each sample

nPts = size(x_n,1);
sPts = size(x_s,1);
nDim = size(x_n,2);

if nDim>nPts
  error('DIM is bigger than number of nodes, maybe you needs transpose your data')
end
if nargin<3
  error('Too few arguments, the use of this function is: SamplesAdjacency(x_n,x_s,range_n)')
end

if nDim==1
  x_n=[x_n, zeros(nPts,1)];
  x_s=[x_s, zeros(sPts,1)];
end


t=cputime;
atria   = nn_prepare(x_s);
fprintf(1,'\tComputing ATRIA tree in        %7.2fs\n',cputime-t);


%% =======================================================
t=cputime;

in_near = {[]};
dn_near = {[]};
scount  = zeros(sPts,1);

diag_box =norm(max(x_n)-min(x_n));

%if the variation isn't so big then efective support is considered fix
if std(range_n/diag_box)<1e-6 
  range_ = max(range_n);
  [ncount,nd_near] = range_search(x_s, atria, x_n, range_);
  for i=1:nPts
    if ncount(i) > 0
      n_near = nd_near{i,1};
      d_near = nd_near{i,2};
      
      in_near{i} = n_near;
      dn_near{i} = d_near;
      scount(n_near) = scount(n_near) + 1;
    end
  end
else
  ncount  = zeros(nPts,1);
  for i=1:nPts
    [ncount(i),nd_near] = range_search(x_s, atria, x_n(i,:), range_n(i));
    
    if ncount(i) > 0
      n_near = nd_near{1,1};
      d_near = nd_near{1,2};
      
      in_near{i} = n_near;
      dn_near{i} = d_near;
      scount(n_near) = scount(n_near) + 1;
    end
  end
end
clear atria nd_near
fprintf(1,'\tComputing AdjNears part A in   %7.2fs\n',cputime-t);

%% =======================================================
t=cputime;

is_near = {[]};
ds_near = {[]};

for k=1:sPts
  is_near{k} = zeros(1,scount(k));
  ds_near{k} = zeros(1,scount(k));
  scount(k)  = 1;
end
for i=1:nPts
  if ncount(i) > 0
    n_near = in_near{i};
    d_near = dn_near{i};
    
    for k=1:length(n_near)
      s   = n_near(k);
      d   = d_near(k);
      ipos = scount(s);
      is_near{s}(ipos) = i;
      ds_near{s}(ipos) = d;
      scount(s) = ipos + 1;
    end
  end
end
fprintf(1,'\tComputing AdjNears part B in   %7.2fs\n',cputime-t);

t=cputime;
for k=1:sPts
  % the adjacency list is rearanged such that the first indexes are related
  % with the closest nodes to each sample point.
  [dist, ids] = sort(ds_near{k},'ascend');
  is_near{k}  = is_near{k}(ids);
  ds_near{k}  = dist;
end
fprintf(1,'\tComputing AdjNears sorting in  %7.2fs\n',cputime-t);
