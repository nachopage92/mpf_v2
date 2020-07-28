function [h_n,iters,h_n0]= NodalSpacing(x_n,options,id_bd)
% function [h_n,iters,h_n0]= NodalSpacing(x_n,options,id_bd)
% Computation of the beta parameter and the range for searching nearest
% neighbours
% INPUT:
%   x_n     : pointset, (NxDim)
%   options : this structure must be contain the main settings for computing
%             the nearest neighbors, such that gamma, knn...
%   id_bd   : indices for the nodes on the boundary (or close to it), these will have
%             a special treatment

if nargin < 2
  options.exist = 1;
end
if nargin<3
  id_bd = [];
end

nPts = size(x_n,1);

if isfield(options,'knn')
  knn = options.knn;
else
  knn = 0;
end

if isfield(options,'dim')
  dim  = options.dim;
else
  error('dim is not defined');
end


if isfield(options,'gamma')
  gammaB = options.gamma;
else
  gammaB = 2;
end


if knn>0
  h_n0 = initialize_h(x_n,options,id_bd);
  if dim==1
    [dist,nn_id] = nn_search1D(x_n, x_n, knn);
    dist =[zeros(nPts,1) dist];
    nn_id=[(1:nPts)' nn_id];
  else
    atria=nn_prepare(x_n);
    [nn_id,dist] = nn_search(x_n, atria, x_n, knn+1);
  end

  TolB  = 1.e-3;
  iters = zeros(nPts,1);
  h_n   = zeros(nPts,1);
  for i=1:nPts
    hi0  = h_n0(i);
    %di   = dist{i};
    %nn_i = nn_id{i};
    di   = dist(i,:);
    nn_i = nn_id(i,:);
    h_ni = h_n0(nn_i);
    dbeta= 1;
    it   = 0;
    him  = hi0;
    
    while dbeta>TolB && it<100
      betai= gammaB/him^2;
      %expb = exp(-betai*(di-hi0).^2);
      expb = exp(-betai*(di).^2);
      wi   = expb/sum(expb);
      hi   = sum(h_ni.*wi);
      % hi   = sum(di.*wi);
      dbeta= betai-gammaB/hi^2;
      dbeta= abs(dbeta/betai);
      it   = it+1;
      
      if it>95
          fprintf(1,'i=%2d  it=%3d  h=[%5.4f %5.4f %5.4f]\n', i, it, hi0, him, hi);
      end
      
      if it>20
        %w  = exp(-betai*(hi/hi0-1).^2);
        w  = 0.5;
      else
        w  = 1.0;
      end
      him= hi*w + him*(1-w);
    end
    h_n(i)    = him;
    iters(i)  = it;
  end
elseif knn == 0
  if isfield(options,'spacing')
    h_n = options.spacing * ones(nPts,1);
  else
    error('SPACING is needed (knn=0)')
  end
else
  error('ERROR')
end

%%
function [h_n,dist,nn_id] = initialize_h(x_n,options,id_bd)
% interior node indices
nPts = length(x_n);
bPts = length(id_bd);
id_in = (1:nPts);
id_in(id_bd)=[];
iPts = length(id_in);
dim  = options.dim;
h_n  = zeros(1,nPts);
dist = {zeros(nPts,1)};
nn_id= {zeros(nPts,1)};


if dim==1
  %boundary nodes
  [dj_bd,nn_bd] = nn_search1D(x_n, x_n(id_bd), 1);
  for j=1:bPts
    i     = id_bd(j);
    di    = dj_bd(j,:);
    h_n(i)= di(1);
    %     h_n(i)= 0.5*(di(1)+0.5*di(2));
    dist{i}= di;
    nn_id{i}=nn_bd(j,:);
  end
  
  %interior nodes
  [dj_in,nn_in] = nn_search1D(x_n, x_n(id_in), 2);
  for j=1:iPts
    i     = id_in(j);
    di    = dj_in(j,:);
    h_n(i)= 1/2*(di(1)+di(2));
    %     h_n(i)= 1/3*(di(1)+di(2)+di(3)/2);
    %     h_n(i)= 1/4*(di(1)+di(2)+di(3)/2+di(4)/2);
    %     h_n(i)= 1/3*(0.5*(di(1)+di(2)) + 0.25*(di(3)+di(4)) + di(5)/3);
    dist{i}= di;
    nn_id{i}=nn_in(j,:);
  end
  
elseif dim==2
  %Boundary nodes are split in C:corners and E:edges
  kC = 1; kE=1;
  if ~isempty(id_bd)
    x_bd = x_n(id_bd,:);
    atriaBD= nn_prepare(x_bd);
    id = (1:bPts);
    nn_bd  = nn_search(x_bd, atriaBD, id, 2, 0);
    id_bdC = zeros(1,bPts);
    id_bdE = zeros(1,bPts);
    for j=1:bPts
      i  = id_bd(j);
      jd = nn_bd(j,:);
      v1 = x_bd(jd(1),:)-x_bd(j,:);
      v2 = x_bd(jd(2),:)-x_bd(j,:);
      sp = v1*v2';
      if sp<0
        id_bdE(kE) = i;
        kE = kE+1;
      else
        id_bdC(kC) = i;
        kC = kC+1;
      end
    end
  end
  bPtsC = kC-1;
  bPtsE = kE-1;
  
  %% spacing stimation
  atria = nn_prepare(x_n);
  %boundary nodes: CORNERS
  if bPtsE>0
    id_bdC= id_bdC(1:bPtsC);
    [nn_bd,dj_bdC] = nn_search(x_n, atria, id_bdC, 5, 0);
    for j=1:bPtsC
      i      = id_bdC(j);
      dj     = dj_bdC(j,:);
      h_n(i) = 1/5*( 2*mean(dj(1:2)) + 0.5*(1+1/sqrt(2))*dj(3) + sum(dj(4:5))*0.5 );
      dist{i}= dj;
      nn_id{i}=nn_bd(j,:);
    end
  end
  
  %boundary nodes: EDGES
  if bPtsE>0
    id_bdE= id_bdE(1:bPtsE);
    [nn_bd,dj_bdE] = nn_search(x_n, atria, id_bdE, 5, 0);
    for j=1:bPtsE
      i       = id_bdE(j);
      dj      = dj_bdE(j,:);
      h_n(i)  = 1/5*( sum(dj(1:3)) + 0.5*(1+1/sqrt(2))*dj(4) + 0.5*(1/sqrt(2)+0.5)*dj(5));
      dist{i} = dj;
      nn_id{i}=nn_bd(j,:);
    end
  end
  
  %[nn_in dj_in] = nn_search(x_n, atria, id_in, 5, 0);
  [nn_in,dj_in] = nn_search(x_n, atria, id_in, 4, 0);
  %interior nodes
  for j=1:iPts
    i       = id_in(j);
    dj      = dj_in(j,:);
    h_n(i)  = mean(dj);  %mean(dj(1:4));
    %     h_n(i) = 1/10*( sum(dj(1:4)) + (dj(5)+dj(6))*0.5*(1+1/sqrt(2)) + ...
    %       sum(dj(7:8))*0.5*(1/2+1/sqrt(2)) + sum(dj(9:10))*0.5);
    dist{i} = dj;
    nn_id{i}= nn_in(j,:);
  end
else
  knn_ = min([2*dim options.knn]);
  atria=nn_prepare(x_n);
  [near_tmp,dist_tmp]= nn_search(x_n,atria,id_in,knn_,0);
  h_n   = mean(dist_tmp,2)';
  %h_min= prctile(h_n,10);
  %h_max= prctile(h_n,90);
  %h_n(h_n<h_min) = h_min;
  %h_n(h_n>h_max) = h_max;
  for i=1:nPts
    dist{i} = dist_tmp(i,:);
    nn_id{i}= near_tmp(i,:);
  end
  
  %hcut = quantile(h_n,0.05);
  %h_n(h_n<hcut)=hcut;
%   near_tmp(1,:)-1
%   dist_tmp(1,:)
end



%% This funtion computes for the point x_a its k-nearest neihbors in x_s
function [dist,near] = nn_search1D(x_s,x_a,knn)
aPts = length(x_a);
dist = zeros(aPts,knn);
near = zeros(aPts,knn);
for i=1:aPts
  di_a = (x_s-x_a(i)).^2;
  [di_a,nn_a] = sort(di_a);
  knn_  = min([length(di_a) knn+1]);
  dist(i,:) = sqrt(di_a(2:knn_)');
  near(i,:) = nn_a(2:knn_)';
end

%% This funtion computes for the x_a point its nearest neihbors in x_s such
%% that they are inside of its range
function [count,near] = range_search1D(x_s,x_a,range_a)
% Naive construction of the neighbor list
sPts = length(x_s);
all_ = (1:sPts);

%Naive find in range
dist = (x_s-x_a).^2;
dist = sqrt(dist);
near_= all_(dist<range_a);
count= length(near_);
dist_= dist(near_);
near = {near_,dist_};