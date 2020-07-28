function a_near = QueriesAdjacency(x_n,x_a,range_a)
% function a_near = QueriesAdjacency(x_n,x_a,range_a)
% This funtion computes to each query point for x_a its neihbors in x_n

[n_a dim]= size(x_a);
if dim>1
  atria    = nn_prepare(x_n);
end

a_near  = {[]};
for i=1:n_a
    a_near{i} = [];
end

for i=1:n_a
  if dim>1
    [count, a_near0] = range_search(x_n, atria, x_a(i,:), range_a(i));
  else
    [count, a_near0] = range_search1D(x_n, x_a(i), range_a(i));
  end
  
  if count > 0
    a_nearI  = sort(a_near0{1,1});
    a_near{i}= a_nearI(a_nearI>0);
  end
end

function [count near] = range_search1D(x_s,x_a,range_a)
% This funtion computes to each x_b point its neihbors in x_a

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