function is_near = SamplesAdjacencyFEM(tri,gPts)
% function is_near = SamplesAdjacencyFEM(tri,gPts)
% This funtion computes to each sample point its neihbors in x_n given an input
% triangulation tri
%
% gPts is the number of Gauss points in each element



nElem= size(tri,1);
sPts = nElem*gPts;


is_near  = {zeros(sPts,1)};
for k=1:sPts
  is_near{k}=[];
end

for e=1:nElem
  eD = (e-1)*gPts;
  for k=1:gPts
    is_near{eD+k} = tri(e,:);
  end
end

