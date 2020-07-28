function [x_nodes,x_nodes0,id_bd] = SquareNodes(ghost,N, L)
% function [x_nodes x_nodes0 id_bd] = SquarePlateNodes(Nx, L, meshtype)
% This function returns the set of nodes for a square plate of length L and its
% extended set of nodes + ghost
%
% INPUT:
% <> N: number of nodes in the x axis (horizontal)
% <> L: length of the side
%    UniformGrid(N, N, L, L, centre)
%
% OUTPUT:
% <> x_nodes: full arrangement of
% <> x_nodes0: interior nodes
% <> id_bd: index identifiers for the boundary (physical).

if ghost==0
  gh = 0;
else
  gh = 1;
end

% centre of the plate
centre  = [0 0];

h = L/(N-1);

x_nodes0 = UniformGrid2D(N, N, L, L, centre);

if gh == 1
  x_nodes = UniformGrid2D(N+2, N+2, L+2*h, L+2*h, centre);
  nPts    = length(x_nodes);
  ids     = (1:nPts);
  id_bd.G{1} = ids(abs(x_nodes(:,2)+(L/2+h))<1.e-6); % (y=-L/2-h)
  id_bd.G{2} = ids(abs(x_nodes(:,1)-(L/2+h))<1.e-6); % (x=L/2+h)
  id_bd.G{3} = ids(abs(x_nodes(:,2)-(L/2+h))<1.e-6); % (y=L/2+h)
  id_bd.G{4} = ids(abs(x_nodes(:,1)+(L/2+h))<1.e-6); % (x=-L/2-h)
  id_bd.G = unique([id_bd.G{:}]);
else
  x_nodes = x_nodes0;
  nPts    = length(x_nodes);
  ids     = (1:nPts);
  id_bd.G = [];
end

%Select the nodes that have Dirichlet BC
id_bd.D = ids(x_nodes(:,1) >= 0 & abs(x_nodes(:,2)) < 1e-06);