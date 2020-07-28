function [x_n i_bd i_in]= HexagonalGrid(Nx, Ny, Lx, Ly, centre)
% [x_n i_bd i_in]= HexagonalGrid(Nx, Ny, Lx, Ly, centre)
% function to build points in an 'hexagonal grid' arrangment, but in a very
% easy fashion.
% 
% Nx    : number of nodes in x's direction
% Ny    : number of nodes in y's direction
% Lx    : the length of the mesh in x.
% Ly    : the length of the mesh in y.
% centre: coordinates of the mesh geometric centre

%the origin is the first point (to avoid duplication is outside of the loop)

nPts = Nx*Ny;
bPts = 2*Nx + 2*(Ny-2);
    
if mod(Nx-1,2) == 0  %if the number of nodes in the y-axis is odd
  nPts = nPts + (Nx-1)/2;
else
  nPts = nPts + Nx/2;
  bPts = bPts + 1;
end

iPts = nPts - bPts;
	
x_n  = zeros(nPts,2);   % node coordinates
i_bd = zeros(bPts);     % boundary node indexes
i_in = zeros(iPts);     % interior node indexes


dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

lx = 0.5*Lx - centre(1);
ly = 0.5*Ly - centre(2);

ptId=1;

for i=1:Ny
  y = (i-1)*dy - ly;
  if mod(i-1,2) == 0
    for j=1:Nx
      x_n(ptId,1) = (j-1)*dx - lx;
      x_n(ptId,2) = y;
      ptId = ptId+1;
    end
  else
    x_n(ptId,1) = -lx;
    x_n(ptId,2) =  y;
    ptId = ptId+1;
    for j=2:Nx
      x_n(ptId,1) = (j-1.5)*dx - lx;
      x_n(ptId,2) = y;
      ptId = ptId+1;
    end
    x_n(ptId,1) = Lx - lx;
    x_n(ptId,2) = y;
    ptId = ptId+1;
  end
end
x_n = unique(x_n,'rows');