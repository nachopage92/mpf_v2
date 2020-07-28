function [x_n,grid]= UniformGrid2D(Nx, Ny, Lx, Ly, centre)
% [x_n,grid]= UniformGrid2D(Nx, Ny, Lx, Ly, centre)
% function to build points in an 'uniform grid' arrangment, but in a very easy fashion.
% 
% INPUT:
%   Nx    : number of nodes in x's direction
%   Ny    : number of nodes in y's direction
%   Lx    : the length of the mesh in x.
%   Ly    : the length of the mesh in y.
%   centre: coordinates of the mesh geometric centre
%
% OUTPUT:
%   x_n       : nodes coordinates
%   grid.quad : quadrilateral connectivity (mesh)
%   grid.tri  : triangular connectivity (mesh)

%the origin is the first point (to avoid duplication is outside of the loop)

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

Lx = 0.5*Lx - centre(1);
Ly = 0.5*Ly - centre(2);

x_n=zeros(Nx*Ny,2);

ptId = 1;
for i=0:Ny-1
  y= i*dy - Ly;
  for j=0:Nx-1
    x_n(ptId,1) = j*dx - Lx;
    x_n(ptId,2) = y;
    ptId = ptId+1;
  end  
end

%Quadrilateral connectivity
nElem=(Nx-1)*(Ny-1);
quad = zeros(nElem,4);
e=1;
for i=1:Ny-1
  for j=1:Nx-1
    n1 = j+(i-1)*Nx;
    n2 = j+i*Nx;
    quad(e,:) = [n1 n1+1 n2+1 n2];
    e = e + 1;
  end  
end
grid.quad = quad;


%Triangular mesh connectivity: interlaced
nElem= (Nx-1)*(Ny-1)*2;
tri  = zeros(nElem,3);
e=1;
for i=1:Ny-1
  for j=1:Nx-1
    n1 = j+(i-1)*Nx;
    n2 = n1+1;
    n3 = j+i*Nx+1;
    n4 = n3-1;
    
    if mod(i+j,2)==0
      tri(e,:) = [n1 n2 n3];
      e = e + 1;
      tri(e,:) = [n1 n3 n4];
      e = e + 1;
    else
      tri(e,:) = [n1 n2 n4];
      e = e + 1;
      tri(e,:) = [n2 n3 n4];
      e = e + 1;
    end
  end
end
grid.tri = tri;
