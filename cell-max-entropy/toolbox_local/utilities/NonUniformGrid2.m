function x_a = NonUniformGrid2(D,N,L)
% function x_a = NonUniformGrid2(D,N,L)
% function to build points in a 'non uniform grid' arrangment.
%
% 
% D: dimension, could be 2 or 3. Default value is 2.
% N: number of samples in r and theta N=[Nx Ny Nz] (Default: [10 10 10])
% L: is the length of the square/cube side L=[Lx Ly Lz] (Default: [1 1 1])

if nargin<3
  L = [1 1 1];
elseif length(L)<2
  L = L*[1 1 1];
end

if nargin<2
  N = [10 10 10];
elseif length(N)<2
  N = N*[1 1 1];
end
if nargin<1
  D = 2;
end
if nargin>4
  error('Only you can set for input arguments: D, N, L, type')
end


if D==2
  x_a = nug2D(N,L);
elseif D==3
  x_a = [];
  bz  = 0.5*pi/(N(3)-1);
  xy  = nug2D(N,L);
  Nxy = N(1)*N(2);
  for i=0:N(3)-1
    z = L(3)*sin(bz*i)*sin(bz*i)*ones(N(1)*N(2), 1);
    x_a(i*Nxy+1:(i+1)*Nxy,:) = [xy z];
  end
else
  error('Only you can set D=2,3')
end


%Loop to create the points in a rectangular no uniform grid in 2D
function x_a = nug2D(N,L)
k=1;
x_a = zeros(N(1)*N(2),2);
bx = 0.5*pi/(N(1)-1);
by = 0.5*pi/(N(2)-1);
for j=0:N(2)-1
  y = L(2)*sin(by*j)*sin(by*j);
  for i=0:N(1)-1
    x_a(k,:) = [L(1)*sin(bx*i)*sin(bx*i) y];
    k=k+1;
  end
end
