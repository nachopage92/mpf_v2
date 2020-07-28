function [x_a r]= NonUniformGrid(D,N,L,type)
% function [x_a r]= NonUniformGrid(D,N,L,type)
% function to build points in a 'non uniform grid' arrangment. The points are put in 5
% diferent configurations. All of them are build, let say, in a radial way... check the
% code or plot same examples for more details.
%
% D: dimension, D=2 is the default value
% N: number of samples in r and theta, N=10
% L: is the length of the square side, L=1
% type: set the mesh type, type=1,..,5.

if nargin<4
  type = 1;
end
if nargin<3
  L = 1;
end
if nargin<2
  N = 10;
end
if nargin<1
  D = 2;
end
if nargin>4
  error('Only you can set for input arguments: D, N, L, type')
end

r = zeros(N-1,1);
switch (type)
  case 1
    a = L*sqrt(2)/(N-1);
    for i=1:N-1
      r(i) = a*i;
    end
  case 2
    a = L*sqrt(2)/((N-1)*(N-1));
    for i=1:N-1
      r(i) = a*i*i;
    end
  case 3
    a = L*sqrt(2)/((N-1)*(N-1)*(N-1));
    for i=1:N-1
      r(i) = a*i*i*i;
    end
  case 4
    a = L*sqrt(2);
    b = 0.5*pi/(N-1);
    for i=1:N-1
      r(i) = a*sin(b*i)*sin(b*i);
    end
  case 5
    a = L*sqrt(2)/(N-1);
    b = 0.5*pi/(N-1);
    for i=1:N-1
      r(i) = a*i*sin(b*i)*sin(b*i);
    end
end


if D==2
  %the origin is the first point (to avoid duplication this is put outside of the loop)
  x_a(1,:)=[0 0];
  x = nug2D(N,L,r);
  x_a=[x_a; x];
elseif D==3
  x_a(1,:) = [0 0 0];
  xy = nug2D(N,L,r);
  z  = zeros(N*N-1, 1);
  x_a = [x_a; xy z];
  for i=1:N-1
    x_a(i*N*N+1,:) = [0 0 r(i)];
    z = r(i)*ones(N*N-1, 1);
    x_a = [x_a; xy z];
  end
else
  error('Only you can set D=2,3')
end


%Loop to create the points in a rectangular no uniform grid in 2D
function x_a = nug2D(N,L,r)

j=1;
Th = pi/4;
x_a=zeros((N-1)*(N-1),2);
for ir=1:N-1
  rri=r(ir)*cos(Th);
  for it=1:ir+1
    thi= (it-1)*Th/(ir);
    rr = rri/cos(thi);
    x  = rr*cos(thi);
    y  = rr*sin(thi);
    x_a(j,:)=[x y];
    j=j+1;
  end

  rri=r(ir)*sin(Th);
  for it=1:ir
    thi= it*Th/(ir)+Th;
    rr= rri/sin(thi);
    x = rr*cos(thi);
    y = rr*sin(thi);
    x_a(j,:)=[x y];
    j=j+1;
  end
end

nPts = length(x_a);
ids  = (1:nPts);
id_x0= ids(abs(x_a(:,1))<1e-6);
id_y0= ids(abs(x_a(:,2))<1e-6);
id_xL= ids(abs(x_a(:,1)-L)<1e-6);
id_yL= ids(abs(x_a(:,2)-L)<1e-6);
x_a(id_x0,1) = 0;
x_a(id_xL,1) = L;
x_a(id_y0,2) = 0;
x_a(id_yL,2) = L;