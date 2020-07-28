function [x_s, w_g, tri] = MakeGLSamples2D(x_a, options)
% function [x_s, w_g, tri] = MakeGLSamples2D(x_a, options)
%
% Given a input point data of size Nx2 this function performs the Delaunay
% triangulation and creates the Gauss-Legendre points for numerical
% integration
%
% INPUT
% x_a: input 2D points
% options.quality: quality mesure for the triangles (default 0.5)
% options.orderGL: order of the cubature to be used 
%
% OUTPUT
% x_s: sample gauss points in the physical space
% w_g: weights for the numerical integration at each sample point
% tri: connectivity of the triangulation of Delaunay 

%GAUSS POINTS in the reference element
order    = options.orderGL;
[xi,wei] = GaussLegendreCubature(order,2);
gPts = length(wei);

%DELAUNAY TRIANGULATION
if isfield(options,'tri')
  tri = options.tri;
else
  tri = delaunay(x_a(:,1),x_a(:,2));%,{'Qt','Qbb','Qc','Qz'});
end
nTri= length(tri);
x_s = zeros(nTri*gPts, 2);
w_g = zeros(nTri*gPts, 1);


%quality mesure for the triangles
if isfield(options,'quality')
  quality = options.quality;
else
  quality = 0.5;
end

% triangular elements are made
ig = 0;
ie = 0;
id_e = [];
for i=1:nTri
  Xa    = x_a(tri(i,:),:); %vertices
  [q,A] = TriangleQuality(Xa);
  if q > quality
    wg = A*wei;
    for j=1:gPts %for each Gauss point 
      ig = ig+1;
      %triangular coordinates
      w1 = 1-xi(j,1)-xi(j,2);
      w2 = xi(j,1);
      w3 = xi(j,2);
      x_s(ig,:)= w1*Xa(1,:) + w2*Xa(2,:) + w3*Xa(3,:);  %Gauss point in the physical coordinates
      w_g(ig)  = wg(j); 
    end
  else
    ie = ie + 1;
    id_e(ie) = i;
  end
end
id_g = (ig+1:nTri*gPts);
x_s(id_g,:) = [];
w_g(id_g)   = [];
tri(id_e,:) = [];


% returns a triangle quality measure q and the triangle area a
function [q,a] = TriangleQuality(p)

% Triangle sides
r23x = p(3,1) - p(2,1);
r23y = p(3,2) - p(2,2);
r31x = p(1,1) - p(3,1);
r31y = p(1,2) - p(3,2);
r21x = p(1,1) - p(2,1);
r21y = p(1,2) - p(2,2);

% Area
a = abs(r31x*r23y-r31y*r23x) * 0.5;

h1= sqrt(r21x^2+r21y^2);
h2= sqrt(r23x^2+r23y^2);
h3= sqrt(r31x^2+r31y^2);

q = 4*a*sqrt(3)/(h1^2 + h2^2 + h3^2);