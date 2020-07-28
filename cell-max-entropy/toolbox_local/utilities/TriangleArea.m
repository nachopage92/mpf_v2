function ar = TriangleArea(p)
% this function computes the triangle area
% p: Corner points of the triangle

% Triangle sides
r23x = p(3,1) - p(2,1);
r23y = p(3,2) - p(2,2);
r31x = p(1,1) - p(3,1);
r31y = p(1,2) - p(3,2);

% Area
ar = abs(r31x*r23y-r31y*r23x) * 0.5;