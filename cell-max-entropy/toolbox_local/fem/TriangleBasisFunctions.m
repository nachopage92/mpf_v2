function [Na,dNa]= TriangleBasisFunctions(Xa,x)
% this function computes the triangular basis functions and their gradients at x (P1 element)
% Xa: Corner points of the triangle

A = TriangleArea(Xa);

A = 0.5/A;

xx=x(1);
yy=x(2);

x1=Xa(1,1);    y1=Xa(1,2);
x2=Xa(2,1);    y2=Xa(2,2);
x3=Xa(3,1);    y3=Xa(3,2);

N1 = A*((x2*y3-x3*y2) + (y2-y3)*xx + (x3-x2)*yy);

N2 = A*((x3*y1-x1*y3) + (y3-y1)*xx + (x1-x3)*yy);

N3 = A*((x1*y2-x2*y1) + (y1-y2)*xx + (x2-x1)*yy);

Na = [N1;
      N2; 
      N3];

    
if nargout==1
  return;
end

%Gradient of Na
dN1 = A * [(y2-y3), (x3-x2)];

dN2 = A * [(y3-y1), (x1-x3)];

dN3 = A * [(y1-y2), (x2-x1)];

dNa = [dN1; 
       dN2;
       dN3];