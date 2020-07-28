function [u_sol, du_sol] = AnalyticalSolution(bc_homog, x_n)
% function [u_sol, du_sol] = AnalyticalSolution(bc_homog, x_n)
%
% CODINA'S Heat Equation problem
% bc_homog   0: homogeneous; 1: non-homogeneous
  
% val;			Function value for specific x and y
% dval_x;		X Derivative of the function for specific x and y
% dval_y;		y Derivative of the function for specific x and y
	
N   = 200;
epi = 8.0/(pi*pi*pi);
iks = zeros(N,1);
for k=1:2:N %even values only
  iks(k) = 1/ ( k*k*k * sinh( k*pi ) );
end

nPts  = size(x_n,1);
u_sol = zeros(nPts,1);
du_sol= zeros(nPts,2);

for i=1:nPts
  x = x_n(i,1);
  y = x_n(i,2);

  val    = 0.0;
  dval_x = 0.0;
  dval_y = 0.0;

  for k=1:2:N %even values only
    ik     = iks(k);
    val    = val    + sin( k*pi*x ) * sinh( k*pi*y ) * ik;
    dval_x = dval_x + k*pi * cos( k*pi*x ) * sinh( k*pi*y ) * ik;
    dval_y = dval_y + k*pi * sin( k*pi*x ) * cosh( k*pi*y ) * ik;
  end

  % 0: homogeneous; 1: non-homogeneous
  u_sol(i)    = epi*val;
  du_sol(i,:) = epi*[dval_x, dval_y];

  if bc_homog == 0
    u_sol(i)    = u_sol(i) - x*y*(1.0-x);
    du_sol(i,:) = du_sol(i,:) - [y*(1.0-2.0*x), x*(1.0-x)];
  end
end
