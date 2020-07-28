function [u_sol, du_sol, lu_sol] = Unstr_AnalyticalSolution(x_n)
% function [u_sol, du_sol] = AnalyticalSolution(x_n)
%
% u_sol;			Function value for specific x and y
% du_sol_x;		X Derivative of the function for specific x and y
% du_sol_y;		y Derivative of the function for specific x and y
% lu_sol      Laplacian of u

A=[10 50 100 50]';
B=[180 450 800 1000]';
C=[0.51 0.52;
  0.31 0.34;
  0.73 0.71;
  0.28 0.72]; %locaction of sources

nPts   = size(x_n,1);

u_sol  = zeros(nPts,1);
du_sol = zeros(nPts,2);
lu_sol = zeros(nPts,1);

for i=1:nPts
  dx  = ones(4,1)*x_n(i,:) - C;
  
  val = A.*exp(-B.*sum(dx.^2,2));
  
  u_sol(i)    = sum(val);
  du_sol(i,1) = sum(val.*(-2*B.*dx(:,1)));
  du_sol(i,2) = sum(val.*(-2*B.*dx(:,2)));
  
  lu_sol(i)   = 4*sum(B.*val.*(B.*sum(dx.^2,2)-ones(4,1)));
end
