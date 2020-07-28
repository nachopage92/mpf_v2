function [u_exact,sig_exact]=Beam_AnalyticalSolution(pos,P,D,L,E,nu)
% function [u_exact,sig_exact]=Beam_AnalyticalSolution(pos,P,D,L,E,nu)
% 
% Analytical solution for the Timoshenko's cantilever beam is computed at x
% INPUT PARAMETERS:
% Boundary:
%   P: maximum load, the distribution is parabolic.
% Geometrical: 
%   D:diameter,
%   L:length
% Material: 
%   E:Young's modulus, 
%   nu:Poisson


u_exact  = zeros(size(pos,1),2);
sig_exact= zeros(size(pos,1),3);

%plane strain
E_hat = E/(1-nu^2);
nu_hat= nu/(1-nu);
I = D^3/12;

x = pos(:,1);
y = pos(:,2);

u_exact(:,1) = -P/6/E_hat/I*y.*( ...
                                 (6*L-3*x).*x + ...
                                 (2+nu_hat)*(y.^2 - D^2/4)  );

u_exact(:,2) = P/6/E_hat/I*( ...
                             3*nu_hat*(y.^2).*(L-x) + ...
                             (4+5*nu_hat)*D^2/4*x + ...
                             (3*L-x).*(x.^2)   );

sig_exact(:,1) = -P/I*(L-x).*y;
sig_exact(:,2) = 0;
sig_exact(:,3) =  P/2/I*( (D/2)^2 - y.^2 );
