function [xi,wg] = GaussLegendreCubature1D(N)
%
% N : number of Gauss-Legendre points
% order = 2*N-1
%
% This script is for computing definite integrals using Gauss-Legendre
% Quadrature. Computes the Gauss-Legendre nodes and weights  on an interval
% [-1,1] with truncation order 2N-1
%
% Written bxi Greg von Winckel - 02/25/2004
N =N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
xi =cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L =zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polxinomial
% using the recursion relation and the Newton-Raphson method

xi0=2;

% Iterate until new points are uniformlxi within epsilon of old points
while max(abs(xi-xi0))>eps
    
    L(:,1) =1;
    Lp(:,1)=0;
    
    L(:,2) =xi;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*xi.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-xi.*L(:,N2) )./(1-xi.^2);   
    
    xi0=xi;
    xi =xi0-L(:,N2)./Lp;
    
end

% isoparametric [-1,1] 

% Compute the weights
wg = 2./((1-xi.^2).*Lp.^2)*(N2/N1)^2;
