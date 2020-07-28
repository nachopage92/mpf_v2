function [ x, w ] = fejer1_compute ( n )

%% FEJER1_COMPUTE computes a Fejer type 1 quadrature rule.
%
%  Discussion:
%
%    This method uses a direct approach.  The paper by Waldvogel
%    exhibits a more efficient approach using Fourier transforms.
%
%  Modified:
%
%    03 March 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Philip Davis, Philip Rabinowitz,
%    Methods of Numerical Integration,
%    Second Edition,
%    Dover, 2007,
%    ISBN: 0486453391,
%    LC: QA299.3.D28.
%
%    Walter Gautschi,
%    Numerical Quadrature in the Presence of a Singularity,
%    SIAM Journal on Numerical Analysis,
%    Volume 4, Number 3, 1967, pages 357-362.
%
%    Joerg Waldvogel,
%    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
%    BIT Numerical Mathematics,
%    Volume 43, Number 1, 2003, pages 1-18.
%
%  Parameters:
%
%    Input, integer N, the order of the rule.
%
%    Output, real X(N), W(N), the abscissas and weights of the rule.
%
  theta(1:n) = ( 2*n-1 : -2 : 1 ) * pi / ( 2 * n );
  x(1:n) = cos ( theta(1:n) );

  for i = 1 : n
    w(i) = 1;
    for j = 1 : floor ( n / 2 )
      w(i) = w(i) - 2 * cos ( 2 * j * theta(i) ) / ( 4 * j * j - 1 );
    end
  end

  w(1:n) = 2 * w(1:n) / n;
