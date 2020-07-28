function [ xtab, weight ] = cheb_u_set ( norder )

%% CHEB_U_SET sets abscissas and weights for Gauss-Chebyshev quadrature.
%
%  Discussion:
%
%    The integration interval is [ -1, 1 ].
%
%    The weight function is w(x) = sqrt ( 1 - x**2 ).
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) SQRT ( 1 - X**2 ) * F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    If NORDER points are used, then Gauss-Chebyshev quadrature
%    will compute the integral exactly, whenever F(X) is a polynomial
%    of degree 2*NORDER-1 or less.
%
%    The abscissas are zeroes of the Chebyshev polynomials
%    of the second kind, U(NORDER)(X).
%
%  Modified:
%
%    12 October 2005
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Daniel Zwillinger, editor,
%    Standard Mathematical Tables and Formulae,
%    30th Edition,
%    CRC Press, 1996.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule,
%    which are all equal to PI / NORDER.
%
  for i = 1 : norder
    angle = i * pi / ( norder + 1 );
    xtab(i) = cos ( angle );
    weight(i) = pi * ( sin ( angle ) )^2 / ( norder + 1 );
  end
