function [ xtab, weight ] = legendre_set_sqrtx2_01 ( norder )

%% LEGENDRE_SET_SQRTX2_01 sets a Gauss-Legendre rule for F(X) / SQRT(X) on [0,1].
%
%  Discussion:
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function is w(x) = 1 / sqrt ( x ).
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) F(X) / SQRT ( X ) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%  Modified:
%
%    13 October 2005
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Daniel Zwillinger, editor,
%    CRC Standard Mathematical Tables and Formulae,
%    CRC Press, 30th Edition, 2000, page 696.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%
  norder2 = 2 * norder + 1;

  [ xtab2, weight2 ] = legendre_set ( norder2 );

  xtab(1:norder) = xtab2(norder+2:2*norder+1).^2;

  weight(1:norder) = 2.0 * weight2(norder+2:2*norder+1);
