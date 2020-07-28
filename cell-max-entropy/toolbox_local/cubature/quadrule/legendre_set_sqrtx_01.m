function [ xtab, weight ] = legendre_set_sqrtx_01 ( norder )

%% LEGENDRE_SET_SQRTX_01 sets a Gauss-Legendre rule for SQRT(X) * F(X) on [0,1].
%
%  Discussion:
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function is w(x) = sqrt(x).
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) SQRT ( X ) * F(X) dX =
%      Integral ( 0 <= Y <= 1 ) 2 * Y**2 * F(Y**2) dY.
%      (using Y = SQRT(X) )
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%  Modified:
%
%    15 October 2005
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

  weight(1:norder) = 2.0 * weight2(norder+2:2*norder+1) .* xtab(1:norder);
