function [ xtab, weight ] = cheb_tc_set ( norder )

%% CHEB_TC_SET sets up closed Gauss-Chebyshev (first kind) quadrature.
%
%  Discussion:
%
%    The integration interval is [ -1, 1 ].
%
%    The weight function is w(x) = 1 / SQRT ( 1 - X**2 )
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) F(X) / SQRT ( 1 - X**2 ) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    If NORDER points are used, then Gauss-Chebyshev quadrature
%    will compute the integral exactly, whenever F(X) is a polynomial
%    of degree 2*NORDER-3 or less.
%
%    The abscissas include -1 and 1.
%
%    If the order is doubled, the abscissas of the new rule include
%    all the points of the old rule.  This fact can be used to
%    efficiently implement error estimation.
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
%    Input, int NORDER, the order of the rule, which must be at least 2.
%
%    Output, double XTAB(NORDER), the abscissas of the rule.
%
%    Output, double WEIGHT(NORDER), the weights of the rule.
%    The first and last weights are 0.5 * PI / ( NORDER - 1),
%    and all other weights are PI / ( NORDER - 1 ).
%
  if ( norder < 2 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'CHEB_TC_SET - Fatal error!\n' );
    fprintf ( 1, '  NORDER must be at least 2.\n' );
    fprintf ( 1, '  The input value was NORDER = %d\n', norder );
    error ( 'CHEB_TC_SET - Fatal error!' );
  end

  for i = 1 : norder

    angle = ( i - 1 ) * pi / ( norder - 1 );
    xtab(i) = cos ( angle );

  end

  weight(1) = pi / ( 2 * ( norder - 1 ) );
  weight(2:norder-1) = pi / ( norder - 1 );
  weight(norder) = pi / real ( 2 * ( norder - 1 ) );
