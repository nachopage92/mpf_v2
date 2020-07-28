function result = laguerre_sum ( func, a, norder, xtab, weight )

%% LAGUERRE_SUM carries out Laguerre quadrature over [ A, +Infinity ).
%
%  Discussion:
%
%    The simplest Laguerre integral to approximate is the
%    integral from 0 to INFINITY of EXP(-X) * F(X).  When this is so,
%    it is easy to modify the rule to approximate the integral from
%    A to INFINITY as well.
%
%    Another common Laguerre integral to approximate is the
%    integral from 0 to Infinity of EXP(-X) * X**ALPHA * F(X).
%    This routine may be used to sum up the terms of the Laguerre
%    rule for such an integral as well.  However, if ALPHA is nonzero,
%    then there is no simple way to extend the rule to approximate the
%    integral from A to INFINITY.  The simplest procedures would be
%    to approximate the integral from 0 to A.
%
%    The integration interval is [ A, +Infinity ) or [ 0, +Infinity ).
%
%    The weight function is w(x) = EXP ( - X ) or EXP ( - X ) * X**ALPHA
%
%    The integral to approximate:
%
%      Integral ( A <= X <= +Infinity ) EXP ( -X ) * F(X) dX 
%    or
%      Integral ( 0 <= X <= +Infinity ) EXP ( -X ) * X**ALPHA * F(X) dX
%
%    The quadrature rule:
%
%      EXP ( - A ) * Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) + A )
%
%    or
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
%    Abramowitz, Stegun,
%    Handbook of Mathematical Functions,
%    National Bureau of Standards, 1964.
%
%    Daniel Zwillinger, editor,
%    Standard Mathematical Tables and Formulae,
%    30th Edition,
%    CRC Press, 1996.
%
%  Parameters:
%
%    Input, external FUNC, the name of the FORTRAN function which
%    evaluates the integrand.  The function must have the form
%      function func ( x ).
%
%    Input, real A, the beginning of the integration interval.
%
%    Input, integer NORDER, the order of the rule.
%
%    Input, real XTAB(NORDER), the abscissas of the rule.
%
%    Input, real WEIGHT(NORDER), the weights of the rule.
%
%    Output, real RESULT, the approximate value of the integral.
%
  if ( norder < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'LAGUERRE_SUM - Fatal error!\n' );
    fprintf ( 1, '  Nonpositive NORDER = %d\n', norder );
    error ( 'LAGUERRE_SUM - Fatal error!' );
  end

  result = 0.0;
  for i = 1 : norder
    result = result + weight(i) * func ( xtab(i) + a );
  end

  result = exp ( - a ) * result;
