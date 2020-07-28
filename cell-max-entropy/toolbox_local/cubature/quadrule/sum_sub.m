function result = sum_sub ( func, a, b, nsub, norder, xlo, xhi, ...
  xtab, weight )

%% SUM_SUB carries out a composite quadrature rule.
%
%  Discussion:
%
%    SUM_SUB assumes the original rule was written for [XLO,XHI].
%
%    The integration interval is [A,B].
%
%    The integral to approximate:
%
%      Integral ( A <= X <= B ) F(X) dX
%
%  Modified:
%
%    15 October 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, external FUNC, the name of the FORTRAN function which
%    evaluates the integrand.  The function must have the form
%      function func ( x ).
%
%    Input, real A, B, the lower and upper limits of integration.
%
%    Input, integer NSUB, the number of equal subintervals into
%    which the finite interval (A,B) is to be subdivided for
%    higher accuracy.  NSUB must be at least 1.
%
%    Input, integer NORDER, the order of the rule.
%    NORDER must be at least 1.
%
%    Input, real XLO, XHI, the left and right endpoints of
%    the interval over which the quadrature rule was defined.
%
%    Input, real XTAB(NORDER), the abscissas of a quadrature
%    rule for the interval [XLO,XHI].
%
%    Input, real WEIGHT(NORDER), the weights of the
%    quadrature rule.
%
%    Output, real RESULT, the approximate value of the integral.
%
  if ( a == b )
    result = 0.0;
    return
  end

  if ( norder < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SUM_SUB - Fatal error!\n' );
    fprintf ( 1, '  Nonpositive value of NORDER = %d\n', norder );
    error ( 'SUM_SUB - Fatal error!' );
  end

  if ( nsub < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SUM_SUB - Fatal error!\n' );
    fprintf ( 1, '  Nonpositive value of NSUB = %d\n', nsub );
    error ( 'SUM_SUB - Fatal error!' );
  end

  if ( xlo == xhi )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SUM_SUB - Fatal error!\n' );
    fprintf ( 1, '  XLO = XHI.\n' );
    error ( 'SUM_SUB - Fatal error!' );
  end

  volume = 0.0;
  result = 0.0;

  for j = 1 : nsub

    a_sub = ( ( nsub - j + 1 ) * a   ...
            + (        j - 1 ) * b ) ...
            / ( nsub         );

    b_sub = ( ( nsub - j )     * a   ...
            + (        j )     * b ) ...
            / ( nsub     );

    quad_sub = 0.0;
    for i = 1 : norder
      x = ( ( xhi - xtab(i)       ) * a_sub   ...
          + (       xtab(i) - xlo ) * b_sub ) ...
          / ( xhi           - xlo );
      quad_sub = quad_sub + weight(i) * func ( x );
    end

    volume_sub = ( b - a ) / ( ( xhi - xlo ) * nsub );
    result_sub = quad_sub * volume_sub;

    volume = volume + volume_sub;
    result = result + result_sub;

  end
