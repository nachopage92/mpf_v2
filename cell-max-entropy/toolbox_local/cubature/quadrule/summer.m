function result = summer ( func, norder, xtab, weight )

%% SUMMER carries out a quadrature rule over a single interval.
%
%  Discussion:
%
%    RESULT = sum ( 1 <= I <= NORDER ) WEIGHT(I) * FUNC ( XTAB(I) )
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
    fprintf ( 1, 'SUMMER - Fatal error!\n' );
    fprintf ( 1, '  NORDER must be at least 1.\n' );
    fprintf ( 1, '  The input value was NORDER = %d\n', norder );
    error ( 'SUMMER - Fatal error!' );
  end

  result = 0.0;
  for i = 1 : norder
    result = result + weight(i) * func ( xtab(i) );
  end
