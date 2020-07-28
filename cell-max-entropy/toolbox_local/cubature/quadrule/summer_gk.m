function [ resultg, resultk ] = summer_gk ( func, norderg, weightg, norderk, ...
  xtabk, weightk )

%% SUMMER_GK carries out Gauss-Kronrod quadrature over a single interval.
%
%  Discussion:
%
%    The abscissas for the Gauss-Legendre rule of order NORDERG are
%    not required, since they are assumed to be the even-indexed
%    entries of the corresponding Kronrod rule.
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
%    Input, integer NORDERG, the order of the Gauss-Legendre rule.
%
%    Input, real WEIGHTG(NORDERG), the weights of the
%    Gauss-Legendre rule.
%
%    Input, integer NORDERK, the order of the Kronrod rule.  NORDERK
%    must equal 2 * NORDERG + 1.
%
%    Input, real XTABK(NORDERK), the abscissas of the Kronrod rule.
%
%    Input, real WEIGHTK(NORDERK), the weights of the Kronrod rule.
%
%    Output, real RESULTG, the approximate value of the
%    integral, based on the Gauss-Legendre rule.
%
%    Output, real RESULTK, the approximate value of the integral,
%    based on the Kronrod rule.
%
  if ( norderk ~= 2 * norderg + 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SUMMER_GK - Fatal error!\n' );
    fprintf ( 1, '  NORDERK must equal 2 * NORDERG + 1.\n' );
    fprintf ( 1, '  The input value was NORDERG = %d\n', norderg );
    fprintf ( 1, '  The input value was NORDERK = %d\n', norderk );
    error ( 'SUMMER_GK - Fatal error!' );
  end

  resultg = 0.0;
  resultk = 0.0;

  for i = 1 : norderk

    fk = func ( xtabk(i) );

    resultk = resultk + weightk(i) * fk;

    if ( mod ( i, 2 ) == 0 )
      resultg = resultg + weightg(i/2) * fk;
    end

  end
