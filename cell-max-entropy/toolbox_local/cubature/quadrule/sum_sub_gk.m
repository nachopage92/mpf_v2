function [ resultg, resultk, error_est ] = sum_sub_gk ( func, a, b, nsub, ...
  norderg, weightg, norderk, xtabk, weightk )

%% SUM_SUB_GK carries out a composite Gauss-Kronrod rule.
%
%  Discussion:
%
%    The integration interval is [ A, B ].
%
%    The integral to approximate:
%
%      Integral ( A <= X <= B ) F(X) dX
%
%    The quadrature rule:
%
%      H = ( B - A ) / NSUB
%      XMID(J) = A + 0.5 * H * ( 2 * J - 1 )
%
%      Sum ( 1 <= J <= NSUB )
%        Sum ( 1 <= I <= NORDERK )
%          WEIGHTK(I) * F ( XMID(J) + 0.5 * H * XTABK(I) )
%
%    The Gauss-Legendre weights should be computed by LEGCOM or LEGSET.
%    The Kronrod abscissas and weights should be computed by KRONSET.
%
%    The orders of the Gauss-Legendre and Kronrod rules must satisfy
%    NORDERK = 2 * NORDERG + 1.
%
%    The Kronrod rule uses the abscissas of the Gauss-Legendre rule,
%    plus more points, resulting in an efficient and higher order estimate.
%
%    The difference between the Gauss-Legendre and Kronrod estimates
%    is taken as an estimate of the error in the approximation to the
%    integral.
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
%    Input, integer NORDERG, the order of the Gauss-Legendre rule.
%    NORDERG must be at least 1.
%
%    Input, real WEIGHTG(NORDERG), the weights of the
%    Gauss-Legendre rule.
%
%    Input, integer NORDERK, the order of the Kronrod rule.
%    NORDERK must be at least 1.
%
%    Input, real XTABK(NORDERK), the abscissas of the
%    Kronrod rule.
%
%    Input, real WEIGHTK(NORDERK), the weights of the
%    Kronrod rule.
%
%    Output, real RESULTG, the approximate value of the
%    integral based on the Gauss-Legendre rule.
%
%    Output, real RESULTK, the approximate value of the
%    integral based on the Kronrod rule.
%
%    Output, real ERROR_EST, an estimate of the approximation
%    error.  This is computed by taking the sum of the absolute values of
%    the differences between the Gauss-Legendre and Kronrod rules
%    over each subinterval.  This is usually a good estimate of
%    the error in the value RESULTG.  The error in the Kronrod
%    estimate RESULTK is usually much smaller.
%
  resultg = 0.0;
  resultk = 0.0;
  error_est = 0.0;

  if ( a == b )
    return
  end

  if ( norderk ~= 2 * norderg + 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'SUM_SUB_GK - Fatal error!\n' );
    fprintf ( 1, '  NORDERK must equal 2 * NORDERG + 1.\n' );
    fprintf ( 1, '  The input value was NORDERG = %d\n', norderg );
    fprintf ( 1, '  The input value was NORDERK = %d\n', norderk );
    error (  'SUM_SUB_GK - Fatal error!' );
  end

  h = ( b - a ) / nsub;

  for j = 1 : nsub

    xmid = a + 0.5 * h * ( 2 * j - 1 );

    partg = 0.0;
    partk = 0.0;

    for i = 1 : norderk

      x = xmid + 0.5 * h * xtabk(i);
      fk = func ( x );
      partk = partk + 0.5 * h * weightk(i) * fk;

      if ( mod ( i, 2 ) == 0 )
        partg = partg + 0.5 * h * weightg(i/2) * fk;
      end

    end

    resultg = resultg + partg;
    resultk = resultk + partk;
    error_est = error_est + abs ( partk - partg );

  end
