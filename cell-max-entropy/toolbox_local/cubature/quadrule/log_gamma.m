function value = log_gamma ( x )

%% LOG_GAMMA calculates the natural logarithm of GAMMA(X).
%
%  Discussion:
%
%    The method uses Stirling's approximation, and is accurate to about
%    12 decimal places.
%
%  Modified:
%
%    12 October 2005
%
%  Reference:
%
%    Arthur Stroud and Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966.
%
%  Parameters:
%
%    Input, double X, the evaluation point.  The routine
%    will fail if GAMMA(X) is not positive.  X should be greater than 0.
%
%    Output, double VALUE, the natural logarithm of the
%    gamma function of X.
%
  if ( x < 0.5 )

    m = 1;
    x2 = 1.0 - x;

  else

    m = 0;
    x2 = x;

  end

  k = -1;

  while ( 1 )

    k = k + 1;

    if ( 6.0 < x2 + k )
      break
    end

  end

  z = x2 + k;

  y = ( z - 0.5 ) * log ( z ) - z + 0.9189385332047 + ...
       ( ( ( ( ( ...
       - 4146.0 / z^2 ...
       + 1820.0 ) / z^2 ...
       - 1287.0 ) / z^2 ...
       + 1716.0 ) / z^2 ...
       - 6006.0 ) / z^2 ...
       + 180180.0 ) / z / 2162160.0;

  if ( 0 < k )

    for i = 1 : k
      y = y - log ( x2 + k - i );
    end

  end

  if ( m ~= 0 )

    p = pi / sin ( pi * ( 1.0 - x2 ) );

    if ( p <= 0.0 )

      fprintf ( 1, '\n' );
      fprintf ( 1, 'LOG_GAMMA - Fatal error!\n' );

      error ( 'LOG_GAMMA - Fatal error!' );

    else

      y = log ( p ) - y;

    end

  end

  value = y;
