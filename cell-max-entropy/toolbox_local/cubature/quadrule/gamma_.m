function value = gamma ( x )

%% GAMMA computes the gamma function using Hastings's approximation.
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
%    Input, real X, the argument at which the gamma function
%    is to be evaluated.  X must be greater than 0, and less than 70.
%
%    Output, real GAMMA, the gamma function at X.
%
  if ( x <= 0.0 )
    gamma = 0.0;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GAMMA - Fatal error!\n' );
    fprintf ( 1, '  Input argument X <= 0.\n' );
    error ( 'GAMMA - Fatal error!' );
  end

  if ( 70.0 <= x )
    gamma = Inf;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'GAMMA - Fatal error!\n' );
    fprintf ( 1, '  Input argument 70 <= X.\n' );
    error ( 'GAMMA - Fatal error!' );
  end

  if ( x == 1.0 )
    value = 1.0;
    return
  end

  if ( x <= 1.0 )
    value = gamma_hastings ( x ) / x;
    return
  end

  z = x;

  za = 1.0;

  while ( 1 )

    z = z - 1.0;

    if ( z < 1.0 )
      value = za * gamma_hastings ( z );
      break;
    elseif ( z == 1.0 )
      value = za;
      break;
    end

    za = za * z;

  end
