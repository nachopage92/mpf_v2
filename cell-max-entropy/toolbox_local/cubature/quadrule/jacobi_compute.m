function [ xtab, weight ] = jacobi_compute ( norder, alpha, beta )

%% JACOBI_COM computes the abscissa and weights for Gauss-Jacobi quadrature.
%
%  Discussion:
%
%    The integration interval is [ -1, 1 ].
%
%    The weight function is w(x) = 1.0.
%
%    The integral to approximate is:
%
%      Integral ( -1 <= X <= 1 ) (1-X)**ALPHA * (1+X)**BETA * F(X) dX
%
%    The quadrature formula is:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    Thanks to Xu Xiang of Fudan University for pointing out that
%    an earlier implementation of this routine was incorrect!
%
%  Modified:
%
%    14 May 2007
%
%  Author:
%
%    Arthur Stroud, Don Secrest
%
%    MATLAB version by John Burkardt
%
%  Reference:
%
%    Arthur Stroud, Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966,
%    LC: QA299.4G3S7.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the quadrature rule to be computed.
%
%    Input, real ALPHA, BETA, the exponents of (1-X) and
%    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
%    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
%
%    Output, real XTAB(NORDER), the abscissas.
%
%    Output, real WEIGHT(NORDER), the weights.
%

%
%  Check ALPHA and BETA.
%
  if ( alpha <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'JACOBI_COMPUTE - Fatal error!\n' );
    fprintf ( 1, '  -1.0 < ALPHA is required.\n' );
    error ( 'JACOBI_COMPUTE - Fatal error!' );
  end

  if ( beta <= -1.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'JACOBI_COMPUTE - Fatal error!\n' );
    fprintf ( 1, '  -1.0 < BETA is required.\n' );
    error ( 'JACOBI_COMPUTE - Fatal error!' );
  end
%
%  Set the recursion coefficients.
%
  for i = 1 : norder

    if ( alpha + beta == 0.0 | beta - alpha == 0.0 )

      b(i) = 0.0;

    else

      b(i) = ( alpha + beta ) * ( beta - alpha ) / ...
            ( ( alpha + beta + 2 * i ) ...
            * ( alpha + beta + 2 * i - 2 ) );

    end

    if ( i == 1 )

      c(i) = 0.0;

    else

      c(i) = 4.0 * ( i - 1 ) * ( alpha + i - 1 ) * ( beta + i - 1 ) ...
        * ( alpha + beta + i - 1 ) / ( ( alpha + beta + 2 * i - 1 ) ...
        * ( alpha + beta + 2 * i - 2 )^2 * ( alpha + beta + 2 * i - 3 ) );

    end

  end

  delta = exp ( log_gamma ( alpha        + 1.0 ) ...
              + log_gamma (         beta + 1.0 ) ...
              - log_gamma ( alpha + beta + 2.0 ) );

  cc = delta * 2.0^( alpha + beta + 1.0 ) * prod ( c(2:norder) );

  for i = 1 : norder

    if ( i == 1 )

      an = alpha / norder;
      bn = beta / norder;

      r1 = ( 1.0 + alpha ) * ( 2.78 / ( 4.0 + norder * norder ) ...
        + 0.768 * an / norder );

      r2 = 1.0 + 1.48 * an + 0.96 * bn + 0.452 * an^2 + 0.83 * an * bn;

      x = ( r2 - r1 ) / r2;

    elseif ( i == 2 )

      r1 = ( 4.1 + alpha ) / ...
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( norder - 8.0 ) * ( 1.0 + 0.12 * alpha ) / norder;

      r3 = 1.0 + 0.012 * beta * ...
        ( 1.0 + 0.25 * abs ( alpha ) ) / norder;

      x = x - r1 * r2 * r3 * ( 1.0 - x );

    elseif ( i == 3 )

      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( norder - 8.0 ) / norder;

      r3 = 1.0 + 8.0 * beta / ( ( 6.28 + beta ) * norder * norder );

      x = x - r1 * r2 * r3 * ( xtab(1) - x );

    elseif ( i < norder - 1 )

      x = 3.0 * xtab(i-1) - 3.0 * xtab(i-2) + xtab(i-3);

    elseif ( i == norder - 1 )

      r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta );

      r2 = 1.0 / ( 1.0 + 0.639 * ( norder - 4.0 ) ...
        / ( 1.0 + 0.71 * ( norder - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * norder * norder ) );

      x = x + r1 * r2 * r3 * ( x - xtab(i-2) );

    elseif ( i == norder )

      r1 = ( 1.0 + 0.37 * beta ) / ( 1.67 + 0.28 * beta );

      r2 = 1.0 / ( 1.0 + 0.22 * ( norder - 8.0 ) / norder );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / ( ( 6.28 + alpha ) * norder * norder ) );

      x = x + r1 * r2 * r3 * ( x - xtab(i-2) );

    end

    [ x, dp2, p1 ] = jacobi_root ( x, norder, alpha, beta, b, c );

    xtab(i) = x;
    weight(i) = cc / ( dp2 * p1 );

  end
%
%  Reverse the order of the values.
%
  xtab = r8vec_reverse ( norder, xtab );
  weight = r8vec_reverse ( norder, weight );
