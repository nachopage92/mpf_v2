function [ xtab, weight ] = hermite_compute ( norder )

%% HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
%
%  Discussion:
%
%    The abscissas are the zeros of the N-th order Hermite polynomial.
%
%    The integration interval is ( -Infinity, +Infinity ).
%
%    The weight function is w(x) = exp ( - x**2 );
%
%    The integral to approximate:
%
%      Integral ( -INFINITY < X < +INFINITY ) exp ( - X**2 ) * F(X) dX
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
%    Input, integer NORDER, the order of the formula to be computed.
%
%    Output, real XTAB(NORDER), the Gauss-Hermite abscissas.
%
%    Output, real WEIGHT(NORDER), the Gauss-Hermite weights.
%
  cc = 1.7724538509 * gamma ( norder ) / ( 2.0^( norder - 1 ) );

  s = ( 2.0 * norder + 1.0 )^( 1.0 / 6.0 );

  for i = 1 : floor ( ( norder + 1 ) / 2 )

    if ( i == 1 )

      x = s^3 - 1.85575 / s;

    elseif ( i == 2 )

      x = x - 1.14 * ( ( norder )^0.426 ) / x;

    elseif ( i == 3 )

      x = 1.86 * x - 0.86 * xtab(1);

    elseif ( i == 4 )

      x = 1.91 * x - 0.91 * xtab(2);

    else

      x = 2.0 * x - xtab(i-2);

    end

    [ x, dp2, p1 ] = hermite_root ( x, norder );

    xtab(i) = x;
    weight(i) = ( cc / dp2 ) / p1;

    xtab(norder-i+1) = - x;
    weight(norder-i+1) = weight(i);

  end
%
%  Reverse the order of the XTAB values.
%
  xtab  = r8vec_reverse ( norder, xtab )';
  weight= weight';
