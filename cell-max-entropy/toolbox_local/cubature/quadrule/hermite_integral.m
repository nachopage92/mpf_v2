function value = hermite_integral ( n )

%% HERMITE_INTEGRAL returns the value of a Hermite polynomial integral.
%
%  Discussion:
%
%    H(n) = Integral ( -Infinity < x < Infinity ) x^n exp(-x^2) dx
%
%    H(n) is 0 for n odd.
%
%    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
%
%  Modified:
%
%    25 July 2007
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the order of the integral.
%    0 <= N.
%
%    Output, real VALUE, the value of the integral.
%
  if ( n < 0 )

    value = - r8_huge ( );

  elseif ( mod ( n, 2 ) == 1 )

    value = 0.0;

  else

    value = i4_factorial2 ( n - 1 ) * sqrt ( pi ) / 2.0^( n / 2 );

  end

