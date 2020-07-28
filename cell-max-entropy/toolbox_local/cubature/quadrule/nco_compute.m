function [ xtab, weight ] = nco_compute ( order )

%% NCO_COMPUTE computes a Newton-Cotes open quadrature rule.
%
%  Discussion:
%
%    For the interval [X_MIN,X_MAX], the Newton-Cotes open quadrature rule
%    estimates
%
%      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
%
%    using ORDER equally spaced abscissas XTAB(I) and a weight vector
%    WEIGHT(I):
%
%      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
%
%    For the OPEN rule, the abscissas do not include the end points.
%
%  Modified:
%
%    26 May 2007
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer ORDER, the order of the rule.
%
%    Output, real XTAB(ORDER), the abscissas.
%
%    Output, real WEIGHT(ORDER), the weights.
%
  x_min = -1.0;
  x_max =  1.0;

  for i = 1 : order
    xtab(i) = ( ( order - i + 1 ) * x_min   ...
              + (         i     ) * x_max ) ...
              / ( order     + 1 );
  end

  weight = nc_compute ( order, x_min, x_max, xtab );
