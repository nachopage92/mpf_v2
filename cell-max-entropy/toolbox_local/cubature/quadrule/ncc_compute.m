function [ xtab, weight ] = ncc_compute ( order )

%% NCC_COMPUTE computes a Newton-Cotes closed quadrature rule.
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
%      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) ).
%
%    For the CLOSED rule, the abscissas include the end points.
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
%    Output, real XTAB(ORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(ORDER), the weights of the rule.
%

%
%  Compute a closed quadrature rule.
%
  x_min = -1.0;
  x_max =  1.0;

  if ( order == 1 )

    xtab(1) = ( x_max + x_min ) / 2.0;
    weight(1) = x_max - x_min;

  else

    for i = 1 : order
      xtab(i) = ( ( order - i     ) * x_min   ...
                + (         i - 1 ) * x_max ) ...
                / ( order     - 1 );
    end

    weight = nc_compute ( order, x_min, x_max, xtab );

  end


