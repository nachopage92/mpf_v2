function [ xtab, weight ] = ncoh_compute ( order )

%% NCOH_COMPUTE computes a Newton-Cotes "open half" quadrature rule.
%
%  Discussion:
%
%    For the interval [X_MIN,X_MAX], the Newton-Cotes "open half" quadrature 
%    rule estimates
%
%      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
%
%    using ORDER equally spaced abscissas XTAB(I), each of which is
%    the midpoint of one of ORDER equal subintervals,
%    and a weight vector WEIGHT(I):
%
%      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) ).
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
    xtab(i) = ( ( 2 * order - 2 * i + 1 ) * x_min   ...
              + (             2 * i - 1 ) * x_max ) ...
              / ( 2 * order             );
  end

  weight = nc_compute ( order, x_min, x_max, xtab );
