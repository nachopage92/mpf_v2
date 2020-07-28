function [ alpha, beta, gamma ] = bdf_set ( norder )

%% BDF_SET sets weights for backward differentiation ODE weights.
%
%  Discussion:
%
%    GAMMA * Y(N+1) = Sum ( 1 <= I <= NORDER ) ALPHA(I) * Y(N+1-I)
%                     + dX * BETA * Y'(X(N+1),Y(N+1))
%
%    This is equivalent to the backward differentiation corrector formulas.
%
%  Modified:
%
%    12 October 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule, between 1 and 6.
%
%    Output, real ALPHA(NORDER), BETA, GAMMA, the weights.
%
  if ( norder == 1 )
    beta =     1.0;
    gamma =    1.0;
    alpha(1) = 1.0;
  elseif ( norder == 2 )
    beta =       2.0;
    gamma =      3.0;
    alpha(1) =   4.0;
    alpha(2) = - 1.0;
  elseif ( norder == 3 )
    beta =       6.0;
    gamma =     11.0;
    alpha(1) =  18.0;
    alpha(2) = - 9.0;
    alpha(3) =   2.0;
  elseif ( norder == 4 )
    beta =       12.0;
    gamma =      25.0;
    alpha(1) =   48.0;
    alpha(2) = - 36.0;
    alpha(3) =   16.0;
    alpha(4) =  - 3.0;
  elseif ( norder == 5 )
    beta =        60.0;
    gamma =      137.0;
    alpha(1) =   300.0;
    alpha(2) = - 300.0;
    alpha(3) =   200.0;
    alpha(4) =  - 75.0;
    alpha(5) =    12.0;
  elseif ( norder == 6 )
    beta =        60.0;
    gamma =      147.0;
    alpha(1) =   360.0;
    alpha(2) = - 450.0;
    alpha(3) =   400.0;
    alpha(4) = - 225.0;
    alpha(5) =    72.0;
    alpha(6) =  - 10.0;
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, 'BDF_SET - Fatal error!\n' );
    fprintf ( 1, '  Illegal order requested = %d\n', norder );
    error ( 'BDF_SET - Fatal error!' );
  end
