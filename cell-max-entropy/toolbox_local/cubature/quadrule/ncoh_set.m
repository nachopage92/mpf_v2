function [ xtab, weight ] = ncoh_set ( order )

%% NCOH_SET sets abscissas and weights for "open half" Newton-Cotes quadrature.
%
%  Discussion:
%
%    The open Newton-Cotes rules use equally spaced abscissas, and
%    hence may be used with equally spaced data.
%
%    The rules are called "open" because the abscissas do not include
%    the interval endpoints.
%
%    For this uncommon type of open Newton-Cotes rule, the abscissas for
%    rule N are found by dividing the interval into N equal subintervals,
%    and using the midpoint of each subinterval as the abscissa.
%
%    Most of the rules involve negative weights.  These can produce loss
%    of accuracy due to the subtraction of large, nearly equal quantities.
%
%    The integration interval is [ -1, 1 ].
%
%    The weight function is w(x) = 1.0.
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    In Mathematica, the "open half" Newton-Cotes weights and abscissas
%    can be computed by the commands:
%
%      << NumericalMath`NewtonCotes`
%      NewtonCotesWeights [ order, -1, 1, QuadratureType -> Open ]
%
%  Modified:
%
%    29 August 2006
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Milton Abramowitz, Irene Stegun,
%    Handbook of Mathematical Functions,
%    National Bureau of Standards, 1964.
%
%    Stephen Wolfram,
%    The Mathematica Book,
%    Fourth Edition,
%    Wolfram Media / Cambridge University Press, 1999.
%
%    Daniel Zwillinger, editor,
%    Standard Mathematical Tables and Formulae,
%    30th Edition,
%    CRC Press, 1996.
%
%  Parameters:
%
%    Input, integer ORDER, the order of the rule.
%    ORDER must be between 1 and 10.
%
%    Output, real XTAB(ORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(ORDER), the weights of the rule.
%    The weights are rational, symmetric, and should sum to 2.
%    Some weights may be negative.
%
  if ( order == 1 )

    weight(1) = 2.0;

  elseif ( order == 2 )

    weight(1) = 1.0;
    weight(2) = 1.0;

  elseif ( order == 3 )

    d = 4.0;

    weight(1) =   3.0 / d;
    weight(2) =   2.0 / d;
    weight(3) =   3.0 / d;

  elseif ( order == 4 )

    d = 24.0;

    weight(1) = 13.0 / d;
    weight(2) = 11.0 / d;
    weight(3) = 11.0 / d;
    weight(4) = 13.0 / d;

  elseif ( order == 5 )

    d = 576.0;

    weight(1) =  275.0 / d;
    weight(2) =  100.0 / d;
    weight(3) =  402.0 / d;
    weight(4) =  100.0 / d;
    weight(5) =  275.0 / d;

  elseif ( order == 6 )

    d = 640.0;

    weight(1) =   247.0 / d;
    weight(2) =   139.0 / d;
    weight(3) =   254.0 / d;
    weight(4) =   254.0 / d;
    weight(5) =   139.0 / d;
    weight(6) =   247.0 / d;

  elseif ( order == 7 )

    d = 138240.0;

    weight(1) =   49490.0 / d;
    weight(2) =    1764.0 / d;
    weight(3) =  112014.0 / d;
    weight(4) =  -50056.0 / d;
    weight(5) =  112014.0 / d;
    weight(6) =    1764.0 / d;
    weight(7) =   49490.0 / d;

  elseif ( order == 8 )

    d = 967680.0;

    weight(1) =  295627.0 / d;
    weight(2) =   71329.0 / d;
    weight(3) =  471771.0 / d;
    weight(4) =  128953.0 / d;
    weight(5) =  128953.0 / d;
    weight(6) =  471771.0 / d;
    weight(7) =   71329.0 / d;
    weight(8) =  295627.0 / d;

  elseif ( order == 9 )

    d = 2867200.0;

    weight(1) =    832221.0 / d;
    weight(2) =   -260808.0 / d;
    weight(3) =   2903148.0 / d;
    weight(4) =  -3227256.0 / d;
    weight(5) =   5239790.0 / d;
    weight(6) =  -3227256.0 / d;
    weight(7) =   2903148.0 / d;
    weight(8) =   -260808.0 / d;
    weight(9) =    832221.0 / d;

  elseif ( order == 10 )

    d = 18579456.0;

    weight(1) =    4751285.0 / d;
    weight(2) =    -492755.0 / d;
    weight(3) =   12269956.0 / d;
    weight(4) =   -6274220.0 / d;
    weight(5) =    8325190.0 / d;
    weight(6) =    8325190.0 / d;
    weight(7) =   -6274220.0 / d;
    weight(8) =   12269956.0 / d;
    weight(9) =    -492755.0 / d;
    weight(10) =   4751285.0 / d;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'NCOH_SET - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of ORDER = %d\n', order );
    fprintf ( 1, '  Legal values are 1 to 10.\n' );
    error ( 'NCOH_SET - Fatal error!' );

  end
%
%  Set the abscissas.
%
  a = -1.0;
  b = +1.0;

  for i = 1 : order
    xtab(i) = ( ( 2 * order - 2 * i + 1 ) * a   ...
              + (             2 * i - 1 ) * b ) ...
              / ( 2 * order             );
  end
