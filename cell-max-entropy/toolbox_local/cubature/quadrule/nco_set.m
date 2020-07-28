function [ xtab, weight ] = nco_set ( norder )

%% NCO_SET sets abscissas and weights for open Newton-Cotes quadrature.
%
%  Discussion:
%
%    The open Newton-Cotes rules use equally spaced abscissas, and
%    hence may be used with equally spaced data.
%
%    The rules are called "open" because they do not include the interval
%    endpoints.
%
%    Most of the rules involve negative weights.  These can produce loss
%    of accuracy due to the subtraction of large, nearly equal quantities.
%
%    The integration interval is [ -1, 1 ].
%
%    The weight function is w(x) = 1.0;
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    In Mathematica, the open Newton-Cotes weights and abscissas
%    can be computed by the commands:
%
%      << NumericalMath`NewtonCotes`
%      NewtonCotesWeights [ order, -1, 1, QuadratureType -> Open ]
%
%  Modified:
%
%    02 May 2006
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Abramowitz and Stegun,
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
%    Input, integer NORDER, the order of the rule.
%    NORDER must be between 1 and 7, and 9.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%    The weights are rational, symmetric, and should sum to 2.
%    Some weights may be negative.
%
  if ( norder == 1 )

    weight(1) = 2.0;

  elseif ( norder == 2 )

    weight(1) = 1.0;
    weight(2) = 1.0;

  elseif ( norder == 3 )

    d = 3.0;

    weight(1) =   4.0 / d;
    weight(2) = - 2.0 / d;
    weight(3) =   4.0 / d;

  elseif ( norder == 4 )

    d = 12.0;

    weight(1) = 11.0 / d;
    weight(2) =  1.0 / d;
    weight(3) =  1.0 / d;
    weight(4) = 11.0 / d;

  elseif ( norder == 5 )

    d = 10.0;

    weight(1) =   11.0 / d;
    weight(2) = - 14.0 / d;
    weight(3) =   26.0 / d;
    weight(4) = - 14.0 / d;
    weight(5) =   11.0 / d;

  elseif ( norder == 6 )

    d = 1440.0;

    weight(1) =  1222.0 / d;
    weight(2) = - 906.0 / d;
    weight(3) =  1124.0 / d;
    weight(4) =  1124.0 / d;
    weight(5) = - 906.0 / d;
    weight(6) =  1222.0 / d;

  elseif ( norder == 7 )

    d = 945.0;

    weight(1) =    920.0 / d;
    weight(2) = - 1908.0 / d;
    weight(3) =   4392.0 / d;
    weight(4) = - 4918.0 / d;
    weight(5) =   4392.0 / d;
    weight(6) = - 1908.0 / d;
    weight(7) =    920.0 / d;

  elseif ( norder == 9 )

    d = 4536.0;

    weight(1) =    4045.0 / d;
    weight(2) = - 11690.0 / d;
    weight(3) =   33340.0 / d;
    weight(4) = - 55070.0 / d;
    weight(5) =   67822.0 / d;
    weight(6) = - 55070.0 / d;
    weight(7) =   33340.0 / d;
    weight(8) = - 11690.0 / d;
    weight(9) =    4045.0 / d;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'NCO_SET - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    fprintf ( 1, '  Legal values are 1 to 7, and 9.\n' );
    error ( 'NCO_SET - Fatal error!' );

  end
%
%  Set the abscissas.
%
  for i = 1 : norder
    xtab(i) = ( 2 * i - norder - 1 ) / ( norder + 1 );
  end
