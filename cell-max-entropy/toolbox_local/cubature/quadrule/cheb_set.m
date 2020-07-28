function [ xtab, weight ] = cheb_set ( norder )

%% CHEB_SET sets abscissas and weights for Chebyshev quadrature.
%
%  Discussion:
%
%    The integration interval is [ -1, 1 ].
%
%    The weight function w(x) = 1.0.
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    The Chebyshev rule is distinguished by using equal weights.
%
%  Modified:
%
%    12 October 2005
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Milton Abramowitz and Irene Stegun,
%    Handbook of Mathematical Functions,
%    National Bureau of Standards, 1964.
%
%    Hermann Engels,
%    Numerical Quadrature and Cubature,
%    Academic Press, 1980.
%
%    Zdenek Kopal,
%    Numerical Analysis,
%    John Wiley, 1955.
%
%    Daniel Zwillinger, editor,
%    Standard Mathematical Tables and Formulae,
%    30th Edition,
%    CRC Press, 1996.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule.
%    NORDER may only have the values 1, 2, 3, 4, 5, 6, 7 or 9.
%    There are NO other Chebyshev rules with real abscissas.
%
%    Output, real XTAB(NORDER), the abscissas of the rule,
%    which are symmetric in [-1,1].
%
%    Output, real WEIGHT(NORDER), the weights of the rule,
%    which should each equal 2 / NORDER.
%
  if ( norder == 1 )

    xtab(1) = 0.0;

  elseif ( norder == 2 )

    xtab(1) = - 1.0 / sqrt ( 3.0 );
    xtab(2) =   1.0 / sqrt ( 3.0 );

  elseif ( norder == 3 )

    xtab(1) = - 1.0 / sqrt ( 2.0 );
    xtab(2) =   0.0;
    xtab(3) =   1.0 / sqrt ( 2.0 );

  elseif ( norder == 4 )

    xtab(1) =   - sqrt ( ( 1.0 + 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab(2) =   - sqrt ( ( 1.0 - 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab(3) =     sqrt ( ( 1.0 - 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab(4) =     sqrt ( ( 1.0 + 2.0/ sqrt ( 5.0 ) ) / 3.0 );

  elseif ( norder == 5 )

    xtab(1) = - sqrt ( ( 5.0 + sqrt ( 11.0) ) / 12.0 );
    xtab(2) = - sqrt ( ( 5.0 - sqrt ( 11.0) ) / 12.0 );
    xtab(3) =   0.0;
    xtab(4) =   sqrt ( ( 5.0 - sqrt ( 11.0) ) / 12.0 );
    xtab(5) =   sqrt ( ( 5.0 + sqrt ( 11.0) ) / 12.0 );

  elseif ( norder == 6 )

    xtab(1) = - 0.866246818107820591383598;
    xtab(2) = - 0.422518653761111529118546;
    xtab(3) = - 0.266635401516704720331534;
    xtab(4) =   0.266635401516704720331534;
    xtab(5) =   0.422518653761111529118546;
    xtab(6) =   0.866246818107820591383598;

  elseif ( norder == 7 )

    xtab(1) = - 0.883861700758049035704224;
    xtab(2) = - 0.529656775285156811385048;
    xtab(3) = - 0.323911810519907637519673;
    xtab(4) =   0.0;
    xtab(5) =   0.323911810519907637519673;
    xtab(6) =   0.529656775285156811385048;
    xtab(7) =   0.883861700758049035704224;

  elseif ( norder == 9 )

    xtab(1) = - 0.911589307728434473664949;
    xtab(2) = - 0.601018655380238071428128;
    xtab(3) = - 0.528761783057879993260181;
    xtab(4) = - 0.167906184214803943068031;
    xtab(5) =   0.0;
    xtab(6) =   0.167906184214803943068031;
    xtab(7) =   0.528761783057879993260181;
    xtab(8) =   0.601018655380238071428128;
    xtab(9) =   0.911589307728434473664949;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'CHEB_SET - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    fprintf ( 1, '  Legal values are 1 through 7, and 9.\n' );
    error ( 'CHEB_SET - Fatal error!' );

  end

  weight(1:norder) = 2.0 / norder;
