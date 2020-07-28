function [ xtab, weight ] = legendre_set_x0_01 ( norder )

%% LEGENDRE_SET_X0_01 sets a Gauss-Legendre rule for F(X) on [0,1].
%
%  Discussion:
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function w(x) = 1.0.
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) F(X) dX
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
%    John Burkardt
%
%  Reference:
%
%    Abramowitz, Stegun,
%    Handbook of Mathematical Functions,
%    National Bureau of Standards, 1964, page 921.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule.
%    NORDER must be between 1 and 8.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%
  if ( norder == 1 )

    xtab(1) =   0.5;

    weight(1) = 1.0;

  elseif ( norder == 2 )

    xtab(1) = 0.2113248654;
    xtab(2) = 0.7886751346;

    weight(1) = 0.5;
    weight(2) = 0.5;

  elseif ( norder == 3 )

    xtab(1) = 0.1127016654;
    xtab(2) = 0.5000000000;
    xtab(3) = 0.8872983346;

    weight(1) = 5.0 / 18.0;
    weight(2) = 8.0 / 18.0;
    weight(3) = 5.0 / 18.0;

  elseif ( norder == 4 )

    xtab(1) = 0.0694318442;
    xtab(2) = 0.3300094782;
    xtab(3) = 0.6699905218;
    xtab(4) = 0.9305681558;

    weight(1) = 0.1739274226;
    weight(2) = 0.3260725774;
    weight(3) = 0.3260725774;
    weight(4) = 0.1739274226;

  elseif ( norder == 5 )

    xtab(1) = 0.0469100770;
    xtab(2) = 0.2307653449;
    xtab(3) = 0.5000000000;
    xtab(4) = 0.7692346551;
    xtab(5) = 0.9530899230;

    weight(1) = 0.1184634425;
    weight(2) = 0.2393143352;
    weight(3) = 0.2844444444;
    weight(4) = 0.2393143352;
    weight(5) = 0.1184634425;

  elseif ( norder == 6 )

    xtab(1) = 0.0337652429;
    xtab(2) = 0.1693953068;
    xtab(3) = 0.3806904070;
    xtab(4) = 0.6193095930;
    xtab(5) = 0.8306046932;
    xtab(6) = 0.9662347571;

    weight(1) = 0.0856622462;
    weight(2) = 0.1803807865;
    weight(3) = 0.2339569673;
    weight(4) = 0.2339569673;
    weight(5) = 0.1803807865;
    weight(6) = 0.0856622462;

  elseif ( norder == 7 )

    xtab(1) = 0.0254460438;
    xtab(2) = 0.1292344072;
    xtab(3) = 0.2970774243;
    xtab(4) = 0.5000000000;
    xtab(5) = 0.7029225757;
    xtab(6) = 0.8707655928;
    xtab(7) = 0.9745539562;

    weight(1) = 0.0647424831;
    weight(2) = 0.1398526957;
    weight(3) = 0.1909150253;
    weight(4) = 0.2089795918;
    weight(5) = 0.1909150253;
    weight(6) = 0.1398526957;
    weight(7) = 0.0647424831;

  elseif ( norder == 8 )

    xtab(1) = 0.0198550718;
    xtab(2) = 0.1016667613;
    xtab(3) = 0.2372337950;
    xtab(4) = 0.4082826788;
    xtab(5) = 0.5917173212;
    xtab(6) = 0.7627662050;
    xtab(7) = 0.8983332387;
    xtab(8) = 0.9801449282;

    weight(1) = 0.0506142681;
    weight(2) = 0.1111905172;
    weight(3) = 0.1568533229;
    weight(4) = 0.1813418917;
    weight(5) = 0.1813418917;
    weight(6) = 0.1568533229;
    weight(7) = 0.1111905172;
    weight(8) = 0.0506142681;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'LEGENDRE_SET_X0_01 - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    error ( 'LEGENDRE_SET_X0_01 - Fatal error!' );

  end

