function [ xtab, weight ] = legendre_set_x1_01 ( norder )

%% LEGENDRE_SET_X1_01 sets a Gauss-Legendre rule for X * F(X) on [0,1].
%
%  Discussion:
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function w(x) = x.
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) X * F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%  Modified:
%
%    13 October 2005
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

    xtab(1) =   0.6666666667;

    weight(1) = 0.5000000000;

  elseif ( norder == 2 )

    xtab(1) = 0.3550510257;
    xtab(2) = 0.8449489743;

    weight(1) = 0.1819586183;
    weight(2) = 0.3180413817;

  elseif ( norder == 3 )

    xtab(1) = 0.2123405382;
    xtab(2) = 0.5905331356;
    xtab(3) = 0.9114120405;

    weight(1) = 0.0698269799;
    weight(2) = 0.2292411064;
    weight(3) = 0.2009319137;

  elseif ( norder == 4 )

    xtab(1) = 0.1397598643;
    xtab(2) = 0.4164095676;
    xtab(3) = 0.7231569864;
    xtab(4) = 0.9428958039;

    weight(1) = 0.0311809710;
    weight(2) = 0.1298475476;
    weight(3) = 0.2034645680;
    weight(4) = 0.1355069134;

  elseif ( norder == 5 )

    xtab(1) = 0.0985350858;
    xtab(2) = 0.3045357266;
    xtab(3) = 0.5620251898;
    xtab(4) = 0.8019865821;
    xtab(5) = 0.9601901429;

    weight(1) = 0.0157479145;
    weight(2) = 0.0739088701;
    weight(3) = 0.1463888701;
    weight(4) = 0.1671746381;
    weight(5) = 0.0967815902;

  elseif ( norder == 6 )

    xtab(1) = 0.0730543287;
    xtab(2) = 0.2307661380;
    xtab(3) = 0.4413284812;
    xtab(4) = 0.6630153097;
    xtab(5) = 0.8519214003;
    xtab(6) = 0.9706835728;

    weight(1) = 0.0087383108;
    weight(2) = 0.0439551656;
    weight(3) = 0.0986611509;
    weight(4) = 0.1407925538;
    weight(5) = 0.1355424972;
    weight(6) = 0.0723103307;

  elseif ( norder == 7 )

    xtab(1) = 0.0562625605;
    xtab(2) = 0.1802406917;
    xtab(3) = 0.3526247171;
    xtab(4) = 0.5471536263;
    xtab(5) = 0.7342101772;
    xtab(6) = 0.8853209468;
    xtab(7) = 0.9775206136;

    weight(1) = 0.0052143622;
    weight(2) = 0.0274083567;
    weight(3) = 0.0663846965;
    weight(4) = 0.1071250657;
    weight(5) = 0.1273908973;
    weight(6) = 0.1105092582;
    weight(7) = 0.0559673634;

  elseif ( norder == 8 )

    xtab(1) = 0.0446339553;
    xtab(2) = 0.1443662570;
    xtab(3) = 0.2868247571;
    xtab(4) = 0.4548133152;
    xtab(5) = 0.6280678354;
    xtab(6) = 0.7856915206;
    xtab(7) = 0.9086763921;
    xtab(8) = 0.9822200849;

    weight(1) = 0.0032951914;
    weight(2) = 0.0178429027;
    weight(3) = 0.0454393195;
    weight(4) = 0.0791995995;
    weight(5) = 0.1060473594;
    weight(6) = 0.1125057995;
    weight(7) = 0.0911190236;
    weight(8) = 0.0445508044;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'LEGENDRE_SET_X1_01 - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    error ( 'LEGENDRE_SET_X1_01 - Fatal error!' );

  end
