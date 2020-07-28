function [ xtab, weight ] = legendre_set_x2_01 ( norder )

%% LEGENDRE_SET_X2_01 sets a Gauss-Legendre rule for X**2 * F(X) on [0,1].
%
%  Discussion:
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function w(x) = x * x.
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) X*X * F(X) dX
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

    xtab(1) =   0.75;

    weight(1) = 1.0 / 3.0;

  elseif ( norder == 2 )

    xtab(1) = 0.4558481560;
    xtab(2) = 0.8774851773;

    weight(1) = 0.1007858821;
    weight(2) = 0.2325474513;

  elseif ( norder == 3 )

    xtab(1) = 0.2949977901;
    xtab(2) = 0.6529962340;
    xtab(3) = 0.9270059759;

    weight(1) = 0.0299507030;
    weight(2) = 0.1462462693;
    weight(3) = 0.1571363611;

  elseif ( norder == 4 )

    xtab(1) = 0.2041485821;
    xtab(2) = 0.4829527049;
    xtab(3) = 0.7613992624;
    xtab(4) = 0.9514994506;

    weight(1) = 0.0103522408;
    weight(2) = 0.0686338872;
    weight(3) = 0.1434587898;
    weight(4) = 0.1108884156;

  elseif ( norder == 5 )

    xtab(1) = 0.1489457871;
    xtab(2) = 0.3656665274;
    xtab(3) = 0.6101136129;
    xtab(4) = 0.8265196792;
    xtab(5) = 0.9654210601;

    weight(1) = 0.0041138252;
    weight(2) = 0.0320556007;
    weight(3) = 0.0892001612;
    weight(4) = 0.1261989619;
    weight(5) = 0.0817647843;

  elseif ( norder == 6 )

    xtab(1) = 0.1131943838;
    xtab(2) = 0.2843188727;
    xtab(3) = 0.4909635868;
    xtab(4) = 0.6975630820;
    xtab(5) = 0.8684360583;
    xtab(6) = 0.9740954449;

    weight(1) = 0.0018310758;
    weight(2) = 0.0157202972;
    weight(3) = 0.0512895711;
    weight(4) = 0.0945771867;
    weight(5) = 0.1073764997;
    weight(6) = 0.0625387027;

  elseif ( norder == 7 )

    xtab(1) = 0.0888168334;
    xtab(2) = 0.2264827534;
    xtab(3) = 0.3999784867;
    xtab(4) = 0.5859978554;
    xtab(5) = 0.7594458740;
    xtab(6) = 0.8969109709;
    xtab(7) = 0.9798672262;

    weight(1) = 0.0008926880;
    weight(2) = 0.0081629256;
    weight(3) = 0.0294222113;
    weight(4) = 0.0631463787;
    weight(5) = 0.0917338033;
    weight(6) = 0.0906988246;
    weight(7) = 0.0492765018;

  elseif ( norder == 8 )

    xtab(1) = 0.0714910350;
    xtab(2) = 0.1842282964;
    xtab(3) = 0.3304477282;
    xtab(4) = 0.4944029218;
    xtab(5) = 0.6583480085;
    xtab(6) = 0.8045248315;
    xtab(7) = 0.9170993825;
    xtab(8) = 0.9839022404;

    weight(1) = 0.0004685178;
    weight(2) = 0.0044745217;
    weight(3) = 0.0172468638;
    weight(4) = 0.0408144264;
    weight(5) = 0.0684471834;
    weight(6) = 0.0852847692;
    weight(7) = 0.0768180933;
    weight(8) = 0.0397789578;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'LEGENDRE_SET_X2_01 - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    error ( 'LEGENDRE_SET_X2_01 - Fatal error!' );

  end
