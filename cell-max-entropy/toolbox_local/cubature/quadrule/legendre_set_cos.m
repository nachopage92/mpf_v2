function [ xtab, weight ] = legendre_set_cos ( norder )

%% LEGENDRE_SET_COS sets a Gauss-Legendre rule for COS(X) * F(X) on [-PI/2,PI/2].
%
%  Discussion:
%
%    The integration interval is [ -PI/2, PI/2 ].
%
%    The weight function w(x) = COS(X).
%
%    The integral to approximate:
%
%      Integral ( -PI/2 <= X <= PI/2 ) COS(X) * F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    The same rule can be used to approximate
%
%      Integral ( 0 <= X <= PI ) SIN(X) * F(X) dX
%
%    as
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) + PI/2 )
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
%    Gwynne Evans,
%    Practical Numerical Integration,
%    Wiley, 1993, QA299.3E93, page 310.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule.
%    NORDER must be between 1, 2, 4, 8 or 16.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%
  if ( norder == 1 )

    xtab(1) = 0.0;

    weight(1) = 2.0;

  elseif ( norder == 2 )

    xtab(1) = - 0.68366739008990304094;
    xtab(2) =   0.68366739008990304094;

    weight(1) = 1.0;
    weight(2) = 1.0;

  elseif ( norder == 4 )

    xtab(1) = - 1.1906765638948557415;
    xtab(2) = - 0.43928746686001514756;
    xtab(3) =   0.43928746686001514756;
    xtab(4) =   1.1906765638948557415;

    weight(1) = 0.22407061812762016065;
    weight(2) = 0.77592938187237983935;
    weight(3) = 0.77592938187237983935;
    weight(4) = 0.22407061812762016065;

  elseif ( norder == 8 )

    xtab(1) = - 1.4414905401823575701;
    xtab(2) = - 1.1537256454567275850;
    xtab(3) = - 0.74346864787549244989;
    xtab(4) = - 0.25649650741623123020;
    xtab(5) =   0.25649650741623123020;
    xtab(6) =   0.74346864787549244989;
    xtab(7) =   1.1537256454567275850;
    xtab(8) =   1.4414905401823575701;

    weight(1) = 0.027535633513767011149;
    weight(2) = 0.14420409203022750950;
    weight(3) = 0.33626447785280459621;
    weight(4) = 0.49199579660320088314;
    weight(5) = 0.49199579660320088314;
    weight(6) = 0.33626447785280459621;
    weight(7) = 0.14420409203022750950;
    weight(8) = 0.027535633513767011149;

  elseif ( norder == 16 )

    xtab( 1) = - 1.5327507132362304779;
    xtab( 2) = - 1.4446014873666514608;
    xtab( 3) = - 1.3097818904452936698;
    xtab( 4) = - 1.1330068786005003695;
    xtab( 5) = - 0.92027786206637096497;
    xtab( 6) = - 0.67861108097560545347;
    xtab( 7) = - 0.41577197673418943962;
    xtab( 8) = - 0.14003444424696773778;
    xtab( 9) =   0.14003444424696773778;
    xtab(10) =   0.41577197673418943962;
    xtab(11) =   0.67861108097560545347;
    xtab(12) =   0.92027786206637096497;
    xtab(13) =   1.1330068786005003695;
    xtab(14) =   1.3097818904452936698;
    xtab(15) =   1.4446014873666514608;
    xtab(16) =   1.5327507132362304779;

    weight( 1) = 0.0024194677567615628193;
    weight( 2) = 0.014115268156854008264;
    weight( 3) = 0.040437893946503669410;
    weight( 4) = 0.083026647573217742131;
    weight( 5) = 0.13834195526951273359;
    weight( 6) = 0.19741148870253455567;
    weight( 7) = 0.24763632094635522403;
    weight( 8) = 0.27661095764826050408;
    weight( 9) = 0.27661095764826050408;
    weight(10) = 0.24763632094635522403;
    weight(11) = 0.19741148870253455567;
    weight(12) = 0.13834195526951273359;
    weight(13) = 0.083026647573217742131;
    weight(14) = 0.040437893946503669410;
    weight(15) = 0.014115268156854008264;
    weight(16) = 0.0024194677567615628193;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'LEGENDRE_SET_COS - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    error ( 'LEGENDRE_SET_COS - Fatal error!' );

  end
