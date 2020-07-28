function [ xtab, weight ] = bashforth_set ( norder )

%% BASHFORTH_SET sets abscissas and weights for Adams-Bashforth quadrature.
%
%  Discussion:
%
%    Adams-Bashforth quadrature formulas are normally used in solving
%    ordinary differential equations, and are not really suitable for
%    general quadrature computations.  However, an Adams-Bashforth formula
%    is equivalent to approximating the integral of F(Y(X)) between X(M)
%    and X(M+1), using an explicit formula that relies only on known values
%    of F(Y(X)) at X(M-N+1) through X(M).  For this reason, the formulas
%    have been included here.
%
%    Suppose the unknown function is denoted by Y(X), with derivative
%    F(Y(X)), and that approximate values of the function are known at a
%    series of X values, which we write as X(1), X(2), ..., X(M).  We write
%    the value Y(X(1)) as Y(1) and so on.
%
%    Then the solution of the ODE Y'=F(X,Y) at the next point X(M+1) is
%    computed by:
%
%      Y(M+1) = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
%             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+1-I)) approximately.
%
%    In the documentation that follows, we replace F(Y(X)) by F(X).
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function is w(x) = 1.0;
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) F(X) dX.
%
%    The quadrature formula:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( 1 - I ),
%
%    The Adams-Bashforth formulas require equally spaced data.
%
%    Here is how the formula is applied in the case with non-unit spacing:
%
%      Integral ( A <= X <= A+H ) F(X) dX =
%      H * Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( A - (I-1)*H ),
%      approximately.
%
%    The reference lists the second coefficient of the order 8 Adams-Bashforth
%    formula as
%      weight(2) =  -1162169.0 / 120960.0
%    but this should be
%      weight(2) =  -1152169.0 / 120960.0
%
%  Modified:
%
%    28 April 2006
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
%    Jean Lapidus and John Seinfeld,
%    Numerical Solution of Ordinary Differential Equations,
%    Academic Press, 1971.
%
%    Daniel Zwillinger, editor,
%    Standard Mathematical Tables and Formulae,
%    30th Edition,
%    CRC Press, 1996.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule.  NORDER should be
%    between 1 and 10, or 12, 14, 16, 18 or 20.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%    WEIGHT(1) is the weight at X = 0, WEIGHT(2) the weight at X = -1,
%    and so on.  The weights are rational, and should sum to 1.  Some
%    weights may be negative.
%
  if ( norder == 1 )

    weight(1) = 1.0;

  elseif ( norder == 2 )

    d = 2.0;

    weight(1) =   3.0 / d;
    weight(2) = - 1.0 / d;

  elseif ( norder == 3 )

    d = 12.0;

    weight(1) =   23.0 / d;
    weight(2) = - 16.0 / d;
    weight(3) =    5.0 / d;

  elseif ( norder == 4 )

    d = 24.0;

    weight(1) =   55.0 / d;
    weight(2) = - 59.0 / d;
    weight(3) =   37.0 / d;
    weight(4) =  - 9.0 / d;

  elseif ( norder == 5 )

    d = 720.0;

    weight(1) =   1901.0 / d;
    weight(2) = - 2774.0 / d;
    weight(3) =   2616.0 / d;
    weight(4) = - 1274.0 / d;
    weight(5) =    251.0 / d;

  elseif ( norder == 6 )

    d = 1440.0;

    weight(1) =   4277.0 / d;
    weight(2) = - 7923.0 / d;
    weight(3) =   9982.0 / d;
    weight(4) = - 7298.0 / d;
    weight(5) =   2877.0 / d;
    weight(6) =  - 475.0 / d;

  elseif ( norder == 7 )

    d = 60480.0;

    weight(1) =    198721.0 / d;
    weight(2) =  - 447288.0 / d;
    weight(3) =    705549.0 / d;
    weight(4) =  - 688256.0 / d;
    weight(5) =    407139.0 / d;
    weight(6) =  - 134472.0 / d;
    weight(7) =     19087.0 / d;

  elseif ( norder == 8 )

    d = 120960.0;

    weight(1) =     434241.0 / d;
    weight(2) =  - 1152169.0 / d;
    weight(3) =    2183877.0 / d;
    weight(4) =  - 2664477.0 / d;
    weight(5) =    2102243.0 / d;
    weight(6) =  - 1041723.0 / d;
    weight(7) =     295767.0 / d;
    weight(8) =    - 36799.0 / d;

 elseif ( order == 9 )
  
    d = 3628800.0;

    weight(1) =   14097247.0 / d;
    weight(2) =  -43125206.0 / d;
    weight(3) =   95476786.0 / d;
    weight(4) = -139855262.0 / d;
    weight(5) =  137968480.0 / d;
    weight(6) =  -91172642.0 / d;
    weight(7) =   38833486.0 / d;
    weight(8) =   -9664106.0 / d;
    weight(9) =    1070017.0 / d;
  
  elseif ( order == 10 )
  
    d = 7257600.0;

    weight( 1) =   30277247.0 / d;
    weight( 2) = -104995189.0 / d;
    weight( 3) =  265932680.0 / d;
    weight( 4) = -454661776.0 / d;
    weight( 5) =  538363838.0 / d;
    weight( 6) = -444772162.0 / d;
    weight( 7) =  252618224.0 / d;
    weight( 8) =  -94307320.0 / d;
    weight( 9) =   20884811.0 / d;
    weight(10) =   -2082753.0 / d;
  
  elseif ( order == 12 )
  
    d = 958003200.0;

    weight( 1) =    4527766399.0 / d;
    weight( 2) =  -19433810163.0 / d;
    weight( 3) =   61633227185.0 / d;
    weight( 4) = -135579356757.0 / d;
    weight( 5) =  214139355366.0 / d;
    weight( 6) = -247741639374.0 / d;
    weight( 7) =  211103573298.0 / d;
    weight( 8) = -131365867290.0 / d;
    weight( 9) =   58189107627.0 / d;
    weight(10) =  -17410248271.0 / d;
    weight(11) =    3158642445.0 / d;
    weight(12) =    -262747265.0 / d;
  
  elseif ( order == 14 )
  
    d = 5230697472000.0;

    weight( 1) =    27511554976875.0 / d;
    weight( 2) =  -140970750679621.0 / d;
    weight( 3) =   537247052515662.0 / d;
    weight( 4) = -1445313351681906.0 / d;
    weight( 5) =  2854429571790805.0 / d;
    weight( 6) = -4246767353305755.0 / d;
    weight( 7) =  4825671323488452.0 / d;
    weight( 8) = -4204551925534524.0 / d;
    weight( 9) =  2793869602879077.0 / d;
    weight(10) = -1393306307155755.0 / d;
    weight(11) =   505586141196430.0 / d;
    weight(12) =  -126174972681906.0 / d;
    weight(13) =    19382853593787.0 / d;
    weight(14) =    -1382741929621.0 / d;
  
  elseif ( order == 16 )
  
    d = 62768369664000.0;

    weight( 1) =     362555126427073.0 / d;
    weight( 2) =   -2161567671248849.0 / d;
    weight( 3) =    9622096909515337.0 / d;
    weight( 4) =  -30607373860520569.0 / d;
    weight( 5) =   72558117072259733.0 / d;
    weight( 6) = -131963191940828581.0 / d;
    weight( 7) =  187463140112902893.0 / d;
    weight( 8) = -210020588912321949.0 / d;
    weight( 9) =  186087544263596643.0 / d;
    weight(10) = -129930094104237331.0 / d;
    weight(11) =   70724351582843483.0 / d;
    weight(12) =  -29417910911251819.0 / d;
    weight(13) =    9038571752734087.0 / d;
    weight(14) =   -1934443196892599.0 / d;
    weight(15) =     257650275915823.0 / d;
    weight(16) =     -16088129229375.0 / d;
  
  elseif ( order == 18 )
  
    d = 64023737057280000.0;

    weight( 1) =     401972381695456831.0 / d;
    weight( 2) =   -2735437642844079789.0 / d;
    weight( 3) =   13930159965811142228.0 / d;
    weight( 4) =  -51150187791975812900.0 / d;
    weight( 5) =  141500575026572531760.0 / d;
    weight( 6) = -304188128232928718008.0 / d;
    weight( 7) =  518600355541383671092.0 / d;
    weight( 8) = -710171024091234303204.0 / d;
    weight( 9) =  786600875277595877750.0 / d;
    weight(10) = -706174326992944287370.0 / d;
    weight(11) =  512538584122114046748.0 / d;
    weight(12) = -298477260353977522892.0 / d;
    weight(13) =  137563142659866897224.0 / d;
    weight(14) =  -49070094880794267600.0 / d;
    weight(15) =   13071639236569712860.0 / d;
    weight(16) =   -2448689255584545196.0 / d;
    weight(17) =     287848942064256339.0 / d;
    weight(18) =     -15980174332775873.0 / d;
  
  elseif ( order == 20 )
  
    d = 102181884343418880000.0;

    weight( 1) =      691668239157222107697.0 / d;
    weight( 2) =    -5292843584961252933125.0 / d;
    weight( 3) =    30349492858024727686755.0 / d;
    weight( 4) =  -126346544855927856134295.0 / d;
    weight( 5) =   399537307669842150996468.0 / d;
    weight( 6) =  -991168450545135070835076.0 / d;
    weight( 7) =  1971629028083798845750380.0 / d;
    weight( 8) = -3191065388846318679544380.0 / d;
    weight( 9) =  4241614331208149947151790.0 / d;
    weight(10) = -4654326468801478894406214.0 / d;
    weight(11) =  4222756879776354065593786.0 / d;
    weight(12) = -3161821089800186539248210.0 / d;
    weight(13) =  1943018818982002395655620.0 / d;
    weight(14) =  -970350191086531368649620.0 / d;
    weight(15) =   387739787034699092364924.0 / d;
    weight(16) =  -121059601023985433003532.0 / d;
    weight(17) =    28462032496476316665705.0 / d;
    weight(18) =    -4740335757093710713245.0 / d;
    weight(19) =      498669220956647866875.0 / d;
    weight(20) =      -24919383499187492303.0 / d;
  
  else

    fprintf ( 1, '\n' );
    fprintf ( 1,'BASHFORTH_SET - Fatal error!n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    fprintf ( 1, '  Legal values are 1 through 10,\n' );
    fprintf ( 1, '  or 12, 14, 16, 18 or 20.\n' );
    error ( 'BASHFORTH_SET - Fatal error!' );

  end

  for i = 1 : norder
    xtab(i) = 1 - i;
  end
