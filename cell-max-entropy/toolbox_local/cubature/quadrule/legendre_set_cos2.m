function [ xtab, weight ] = legendre_set_cos2 ( norder )

%% LEGENDRE_SET_COS2 sets a Gauss-Legendre rule for COS(X) * F(X) on [0,PI/2].
%
%  Discussion:
%
%    The integration interval is [ 0, PI/2 ].
%
%    The weight function is w(x) = COS(X).
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= PI/2 ) COS(X) * F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%  Discussion:
%
%    The same rule can be used to approximate
%
%      Integral ( 0 <= X <= PI/2 ) SIN(X) * F(X) dX
%
%    as
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( PI/2 - XTAB(I) )
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
%    Wiley, 1993, QA299.3E93, page 311.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule.
%    NORDER must be between 2, 4, 8 or 16.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%
  if ( norder == 2 )

    xtab(1) = 0.26587388056307823382;
    xtab(2) = 1.0351526093171315182;

    weight(1) = 0.60362553280827113087;
    weight(2) = 0.39637446719172886913;

  elseif ( norder == 4 )

    xtab(1) = 0.095669389196858636773;
    xtab(2) = 0.45240902327067096554;
    xtab(3) = 0.93185057672024082424;
    xtab(4) = 1.3564439599666466230;

    weight( 1) = 0.23783071419515504517;
    weight( 2) = 0.40265695523581253512;
    weight( 3) = 0.28681737948564715225;
    weight( 4) = 0.072694951083385267446;

  elseif ( norder == 8 )

    xtab(1) = 0.029023729768913933432;
    xtab(2) = 0.14828524404581819442;
    xtab(3) = 0.34531111151664787488;
    xtab(4) = 0.59447696797658360178;
    xtab(5) = 0.86538380686123504827;
    xtab(6) = 1.1263076093187456632;
    xtab(7) = 1.3470150460281258016;
    xtab(8) = 1.5015603622059195568;

    weight( 1) = 0.073908998095117384985;
    weight( 2) = 0.16002993702338006099;
    weight( 3) = 0.21444434341803549108;
    weight( 4) = 0.21979581268851903339;
    weight( 5) = 0.17581164478209568886;
    weight( 6) = 0.10560448025308322171;
    weight( 7) = 0.042485497299217201089;
    weight( 8) = 0.0079192864405519178899;

  elseif ( norder == 16 )

    xtab( 1) = 0.0080145034906295973494;
    xtab( 2) = 0.041893031354246254797;
    xtab( 3) = 0.10149954486757579459;
    xtab( 4) = 0.18463185923836617507;
    xtab( 5) = 0.28826388487760574589;
    xtab( 6) = 0.40870579076464794191;
    xtab( 7) = 0.54176054986913847463;
    xtab( 8) = 0.68287636658719416893;
    xtab( 9) = 0.82729287620416833520;
    xtab(10) = 0.97018212594829367065;
    xtab(11) = 1.1067865150286247873;
    xtab(12) = 1.2325555697227748824;
    xtab(13) = 1.3432821921580721861;
    xtab(14) = 1.4352370549295032923;
    xtab(15) = 1.5052970876794669248;
    xtab(16) = 1.5510586944086135769;

    weight( 1) = 0.020528714977215248902;
    weight( 2) = 0.046990919853597958123;
    weight( 3) = 0.071441021312218541698;
    weight( 4) = 0.092350338329243052271;
    weight( 5) = 0.10804928026816236935;
    weight( 6) = 0.11698241243306261791;
    weight( 7) = 0.11812395361762037649;
    weight( 8) = 0.11137584940420091049;
    weight( 9) = 0.097778236145946543110;
    weight(10) = 0.079418758985944482077;
    weight(11) = 0.059039620053768691402;
    weight(12) = 0.039458876783728165671;
    weight(13) = 0.022987785677206847531;
    weight(14) = 0.011010405600421536861;
    weight(15) = 0.0038123928030499915653;
    weight(16) = 0.00065143375461266656171;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'LEGENDRE_SET_COS2 - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    error ( 'LEGENDRE_SET_COS2 - Fatal error!' );

  end

