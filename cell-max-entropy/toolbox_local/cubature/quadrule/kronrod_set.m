function [ xtab, weight ] = kronrod_set ( norder )

%% KRONROD_SET sets abscissas and weights for Gauss-Kronrod quadrature.
%
%  Discussion:
%
%    The integration interval is [ -1, 1 ].
%
%    The weight function is w(x) = 1.
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%    A Kronrod rule is used in conjunction with a lower order
%    Gauss rule, and provides an efficient error estimation.
%
%    The error may be estimated as the difference in the two integral
%    approximations.
%
%    The efficiency comes about because the Kronrod uses the abscissas
%    of the Gauss rule, thus saving on the number of function evaluations
%    necessary.  If the Kronrod rule were replaced by a Gauss rule of
%    the same order, a higher precision integral estimate would be
%    made, but the function would have to be evaluated at many more
%    points.
%
%    The Gauss Kronrod pair of rules involves an ( NORDER + 1 ) / 2
%    point Gauss-Legendre rule and an NORDER point Kronrod rule.
%    Thus, the 15 point Kronrod rule should be paired with the
%    Gauss-Legendre 7 point rule.
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
%    R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner,
%    QUADPACK, A Subroutine Package for Automatic Integration,
%    Springer Verlag, 1983.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule, which may be
%    15, 21, 31 or 41, corresponding to Gauss-Legendre rules of
%    order 7, 10, 15 or 20.
%
%    Output, real XTAB(NORDER), the abscissas of the rule, which
%    are symmetrically places in [-1,1].
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%    The weights are positive, symmetric, and should sum to 2.
%
  if ( norder == 15 )

    xtab(1) =  - 0.9914553711208126;
    xtab(2) =  - 0.9491079123427585;
    xtab(3) =  - 0.8648644233597691;
    xtab(4) =  - 0.7415311855993944;
    xtab(5) =  - 0.5860872354676911;
    xtab(6) =  - 0.4058451513773972;
    xtab(7) =  - 0.2077849550789850;
    xtab(8) =    0.0;
    xtab(9) =    0.2077849550789850;
    xtab(10) =   0.4058451513773972;
    xtab(11) =   0.5860872354676911;
    xtab(12) =   0.7415311855993944;
    xtab(13) =   0.8648644233597691;
    xtab(14) =   0.9491079123427585;
    xtab(15) =   0.9914553711208126;

    weight(1) =  0.2293532201052922E-01;
    weight(2) =  0.6309209262997855E-01;
    weight(3) =  0.1047900103222502;
    weight(4) =  0.1406532597155259;
    weight(5) =  0.1690047266392679;
    weight(6) =  0.1903505780647854;
    weight(7) =  0.2044329400752989;
    weight(8) =  0.2094821410847278;
    weight(9) =  0.2044329400752989;
    weight(10) = 0.1903505780647854;
    weight(11) = 0.1690047266392679;
    weight(12) = 0.1406532597155259;
    weight(13) = 0.1047900103222502;
    weight(14) = 0.6309209262997855E-01;
    weight(15) = 0.2293532201052922E-01;

  elseif ( norder == 21 )

    xtab(1) =  - 0.9956571630258081;
    xtab(2) =  - 0.9739065285171717;
    xtab(3) =  - 0.9301574913557082;
    xtab(4) =  - 0.8650633666889845;
    xtab(5) =  - 0.7808177265864169;
    xtab(6) =  - 0.6794095682990244;
    xtab(7) =  - 0.5627571346686047;
    xtab(8) =  - 0.4333953941292472;
    xtab(9) =  - 0.2943928627014602;
    xtab(10) = - 0.1488743389816312;
    xtab(11) =   0.0;
    xtab(12) =   0.1488743389816312;
    xtab(13) =   0.2943928627014602;
    xtab(14) =   0.4333953941292472;
    xtab(15) =   0.5627571346686047;
    xtab(16) =   0.6794095682990244;
    xtab(17) =   0.7808177265864169;
    xtab(18) =   0.8650633666889845;
    xtab(19) =   0.9301574913557082;
    xtab(20) =   0.9739065285171717;
    xtab(21) =   0.9956571630258081;

    weight(1) =  0.1169463886737187E-01;
    weight(2) =  0.3255816230796473E-01;
    weight(3) =  0.5475589657435200E-01;
    weight(4) =  0.7503967481091995E-01;
    weight(5) =  0.9312545458369761E-01;
    weight(6) =  0.1093871588022976;
    weight(7) =  0.1234919762620659;
    weight(8) =  0.1347092173114733;
    weight(9) =  0.1427759385770601;
    weight(10) = 0.1477391049013385;
    weight(11) = 0.1494455540029169;
    weight(12) = 0.1477391049013385;
    weight(13) = 0.1427759385770601;
    weight(14) = 0.1347092173114733;
    weight(15) = 0.1234919762620659;
    weight(16) = 0.1093871588022976;
    weight(17) = 0.9312545458369761E-01;
    weight(18) = 0.7503967481091995E-01;
    weight(19) = 0.5475589657435200E-01;
    weight(20) = 0.3255816230796473E-01;
    weight(21) = 0.1169463886737187E-01;

  elseif ( norder == 31 )

    xtab(1) =  - 0.9980022986933971;
    xtab(2) =  - 0.9879925180204854;
    xtab(3) =  - 0.9677390756791391;
    xtab(4) =  - 0.9372733924007059;
    xtab(5) =  - 0.8972645323440819;
    xtab(6) =  - 0.8482065834104272;
    xtab(7) =  - 0.7904185014424659;
    xtab(8) =  - 0.7244177313601700;
    xtab(9) =  - 0.6509967412974170;
    xtab(10) = - 0.5709721726085388;
    xtab(11) = - 0.4850818636402397;
    xtab(12) = - 0.3941513470775634;
    xtab(13) = - 0.2991800071531688;
    xtab(14) = - 0.2011940939974345;
    xtab(15) = - 0.1011420669187175;
    xtab(16) =   0.0;
    xtab(17) =   0.1011420669187175;
    xtab(18) =   0.2011940939974345;
    xtab(19) =   0.2991800071531688;
    xtab(20) =   0.3941513470775634;
    xtab(21) =   0.4850818636402397;
    xtab(22) =   0.5709721726085388;
    xtab(23) =   0.6509967412974170;
    xtab(24) =   0.7244177313601700;
    xtab(25) =   0.7904185014424659;
    xtab(26) =   0.8482065834104272;
    xtab(27) =   0.8972645323440819;
    xtab(28) =   0.9372733924007059;
    xtab(29) =   0.9677390756791391;
    xtab(30) =   0.9879925180204854;
    xtab(31) =   0.9980022986933971;

    weight(1) =  0.5377479872923349E-02;
    weight(2) =  0.1500794732931612E-01;
    weight(3) =  0.2546084732671532E-01;
    weight(4) =  0.3534636079137585E-01;
    weight(5) =  0.4458975132476488E-01;
    weight(6) =  0.5348152469092809E-01;
    weight(7) =  0.6200956780067064E-01;
    weight(8) =  0.6985412131872826E-01;
    weight(9) =  0.7684968075772038E-01;
    weight(10) = 0.8308050282313302E-01;
    weight(11) = 0.8856444305621177E-01;
    weight(12) = 0.9312659817082532E-01;
    weight(13) = 0.9664272698362368E-01;
    weight(14) = 0.9917359872179196E-01;
    weight(15) = 0.1007698455238756;
    weight(16) = 0.1013300070147915;
    weight(17) = 0.1007698455238756;
    weight(18) = 0.9917359872179196E-01;
    weight(19) = 0.9664272698362368E-01;
    weight(20) = 0.9312659817082532E-01;
    weight(21) = 0.8856444305621177E-01;
    weight(22) = 0.8308050282313302E-01;
    weight(23) = 0.7684968075772038E-01;
    weight(24) = 0.6985412131872826E-01;
    weight(25) = 0.6200956780067064E-01;
    weight(26) = 0.5348152469092809E-01;
    weight(27) = 0.4458975132476488E-01;
    weight(28) = 0.3534636079137585E-01;
    weight(29) = 0.2546084732671532E-01;
    weight(30) = 0.1500794732931612E-01;
    weight(31) = 0.5377479872923349E-02;

  elseif ( norder == 41 )

    xtab(1) =  - 0.9988590315882777;
    xtab(2) =  - 0.9931285991850949;
    xtab(3) =  - 0.9815078774502503;
    xtab(4) =  - 0.9639719272779138;
    xtab(5) =  - 0.9408226338317548;
    xtab(6) =  - 0.9122344282513259;
    xtab(7) =  - 0.8782768112522820;
    xtab(8) =  - 0.8391169718222188;
    xtab(9) =  - 0.7950414288375512;
    xtab(10) = - 0.7463319064601508;
    xtab(11) = - 0.6932376563347514;
    xtab(12) = - 0.6360536807265150;
    xtab(13) = - 0.5751404468197103;
    xtab(14) = - 0.5108670019508271;
    xtab(15) = - 0.4435931752387251;
    xtab(16) = - 0.3737060887154196;
    xtab(17) = - 0.3016278681149130;
    xtab(18) = - 0.2277858511416451;
    xtab(19) = - 0.1526054652409227;
    xtab(20) = - 0.7652652113349733E-01;
    xtab(21) =   0.0;
    xtab(22) =   0.7652652113349733E-01;
    xtab(23) =   0.1526054652409227;
    xtab(24) =   0.2277858511416451;
    xtab(25) =   0.3016278681149130;
    xtab(26) =   0.3737060887154196;
    xtab(27) =   0.4435931752387251;
    xtab(28) =   0.5108670019508271;
    xtab(29) =   0.5751404468197103;
    xtab(30) =   0.6360536807265150;
    xtab(31) =   0.6932376563347514;
    xtab(32) =   0.7463319064601508;
    xtab(33) =   0.7950414288375512;
    xtab(34) =   0.8391169718222188;
    xtab(35) =   0.8782768112522820;
    xtab(36) =   0.9122344282513259;
    xtab(37) =   0.9408226338317548;
    xtab(38) =   0.9639719272779138;
    xtab(39) =   0.9815078774502503;
    xtab(40) =   0.9931285991850949;
    xtab(41) =   0.9988590315882777;

    weight(1) =  0.3073583718520532E-02;
    weight(2) =  0.8600269855642942E-02;
    weight(3) =  0.1462616925697125E-01;
    weight(4) =  0.2038837346126652E-01;
    weight(5) =  0.2588213360495116E-01;
    weight(6) =  0.3128730677703280E-01;
    weight(7) =  0.3660016975820080E-01;
    weight(8) =  0.4166887332797369E-01;
    weight(9) =  0.4643482186749767E-01;
    weight(10) = 0.5094457392372869E-01;
    weight(11) = 0.5519510534828599E-01;
    weight(12) = 0.5911140088063957E-01;
    weight(13) = 0.6265323755478117E-01;
    weight(14) = 0.6583459713361842E-01;
    weight(15) = 0.6864867292852162E-01;
    weight(16) = 0.7105442355344407E-01;
    weight(17) = 0.7303069033278667E-01;
    weight(18) = 0.7458287540049919E-01;
    weight(19) = 0.7570449768455667E-01;
    weight(20) = 0.7637786767208074E-01;
    weight(21) = 0.7660071191799966E-01;
    weight(22) = 0.7637786767208074E-01;
    weight(23) = 0.7570449768455667E-01;
    weight(24) = 0.7458287540049919E-01;
    weight(25) = 0.7303069033278667E-01;
    weight(26) = 0.7105442355344407E-01;
    weight(27) = 0.6864867292852162E-01;
    weight(28) = 0.6583459713361842E-01;
    weight(29) = 0.6265323755478117E-01;
    weight(30) = 0.5911140088063957E-01;
    weight(31) = 0.5519510534828599E-01;
    weight(32) = 0.5094457392372869E-01;
    weight(33) = 0.4643482186749767E-01;
    weight(34) = 0.4166887332797369E-01;
    weight(35) = 0.3660016975820080E-01;
    weight(36) = 0.3128730677703280E-01;
    weight(37) = 0.2588213360495116E-01;
    weight(38) = 0.2038837346126652E-01;
    weight(39) = 0.1462616925697125E-01;
    weight(40) = 0.8600269855642942E-02;
    weight(41) = 0.3073583718520532E-02;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'KRONROD_SET - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of NORDER = %d\n', norder );
    fprintf ( 1, '  Legal values are 15, 21, 31 or 41.\n' );
    error ( 'KRONROD_SET - Fatal error!' );

  end 

