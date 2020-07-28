function [ xtab, weight ] = bdfc_set ( norder )

%% BDFC_SET sets weights for backward differentiation corrector quadrature.
%
%  Definition:
%
%    A backward differentiation corrector formula is defined for a set
%    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
%    that the values of the function to be integrated are known at the
%    abscissas, the formula is written in terms of the function value at
%    X(1), and the backward differences at X(1) that approximate the
%    derivatives there.
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function w(x) = 1.0;
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * BD**(I-1) F ( 1 ).
%
%    Here, "BD**(I-1) F ( 1 )" denotes the (I-1)st backward difference
%    of F at X = 1, using a spacing of 1.  In particular,
%
%    BD**0 F(1) = F(1)
%    BD**1 F(1) = F(1) - F(0)
%    BD**2 F(1) = F(1) - 2 * F(0) + F(-1 )
%
%
%    The relationship between a backward difference corrector and the
%    corresponding Adams-Moulton formula may be illustrated for the
%    BDF corrector of order 4:
%
%      BD**0 F(1) - 1/2 * BD**1 F(1) - 1/12 * BD**2 F(1) - 1/24 * BDF**3 F(1)
%      =            F(1)
%        -  1/2 * ( F(1) -         F(0) )
%        - 1/12 * ( F(1) - 2     * F(0) +        F(-1) )
%        - 1/24 * ( F(1) - 3     * F(0) + 3    * F(-1) -        F(-2) )
%      =   9/24 *   F(1) + 19/24 * F(0) - 5/24 * F(-1) + 1/24 * F(-2)
%
%    which is the Adams-Moulton formula of order 4.
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
%    Simeon Fatunla,
%    Numerical Methods for Initial Value Problems in Ordinary Differential
%      Equations,
%    Academic Press, 1988.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the rule, which can be
%    any value from 1 to 19.
%
%    Output, real XTAB(NORDER), the abscissas of the rule.
%
%    Output, real WEIGHT(NORDER), the weights of the rule.
%
  order_max = 19;

  w(1) =                 1.0;
  w(2) =               - 1.0 /               2.0;
  w(3) =               - 1.0 /              12.0;
  w(4) =               - 1.0 /              24.0;
  w(5) =              - 19.0 /             720.0;
  w(6) =               - 3.0 /             160.0;
  w(7) =             - 863.0 /           60480.0;
  w(8) =             - 275.0 /           24792.0;
  w(9) =           - 33953.0 /         3628800.0;
  w(10) =           - 8183.0 /         1036800.0;
  w(11) =        - 3250433.0 /       479001600.0;
  w(12) =           - 4671.0 /          788480.0;
  w(13) =    - 13695779093.0 /   2615348736000.0;
  w(14) =     - 2224234463.0 /    475517952000.0;
  w(15) =   - 132282840127.0 /  31384184832000.0;
  w(16) =     - 2639651053.0 /    689762304000.0;
  w(17) =  111956703448001.0 /   3201186852864.0;
  w(18) =         50188465.0 /     15613165568.0;
  w(19) = 2334028946344463.0 / 786014494949376.0;

  if ( order_max < norder )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'BDFC_SET - Fatal error!\n' );
    fprintf ( 1, '  Input order %d exceeds maximum %d\n', norder, order_max );
    error ( 'BDFC_SET - Fatal error!' );
  end

  for i = 1 : norder
    weight(i) = w(i);
  end

  for i = 1 : norder
    xtab(i) = ( 2 - i );
  end
