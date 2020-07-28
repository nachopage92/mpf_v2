function [ xtab, weight ] = bdfp_set ( norder )

%% BDFP_SET sets weights for backward differentiation predictor quadrature.
%
%  Discussion:
%
%    A backward differentiation predictor formula is defined for a set
%    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
%    that the values of the function to be integrated are known at the
%    abscissas, the formula is written in terms of the function value at
%    X(2), and the backward differences at X(2) that approximate the
%    derivatives there.  A backward differentiation predictor formula
%    is equivalent to an Adams-Bashforth formula of the same order.
%
%    The integration interval is [ 0, 1 ].
%
%    The weight function is w(x) = 1.0;
%
%    The integral to approximate:
%
%      Integral ( 0 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * BD**(I-1) F ( 0 ),
%
%    Here, "BD**(I-1) F ( 0 )" denotes the (I-1)st backward difference
%    of F at X = 0, using a spacing of 1.  In particular,
%
%    BD**0 F(0) = F(0)
%    BD**1 F(0) = F(0) - F(-1)
%    BD**2 F(0) = F(0) - 2 * F(-1) + F(-2 )
%
%    The relationship between a backward difference predictor and the
%    corresponding Adams-Bashforth formula may be illustrated for the
%    BDF predictor of order 3:
%
%      BD**0 F(0) + 0.5 * BD**1 F(0) + 5/12 * BD**2 F(0)
%      =            F(0)
%        + 1/2  * ( F(0) -         F(1) )
%        + 5/12 * ( F(0) - 2     * F(-1) +      F(-2) )
%      =  23/12 *   F(0) - 16/12 * F(-1) + 5/12 F(-2)
%
%    which is the Adams-Bashforth formula of order 3.
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
%    Output, real WEIGHT(NORDER), the weight of the rule.
%
  order_max = 19;

  w(1) =                       1.0;
  w(2) =                       1.0 /                2.0;
  w(3) =                       5.0 /               12.0;
  w(4) =                       3.0 /                8.0;
  w(5) =                     251.0 /              720.0;
  w(6) =                      95.0 /              288.0;
  w(7) =                   19087.0 /            60480.0;
  w(8) =                    5257.0 /            17280.0;
  w(9) =                 1070017.0 /          3628800.0;
  w(10) =                  25713.0 /            89600.0;
  w(11) =               26842253.0 /         95800320.0;
  w(12) =                4777223.0 /         17418240.0;
  w(13) =           703604254357.0 /    2615348736000.0;
  w(14) =           106364763817.0 /     402361344000.0;
  w(15) =          1166309819657.0 /    4483454976000.0;
  w(16) =               25221445.0 /         98402304.0;
  w(17) =       8092989203533249.0 /    3201186852864.0;
  w(18) =         85455477715379.0 /      34237292544.0;
  w(19) =   12600467236042756559.0 / 5109094217170944.0;

  if ( order_max < norder )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'BDFP_SET - Fatal error!\n' );
    fprintf ( 1, '  Input value NORDER = %d exceeds ORDER_MAX = %d\n', ...
      norder, order_max );
    error ( 'BDFP_SET - Fatal error!' );
  end

  weight(1:norder) = w(1:norder);

  for i = 1 : norder
    xtab(i) = 1 - i;
  end
