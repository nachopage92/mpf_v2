function weight = nc_compute ( order, x_min, x_max, xtab )

%% NC_COMPUTE computes a Newton-Cotes quadrature rule.
%
%  Discussion:
%
%    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule 
%    estimates
%
%      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
%
%    using ORDER abscissas XTAB(I) and a weight vector WEIGHT(I):
%
%      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) ).
%
%    For the CLOSED rule, the equally spaced abscissas include A and B.
%    For the OPEN rule, the equally spaced abscissas do not include A and B.
%
%  Modified:
%
%    26 May 2007
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer ORDER, the order of the rule.
%
%    Input, real X_MIN, X_MAX, the endpoints of the interval..
%
%    Input, real XTAB(ORDER), the abscissas.
%
%    Output, real WEIGHT(ORDER), the weights.
%
  for i = 1 : order
%
%  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
%  and zero at the other nodes.
%
    diftab(1:order) = 0.0;
    diftab(i) = 1.0;

    for j = 2 : order
      for k = j : order
        diftab(order+j-k) = ( diftab(order+j-k-1) - diftab(order+j-k) ) ...
          / ( xtab(order+1-k) - xtab(order+j-k) );
      end
    end

    for j = 1 : order-1
      for k = 1 : order-j
        diftab(order-k) = diftab(order-k) ...
          - xtab(order-k-j+1) * diftab(order-k+1);
      end
    end
%
%  Evaluate the antiderivative of the polynomial at the left and
%  right endpoints.
%
    yvala = diftab(order) / order;
    for j = order-1 : -1 : 1
      yvala = yvala * x_min + diftab(j) / j;
    end
    yvala = yvala * x_min;

    yvalb = diftab(order) / order;
    for j = order-1 : -1 : 1
      yvalb = yvalb * x_max + diftab(j) / j;
    end
    yvalb = yvalb * x_max;

    weight(i) = yvalb - yvala;

  end
