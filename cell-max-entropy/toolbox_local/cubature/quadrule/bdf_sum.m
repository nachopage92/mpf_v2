function result = bdf_sum ( func, norder, xtab, weight )

%% BDF_SUM carries out an explicit backward difference quadrature rule for [0,1].
%
%  Discussion:
%
%    The integral to approximate is
%
%      Integral ( 0 <= X <= 1 ) F(X) dX
%
%    The quadrature formula is:
%
%      RESULT = Sum ( 1 <= I <= NORDER ) WEIGHT(I) * BDF**(I-1) FUNC ( 0 )
%
%    The integral from 0 to 1 is approximated using data at X = 0,
%    -1, -2, ..., -NORDER+1.  This is a form of extrapolation, and
%    the approximation can become poor as NORDER increases.
%
%  Modified:
%
%    13 October 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, external FUNC, the name of the function which evaluates
%    the integrand.  The function must have the form
%      function value = func ( x ).
%
%    Input, integer NORDER, the order of the rule.
%
%    Input, real XTAB(NORDER), the abscissas of the rule.
%
%    Input, real WEIGHT(NORDER), the weights of the rule.
%
%    Output, real RESULT, the approximate value of the integral.
%
  for i = 1 : norder
    diftab(i) = func ( xtab(i) );
  end

  for i = 2 : norder
    for j = i : norder
      diftab(norder+i-j) = ( diftab(norder+i-j-1) - diftab(norder+i-j) );
    end
  end

  result = weight(1:norder) * diftab(1:norder)';
