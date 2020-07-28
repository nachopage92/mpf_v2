function [ x, w ] = clenshaw_curtis_compute ( n )

%% CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
%
%  Discussion:
%
%    This method uses a direct approach.  The paper by Waldvogel
%    exhibits a more efficient approach using Fourier transforms.
%
%    The integration interval is [ -1, 1 ].
%
%    The integral to approximate:
%
%      Integral ( -1 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
%
%    The abscissas for the rule of order ORDER can be regarded 
%    as the cosines of equally spaced angles between 180 and 0 degrees:
%
%      X(I) = cos ( ( I - 1 ) * PI / ( ORDER - 1 ) )
%
%    except for the basic case ORDER = 1, when
%
%      X(1) = 0.
%
%    A Clenshaw-Curtis rule that uses ORDER points will integrate
%    exactly all polynomials of degrees 0 through ORDER-1.  If ORDER
%    is odd, then by symmetry the polynomial of degree ORDER will
%    also be integrated exactly.
%
%    If the value of ORDER is increased in a sensible way, then
%    the new set of abscissas will include the old ones.  One such
%    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...
%
%  Modified:
%
%    05 March 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Charles Clenshaw, Alan Curtis,
%    A Method for Numerical Integration on an Automatic Computer,
%    Numerische Mathematik,
%    Volume 2, Number 1, December 1960, pages 197-205.
%
%    Philip Davis, Philip Rabinowitz,
%    Methods of Numerical Integration,
%    Second Edition,
%    Dover, 2007,
%    ISBN: 0486453391,
%    LC: QA299.3.D28.
%
%    Joerg Waldvogel,
%    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
%    BIT Numerical Mathematics,
%    Volume 43, Number 1, 2003, pages 1-18.
%
%  Parameters:
%
%    Input, integer N, the order of the rule.
%
%    Output, real X(N), W(N), the abscissas and weights of the rule.
%
  if ( n < 1 )

    x = [];
    w = [];
    return

  elseif ( n == 1 )

    x(1) = 0.0;
    w(1) = 2.0;
    return

  end

  theta(1:n) = ( n-1 : -1 : 0 ) * pi / ( n - 1 );
  x(1:n) = cos ( theta(1:n) );

  for i = 1 : n

    w(i) = 1;

    for j = 1 : floor ( ( n - 1 ) / 2 )

      if ( 2 * j == ( n - 1 ) )
        b = 1;
      else
        b = 2;
      end

      w(i) = w(i) - b * cos ( 2 * j * theta(i) ) / ( 4 * j * j - 1 );

    end
  end

  w(1)     =     w(1)     / ( n - 1 );
  w(2:n-1) = 2 * w(2:n-1) / ( n - 1 );
  w(n)     =     w(n)     / ( n - 1 );
