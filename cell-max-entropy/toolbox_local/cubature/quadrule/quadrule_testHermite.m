%% TEST085 tests HERMITE_COMPUTE against HERMITE_INTEGRAL.
%
%  Modified:
%
%    25 July 2007
%
%  Author:
%
%    John Burkardt
%
  order_max = 10;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST085\n' );
  fprintf ( 1, '  HERMITE_COMPUTE computes a Gauss-Hermite rule\n' );
  fprintf ( 1, '  which is appropriate for integrands of the form\n' );
  fprintf ( 1, '    f(x) * exp(-x**2) from -infinity to infinity.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  HERMITE_INTEGRAL determines the exact value of\n' );
  fprintf ( 1, '  this integal when f(x) = x^n.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, ...
    '         N     Order       Estimate       Exact            Error\n' );

  for n = 0 : 2 : 4

    exact = hermite_integral ( n );

    fprintf ( 1, '\n' );

    for order = 1 : order_max

      [ xtab, weight ] = hermite_compute ( order );

      if ( n == 0 )
        f_vec(1:order) = 1.0;
      else
        f_vec(1:order) = xtab(1:order).^n;
      end
      estimate = weight(1:order) * f_vec(1:order)';

      err = abs ( exact - estimate );

      fprintf ( 1, '  %8d  %8d  %14f  %14f  %14g\n', ...
        n, order, estimate, exact, err );


    end

  end
