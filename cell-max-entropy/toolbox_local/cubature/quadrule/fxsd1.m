function value = fxsd1 ( x )

%% FXSD1 evaluates the function x / sqrt ( 1.1 - x**2 ).
%
%  Modified:
%
%    12 October 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the argument of the function.
%
%    Output, real FXSD1, the value of the function.
%
  value = x / sqrt ( 1.1 - x * x );
