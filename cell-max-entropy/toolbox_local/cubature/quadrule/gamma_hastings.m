function value = gamma_hastings ( x )

%% GAMMA_HASTINGS computes Hastings's function for approximation to the Gamma function.
%
%  Modified:
%
%    12 October 2005
%
%  Reference:
%
%    Arthur Stroud and Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966.
%
%  Parameters:
%
%    Input, real X, the argument at which the function is to be evaluated.
%
%    Output, real VALUE, the Hastings function used in approximation of
%    the Gamma function.
%
  value = ((((((( ...
          0.035868343   * x ...
        - 0.193527818 ) * x ...
        + 0.482199394 ) * x ...
        - 0.756704078 ) * x ...
        + 0.918206857 ) * x ...
        - 0.897056937 ) * x ...
        + 0.988205891 ) * x ...
        - 0.577191652 ) * x + 1.0;
