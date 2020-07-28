function [f, df, t, dt] = trimming(v1,v2,x_s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% =======
%
% Inputs:
% v1, v2: vertices of an edge in ccw (Counter Clock Wise) 
% x_s : [ px(1) py(1); ... ; px(m) py(m) ], the m sample points for which we want 
%      to compute the levet set function and its derivatives
%
% Outputs:
% f  : line distance function = [f_1; ...; f_m] (m x 1 array)
% df : gradient of the line distance function = [df_1; ...; df_m] (m x 2 array)
% t  : trimming function, distance from x to a circunference
% dt : gradient of the trimming function
%
% N. Sukumar: UC Davis, M. Arroyo UPC-BarcelonaTech
% D. Millan UPC-BarcelonaTech, February 2014 

m    = size(x_s,1);

f    = zeros(m,1);
df   = zeros(m,2);

t    = zeros(m,1);
dt   = zeros(m,2);


s  = v2 - v1;
d  = norm(s);
id = 1/d;

xc = (v1+v2)*0.5;

for j = 1:m
  x = x_s(j,:);

  s1 = x - v1;
 
  %distance from x to a line that passes through v1 and v2
  f(j)    = (s1(1)*s(2) - s1(2)*s(1))*id;
  df(j,:) = [s(2), -s(1)]*id;
  
  %trimming function, distance from x to a circunference
  t(j)    = ( (0.5*d)^2 - sum((x-xc).^2) )*id;
  dt(j,:) = -2*(x-xc)*id;
end

return
end
