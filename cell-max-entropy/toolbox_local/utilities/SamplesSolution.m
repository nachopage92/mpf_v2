function s_sol = SamplesSolution(u_num,s_near,p_s)
% Numerical solution on the sample points is computed
%
% u_num  : numerical solution in the nodes (u_num=K\f)
% s_nears: nearest neighbor list
% p_s    : shape functions for each sample point 

dim   = size(u_num,2);
n_s   = length(s_near);
s_sol = zeros(n_s,dim);

for k=1:n_s
  k_near = s_near{k};
  p_k    = p_s{k};
  n_k    = length(k_near);

  u_sk   = 0;
  for ia=1:n_k
    i    = k_near(ia);
    u_sk = u_sk + p_k(ia)*u_num(i,:);
  end
  s_sol(k,:) = u_sk;
end
