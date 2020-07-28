function u_n = BndFitting_lme(x_n, options)
% function u_n = BndFitting_lme(x_n, options)
%
% INPUT:
%     x_n: points in 1D
%
% OUTPUT:
%     u_n: value at the nodes/control points

%% =======================================================================================
%  Input parameters: LME options

optLME = options;

optLME.dim   = 1;
optLME.verb  = 0; % 0:off    1:on
optLME.knn   = 0;
optLME.grad  = 0;           % Computation of the Gradient 0:OFF 1:ON
optLME.hess  = 0;

nPts = length(x_n);
L0   = min(x_n);
L1   = max(x_n);

% -------------------------------------------------------------------
x_s = linspace(L0,L1,5*nPts)';
sPts  = size(x_s,1);  %total number of sample points



%% =======================================================================================
%  LME: basis functions computation

% Nodes thermalization, nodes and samples adjacency structures are computed
[beta_n,range_n] = NodalThermalization(x_n, optLME);
optLME.beta_n = beta_n;
optLME.range_n= range_n;

% The basis functions are computed
% adjacency structure with the nearest neighbors nodes to each sample point
% nodal shape parameter
s_near = SamplesAdjacency(x_n,x_s,range_n);

% Local-max entropy basis functions computation
optLME.s_near = s_near;

outLME  = wrapper_lme(x_n,x_s,optLME);

p_samp  = outLME.p_samp;


%% =======================================================================================
% Least squares system assembly 
u_s  = x_s.*(1-x_s);

nn=0;
for k=1:sPts
  nn = max(nn, length(s_near{k}));
end
nn = min(nn,nPts);

M = spalloc(nPts,nPts,nn*nPts);
f = zeros(nPts,1);

for k=1:sPts
  k_near = s_near{k};
  p_k    = p_samp{k};
  
  M(k_near,k_near) = M(k_near,k_near) + (p_k*p_k');
  f(k_near) = f(k_near) + p_k*u_s(k);
end

u_n = M\f;

u_sh = zeros(sPts,1);
for k=1:sPts
  k_near = s_near{k};
  p_k    = p_samp{k};
  
  u_sh(k) = p_k'*u_n(k_near);
end

figure(10)
plot(x_s(:,1),u_s,'bx-','linewidth',2)
hold on
plot(x_s(:,1),u_sh,'k+-','linewidth',2)
plot(x_n,u_n,'ro','markersize',12)
hold off
fprintf(1,'\t|u - u_h| in L2 is %8.2e\n',norm(u_s-u_sh));
pause(0.1)

