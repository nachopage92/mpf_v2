function R = Unstr_RightHandSide(optAssem)
% function R = Unstr_RightHandSide(optAssem)
%
% This problem compute the Right Hand Side for a Poisson's PDE
%
% Input:
%     nPts     : number of nodes
%     e_nears  : list of neighbors to each element
%     dp_s     : shape functions gradient
%     w_samples: gauss points
%
% Output
%     R = Right Hand Side
%

nPts    = optAssem.nPts   ;
e_nears = optAssem.e_nears;
p_samp  = optAssem.p_samp;
w_samp  = optAssem.w_samp ;
x_samp  = optAssem.x_samp ;

nElem   = length(e_nears);
sPts    = length(w_samp);


%% ------------------------------------------------------------------------
R = zeros(nPts,1);

[~,~,lf] = Unstr_AnalyticalSolution(x_samp);

gPts = sPts/nElem;

k = 1;
for e=1:nElem
  k_near = e_nears{e};
  nn_k   = length(k_near);

  Rs = zeros(nn_k,1);
  for g=1:gPts
    p_k  = p_samp{k};
    Rs   = Rs - w_samp(k)*lf(k)*p_k;
    
    k  = k+1;
  end
  % Assembly of Ks in K  
  R(k_near) = R(k_near) + Rs;
end