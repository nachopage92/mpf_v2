function K = StiffnessMatrix_Assembly(optAssem)
% function K = StiffnessMatrix_Assembly(optAssem)
%
% This problem compute the stiffness matrix k for a Poisson's PDE
%
% Input:
%     nPts     : number of nodes
%     s_nears  : list of neighbors
%     dp_s     : shape functions gradient
%     w_samples: gauss points
%
% Output
%     K: stiffness matrix
%
nPts    = optAssem.nPts   ;
s_near  = optAssem.s_near;
dp_samp = optAssem.dp_samp ;
w_samp  = optAssem.w_samp ;


sPts = length(w_samp);

nn=0;
for k=1:sPts
  nn = max(nn, length(s_near{k}));
end
nn = min(nn,nPts);

K = spalloc(nPts,nPts,nn*nPts);

%% ------------------------------------------------------------------------
for k=1:sPts
  k_near = s_near{k};
  dp_k   = dp_samp{k};
  
  K(k_near,k_near) = K(k_near,k_near) + w_samp(k)*(dp_k*dp_k');
end
