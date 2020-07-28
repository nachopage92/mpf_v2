function K = Unstr_StiffnessMatrix_Assembly(optAssem)
% function K = Unstr_StiffnessMatrix_Assembly(optAssem)
%
% This problem compute the stiffness matrix k for a Poisson's PDE
%
% Input:
%     nPts     : number of nodes
%     e_nears  : list of neighbors to each element
%     dp_s     : shape functions gradient
%     w_samples: gauss points
%
% Output
%     K: stiffness matrix
%

nPts    = optAssem.nPts   ;
e_nears = optAssem.e_nears;
dp_samp = optAssem.dp_samp;
w_samp  = optAssem.w_samp ;

nElem   = length(e_nears);
sPts    = length(w_samp);


%% ------------------------------------------------------------------------
%  The stiffness matrix is assembled 
nn_e = zeros(nElem,1);
for e=1:nElem
  nn_e(e)= length(e_nears{e});
end
nn = max(nn_e);

K = spalloc(nPts,nPts,nn*nPts);

gPts = sPts/nElem;

k = 1;
for e=1:nElem
  k_near = e_nears{e};
  nn_k   = length(k_near);

  Ks = zeros(nn_k,nn_k);
  for g=1:gPts
    dp_k = dp_samp{k};
    Ks   = Ks + w_samp(k)*(dp_k*dp_k');
    
    k  = k+1;
  end
  % Assembly of Ks in K  
  K(k_near,k_near) = K(k_near,k_near) + Ks;
end