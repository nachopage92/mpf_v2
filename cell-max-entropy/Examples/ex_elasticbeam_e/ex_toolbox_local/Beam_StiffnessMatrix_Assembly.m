function K = Beam_StiffnessMatrix_Assembly(optAssem)
% function K = Beam_StiffnessMatrix_Assembly(optAssem)
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

% Material parameters
E       = optAssem.E;
nu      = optAssem.nu;

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

K = spalloc(2*nPts,2*nPts,2*nn*nPts);
C_stiff=E/(1+nu)/(1-2*nu)*[1-nu,   nu,          0 ;...
	                           nu, 1-nu,          0 ;...
	                            0,    0, (1-2*nu)/2];

B_k_ = zeros(3,2*nn);
id1_ = 1:2:2*nn;
id2_ = 2:2:2*nn;

gPts = sPts/nElem;

k = 1;
for e=1:nElem
  k_near = e_nears{e};
  nn_k   = length(k_near);

  B_k = B_k_(:,2*nn_k);
  id1 = id1_(1:nn_k);
  id2 = id2_(1:nn_k);

  Ks = zeros(nn_k*2,nn_k*2);
  for g=1:gPts
    dp_k       = dp_samp{k};
    B_k(1,id1) = dp_k(:,1)';
    B_k(2,id2) = dp_k(:,2)';
    B_k(3,id2) = dp_k(:,1)';
    B_k(3,id1) = dp_k(:,2)';

    Ks = Ks + w_samp(k) * (B_k'*(C_stiff*B_k));
    
    k  = k+1;
  end
  % Assembly of Ks in K
  ipoz = [2*k_near-1;2*k_near];
  
  K(ipoz(:),ipoz(:)) = K(ipoz(:),ipoz(:)) + Ks;
end