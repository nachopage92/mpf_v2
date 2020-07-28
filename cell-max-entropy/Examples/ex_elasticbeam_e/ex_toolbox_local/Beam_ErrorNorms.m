function [normL2,normE] = Beam_ErrorNorms(u_h,x_s,optAssem)
% function [normL2,normE] = Beam_ErrorNorms(u_h,x_s,optAssem)
% normL2: L2 nor of the error
% normE : Energy semi-norm

% Material parameters
E       = optAssem.E;
nu      = optAssem.nu;
L       = optAssem.L;
D       = optAssem.D;
P       = optAssem.P;
C_stiff = optAssem.Cmat;
S_stiff = optAssem.Smat;

nPts    = optAssem.nPts   ;
e_nears = optAssem.e_nears ;
p_samp  = optAssem.p_samp;
dp_samp = optAssem.dp_samp;
w_samp  = optAssem.w_samp ;

sPts = length(w_samp);

nElem  = length(e_nears);
gPts   = sPts/nElem;


u_num  = reshape(u_h,2,nPts)';  %in matrix format: nPts x 2

%Displacement and stress tensor exact solutions for the Timoshenko's cantilever beam
[u_sol,sgm_sol] = Beam_AnalyticalSolution(x_s,P,D,L,E,nu);

% nn_s = zeros(sPts,1);
% for k=1:sPts
%   nn_s(k)= length(s_near{k});
% end
% nn = max(nn_s);

nn_e = zeros(nElem,1);
for e=1:nElem
  nn_e(e)= length(e_nears{e});
end
nn = max(nn_e);

B_k_ = zeros(3,2*nn);
id1_ = 1:2:2*nn;
id2_ = 2:2:2*nn;

errL2 = 0;
errE  = 0;

k = 1;
for e =1:nElem
  k_near = e_nears{e};
  u_k    = u_num(k_near,:)';
  nn_k   = length(k_near);
  
  B_k = B_k_(:,2*nn_k);
  id1 = id1_(1:nn_k);
  id2 = id2_(1:nn_k);

  for g=1:gPts
    p_k    = p_samp{k};
    dp_k   = dp_samp{k};

    u_     = (u_k*p_k)';
    
    %L2 nor of the error
    errL2 = errL2 + sum((u_sol(k,:)-u_).^2) * w_samp(k);
    
    %Energy semi-norm
    B_k(1,id1) = dp_k(:,1)';
    B_k(2,id2) = dp_k(:,2)';
    B_k(3,id2) = dp_k(:,1)';
    B_k(3,id1) = dp_k(:,2)';
    
    eps_k = B_k*u_k(:);
    sgm_k = C_stiff*eps_k;
    sgm_  = sgm_sol(k,:)'-sgm_k;
    eps_  = S_stiff*sgm_-eps_k;      %epsilon = S_stiff*sigma
    errE  = errE  + eps_'*sgm_ * w_samp(k);
    
    k = k + 1;
  end
end
normL2 = sqrt(errL2);
normE  = sqrt(0.5*errE);