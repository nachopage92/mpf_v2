function [normL2,normE] = ErrorNorms(u_num,optAssem)
% function [normL2,normE] = ErrorNorms(u_num,optAssem)
% normL2: L2 nor of the error
% normE : Energy semi-norm

e_nears = optAssem.e_nears ;
p_samp  = optAssem.p_samp;
dp_samp = optAssem.dp_samp;
w_samp  = optAssem.w_samp ;
x_samp  = optAssem.x_samp ;

sPts = length(w_samp);

nElem  = length(e_nears);
gPts   = sPts/nElem;


%Displacement and stress tensor exact solutions for the Timoshenko's cantilever beam
[u_sol,du_sol] = Unstr_AnalyticalSolution(x_samp);

errL2 = 0;
errE  = 0;

k = 1;
for e =1:nElem
  k_near = e_nears{e};
  u_k    = u_num(k_near);
  
  for g=1:gPts
    %L2 nor of the error
    p_k   = p_samp{k};
    u_h   = u_sol(k) - p_k'*u_k;
    errL2 = errL2 + (u_h)^2 * w_samp(k);
    
    %Energy semi-norm
    dp_k  = dp_samp{k};   
    du_h  = dp_k'*u_k-du_sol(k,:)';   
    errE  = errE + (du_h'*du_h) * w_samp(k);
    
    k = k + 1;
  end
end
normL2 = sqrt(errL2);
normE  = sqrt(0.5*errE);