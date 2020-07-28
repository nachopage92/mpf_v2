function [p,dp,q,dq,i_fail]=shapef_elem(x_n,x_s,options)
%function [p,dp,q,dq,i_fail]=shapef_elem(TolNR,enear,bd_segm,bd_segm_out,bd_nodes,x_n,x_s,nDim,apower,bpower)
% Calculation of the local maximum-entropy shape functions with Newton's
% method as described in section 4.2 of [1]
%
% INPUT
% =====
% x_n:     nodal set (nPts, nDim)
% x_s:     sample points (sPts, nDim) in an element where the shape functions are evaluated 
% enear:   contains the list of nodes with non-zero shape funcion affecting this element
% TolNR:   tolerance for Newton-Raphson's iterations
%
% OUTPUT
% ======
% p:        p{i} contains the values of the shape functions corresponding
%           to the nodes in enear at the i-th sample point
% dp:       dp{i} contains the values of the nDim spacial derivatives of
%           the shape functions corresponding to the nodes in enear at
%           the i-th sample point
%
% Reference:
% [1] Marino Arroyo and Michael Ortiz, Local maximum-entropy approximation
%     schemes: a seamless bridge between finite elements and meshfree methods,
%     International Journal for Numerical Methods in Engineering, 65:2167-2202 (2006).
warning('off')

%% Basic input checkings and variable definitions ========================================
if isfield(options,'TolNR')
  TolNR   = options.TolNR;
else
  error('Tolerance TolNR for Newton-Raphson algorithm is not defined')
end

if isfield(options,'enear')
  enear   = options.enear;
else
  error('List of nodes enear affecting the element is not defined')
end

[nPts,nDim] = size(x_n);
[sPts,sDim] = size(x_s);     %number of sample points in the element
if nPts<nDim+1
  error('The number of node points is bad set');
end
if nDim ~= sDim
  error('The node points and the samples have different dimension');
end

eNN      = length(enear);   %number of nearest neighbors affecting this element
i_fail   = zeros(sPts,1);
max_iter = 100;

%% we compute the max-ent basis functions ================================================
p   = zeros(sPts,eNN);
dp  = zeros(sPts,eNN,nDim);
aux = zeros(eNN,nDim);

q  = options.q;
dq = options.dq;

for k = 1:sPts
  
  %[2]-------------------------------------------------------------------------
  %Newton's method
  x  = x_s(k,:);
  lam= zeros(1,nDim);
  R  = 10*ones(1,nDim);
  niter=0;
  
  %[2.1]---------------------------------------------------
  %Newton iteration: Regular iteration
  while (norm(R)>TolNR) && niter <= max_iter
    [~,R,J,p_a]=Gamma_(nDim,x_n,x,lam,enear,q(k,:));
    if (abs(rcond(J)))<1e-8
      %disp('Newton Failed, near to singular matrix ')
      %plot(x(1),x(2),'rx','Markersize',10 ,'LineWidth',2);
      %pause(0.01)
      p(k,:)=NaN([eNN 1]);
      i_fail(k)=1;
      break
    end
    dlam=-J\R';
    if sum(dlam.*R')>0
      %disp('Newton Failed, near to ascent direction')
      %plot(x(1),x(2),'mx','Markersize',10 ,'LineWidth',2);
      %pause(0.01)
      p(k,:)=NaN([eNN 1]);
      i_fail(k)=1;
      break
    end
    lam=lam+dlam';
    niter=niter+1;
  end
  if niter >=max_iter
    %disp('Newton Failed, max num iterations')
    %plot(x(1),x(2),'cx','Markersize',10 ,'LineWidth',2);
    %pause(0.01)
    p(k,:)=NaN([eNN 1]);
    i_fail(k)=1;
  end
  if any(isnan(p_a))
    %disp('Newton Failed, NaNs')
    %plot(x(1),x(2),'cx','Markersize',10 ,'LineWidth',2);
    %pause(0.01)
    p(k,:)=NaN([eNN 1]);
    i_fail(k)=1;
  end
  %disp('Newton iterations')
  %disp(niter)
  
  %[2.2]---------------------------------------------------
  % Newton with Line-Search when things go wrong ...
  if i_fail(k) == 1
    i_fail(k)=0;
    %disp('We try with line-search (LS)')
    lam=zeros(1,nDim);
    R=10*ones(1,nDim);
    niter=0;
    while (norm(R)>TolNR) && niter <= max_iter
      [~,R,J,p_a]=Gamma_(nDim,x_n,x,lam,enear,q(k,:));
      dlam=-J\R';
      %Line-Search
      if sum(dlam.*R')>0
        dlam=-dlam;
      end
      opts=optimset('TolX',max(1e-10,TolNR*100),'MaxIter',5000);
      t = fminbnd(@(t) Gamma(nDim,x_n,x,enear,q(k,:),lam,dlam,t),0,1,opts);
      lam=lam+t*dlam';
      niter=niter+1;
    end
    
    if niter >=max_iter
      fprintf(1,'Newton-LS Failed, max num iterations=%d  |R|=%12.2e  TolNR=%12.2e\n',niter,norm(R),TolNR);
      p(k,:)=NaN([eNN 1]);
      i_fail(k)=1;
      
      figure(50);clf
      hold on
      if nDim==1
        plot(x_n,0*x_n,'ro')
        plot(x_n(enear),0*x_n(enear),'g*',x,0,'mx','Markersize',10 ,'LineWidth',2);
      else
        plot(x_n(:,1),x_n(:,2),'ro')
        plot(x_n(enear,1),x_n(enear,2),'g*',x(1),x(2),'mx','Markersize',10 ,'LineWidth',2);
      end
      hold off
      axis equal
      pause(0.1)
    end
    
    if any(isnan(p_a))
      disp('Newton-LS Failed, NaNs')
      p(k,:)=NaN([eNN 1]);
      i_fail(k)=1;

      figure(60);clf
      hold on
      if nDim==1
        plot(x_n,0*x_n,'ro')
        plot(x_n(enear),0*x_n(enear),'g*',x,0*x,'cx','Markersize',10 ,'LineWidth',2);
      else
        plot(x_n(:,1),x_n(:,2),'ro')
        plot(x_n(enear,1),x_n(enear,2),'g*',x(1),x(2),'cx','Markersize',10 ,'LineWidth',2);
      end
      hold off
      axis equal
      pause(0.1)
    end
    
%     fprintf(1,'Newton-LS iterations = %d\n',niter);
  end
  
  %[3]-------------------------------------------------------------------------
  %Checking fails
  if i_fail(k) == 0
    p(k,:)=p_a;
    x_m_xa=zeros(nDim,eNN);
    for kk = 1:nDim
      x_m_xa(kk,:) = x(kk)-x_n(enear,kk);
    end
    aux0 = exp(lam*x_m_xa)';
    aux0 = aux0/sum(q(k,:)'.*aux0);
    dr_dx = zeros(nDim,nDim);
    for kk = 1:nDim
      for jj = 1:nDim
        dr_dx(kk,jj)= sum(aux0.*dq(k,:,jj)'.*x_m_xa(kk,:)');
      end
    end
    dr_dx=dr_dx+eye(nDim);
    dlamdx = -inv(J)*dr_dx ;
    aux1 = dlamdx'*x_m_xa;
    for kk = 1:nDim
      aux2 = sum(aux0.*dq(k,:,kk)');
      aux(:,kk) = aux0.*dq(k,:,kk)' + p_a.*aux1(kk,:)' - p_a * aux2;
    end
    dp(k,:,:) = aux;
    %dp(k,:,:) = aux*scale;
    %dq(k,:,:) = dq(k,:,:)*scale;
  end
end

