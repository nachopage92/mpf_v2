function [q,dq]=priorf_elem(x_n,x_s,options)
%function [q,dq]=priorf_elem(x_n,x_s,options)
% Calculation of the prior shape functions based in R-functions
%
% INPUT
% =====
% x_n:      nodal set (nPts, nDim)
% x_s  :    sample point where the shape functions are evaluated (in an element)
%
% options: parameters setting
%   options.enear      : contains the list of nearest nodes with non-zero shape funcion affecting this element
%   options.bd_segm    : segments defining the polygon of the last ring, which are not on the boundary
%   options.bd_segm_out: segments defining the polygon of the last ring including segments on the boundary
%   options.bd_nodes   : boundary identifier (0:false, 1:true)
%   options.spow       : w^spow, where w is the approximation to the distance function (R-function)
%   options.mpow       : distance function derivative maximum degree of the approximation
%   options.isConv     : flag indicating if the domain is convex:1 or non convex:0 (for priorfun computation)
%   options.pUnity     : partition of unity is enforced by Shepard's method (1: default,
%                        otherwise is not made)
%
% OUTPUT
% ======
% q:        q{i} contains the values of the prior shape functions corresponding
%           to the nodes in near_{i} at the i-th sample point
% dq:       dp{i} contains the values of the nDim spacial derivatives of
%           the prior shape functions corresponding to the nodes in near_{i} at
%           the i-th sample point
%

%% Basic input checkings and variable definitions ========================================
if isfield(options,'enear')
  enear   = options.enear;
else
  error('List of nodes enear affecting the element is not defined')
end
if isfield(options,'bd_segm') && isfield(options,'bd_segm_out') && isfield(options,'bd_nodes')
  bd_segm    = options.bd_segm;
  bd_segm_out= options.bd_segm_out;
  bd_nodes   = options.bd_nodes;
else
  error('Auxiliary arrays at each node indicating the rings segments and boundary are not defined')
end

if isfield(options,'spow') && isfield(options,'mpow')
  spow     = options.spow;
  mpow     = options.mpow;
else
  error('spow and mpow are not defined')
end

if isfield(options,'isConv')
  isConv   = options.isConv;
else
  error('Flag indicating if the domain is (or not) convex has been not defined.')
end

if isfield(options,'pUnity')
  pUnity   = options.pUnity;
else
  pUnity   = 1;
end


[nPts,nDim] = size(x_n);
[sPts,sDim] = size(x_s);     %number of sample points in the element
if nPts<nDim+1
  error('The number of node points is bad set');
end
if nDim ~= sDim
  error('The node points and the samples have different dimension');
end

eNN = length(enear); %number of nearest neighbors affecting this element


%% we compute the prior ==================================================================
Zq  = zeros(sPts,1);
q   = zeros(sPts,eNN);
dZq = zeros(sPts,nDim);
dq  = zeros(sPts,eNN,nDim);
oneD= ones(1,nDim);

% Here we select the Rfunction to be used (like a pointer to a function)
Rfunction_ = @(mesh,v,p,spow,mpow) Rfunction_equiv_(mesh,v,p,spow,mpow);
%Rfunction_ = @(mesh,v,p,spow,mpow) Rfunction_equiv(mesh,v,p,spow,mpow);
% Rfunction_ = @(mesh,v,p,spow,mpow) Rfunction_conj(mesh,v,p,spow,mpow);

for j=1:eNN  % loop over the nearest nodes affecting the element
  
  bd_segm1 = bd_segm{enear(j)};
  bd_segm2 = bd_segm_out{enear(j)};
  bd_diff  = setdiff(bd_segm2,bd_segm1,'rows');
  
  %   if isempty(bd_diff)==0
  %     bd_segm1
  %     bd_segm2
  %
  %     figure(11);clf
  %     plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','b','Markersize',6 )
  %     hold on
  %     plot(x_s(:,1),x_s(:,2),'bx','Markersize',12,'LineWidth',2)
  %     plot(x_n(enear,1),x_n(enear,2),'bo','MarkerFaceColor','c','Markersize',15,'LineWidth',1)
  %     plot(x_n(enear(j),1),x_n(enear(j),2),'ks','MarkerFaceColor','y','Markersize',15,'LineWidth',2)
  %     for ek=1:size(bd_segm2,1)
  %       plot(x_n(bd_segm2(ek,:),1),x_n(bd_segm2(ek,:),2),'r-','LineWidth',4)
  %     end
  %     for ek=1:size(bd_segm1,1)
  %       plot(x_n(bd_segm1(ek,:),1),x_n(bd_segm1(ek,:),2),'kd-','MarkerFaceColor','g','Markersize',12,'LineWidth',2)
  %     end
  %     hold off
  %     axis equal
  %     pause
  %   end
  
  % if one of the ring edges of an interior node have a segment on the boundary
  if isConv==0 && isempty(bd_diff)==0 && bd_nodes(enear(j))==0
    % prior function computed with Rfunction (equivalence)
    % spow-1 to preserve the "smoothness"
    [phi1,dphi1] = Rfunction_(bd_segm1,x_n,x_s,spow-1,mpow);
    
    [phi2,dphi2] = Rfunction_(bd_segm2,x_n,x_s,1,mpow);
    %[phi2,dphi2] = Rfunction_(bd_diff,x_n,x_s,1,mpow);
    
    phi  = (phi1.*phi2);
    
    if sum(isnan(phi))>0
      phi
      pause
    end
    
    dphi = dphi1.*(phi2*oneD) + (phi1*oneD).*dphi2;
    
    if sum(isnan(dphi(:)))>0
      dphi
      dphi(isnan(dphi))=0;
      pause
    end
  else  % For convex domains or for a full interior node OR for a boundary node in non convexs
    [phi,dphi] = Rfunction_(bd_segm1,x_n,x_s,spow,mpow);
    
    %     x_a = x_n(enear(j),:);
    %     [phi,dphi] = Wfunction_convex(bd_segm1,x_a,x_n,x_s,spow,mpow,options.h);
    
    if sum(isinf(phi(:)))>0
      phi
      figure(21);clf
      plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','b','Markersize',6 )
      hold on
      plot(x_a(1),x_a(2),'kv','Markersize',15,'MarkerFaceColor','m','LineWidth',2)
      plot(x_s(:,1),x_s(:,2),'bx','Markersize',12,'LineWidth',2)
      plot(x_n(enear,1),x_n(enear,2),'bo','MarkerFaceColor','c','Markersize',15,'LineWidth',1)
      plot(x_n(enear(j),1),x_n(enear(j),2),'ks','MarkerFaceColor','y','Markersize',15,'LineWidth',2)
      for ek=1:size(bd_segm1,1)
        plot(x_n(bd_segm1(ek,:),1),x_n(bd_segm1(ek,:),2),'kd-','MarkerFaceColor','g','Markersize',12,'LineWidth',2)
      end
      hold off
      axis equal
      pause
    end
    
    if sum(isnan(dphi(:)))>0
      dphi
      dphi(isnan(dphi))=0;
      pause
    end
  end
  
  Zq        =  Zq +  phi;
  dZq       = dZq + dphi;
  q(:,j)    =  phi;
  dq(:,j,:) = dphi;
end



%% prior functions are made a partition of unity by using Shepard's method
if pUnity == 1
%   disp('Partition of Unity is enforced by using Shepard method')
  aux = zeros(eNN,nDim);
  
  for k = 1:sPts
    q(k,:)=q(k,:)/Zq(k);
    for kk = 1:nDim
      aux(:,kk) = dq(k,:,kk)/Zq(k) - q(k,:)/Zq(k)*dZq(k,kk);
      if sum(isnan(aux(:)))>0
        aux
        pause
      end
    end
    dq(k,:,:)=aux;
  end
end
