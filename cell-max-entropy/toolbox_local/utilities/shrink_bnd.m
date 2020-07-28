function x_b = shrink_bnd(x_a,ebnd,shrink_fct)

if nargin==2
  shrink_fct = 0.01;
end

x_b  = x_a;
ibnd = unique(ebnd);
nPts = size(x_a,1);
bPts = length(ibnd);

%Global to local index of the boundary nodes
glo2loc = zeros(nPts,1);
for ib=1:bPts
  glo2loc(ibnd(ib)) = ib;
end

%dual list to retrieve the boundary line elements for a given point at the boundary
nEdge=size(ebnd,1);
nbnd = zeros(nEdge,2);
iw = ones(bPts,1);
for eb=1:nEdge
  a = glo2loc(ebnd(eb,1));
  b = glo2loc(ebnd(eb,2));
  nbnd(a,iw(a)) = eb;
  nbnd(b,iw(b)) = eb;
  iw(a) = 2;
  iw(b) = 2;
end

for ib=1:bPts
  n_bd = [0,0]; %normal at the boundary node x_a(ibnd(ib),:)
  
  % for each nearest edge of the ib-th boundary node
  for j=1:2
    eb    = nbnd(ib,j);           %node's ID the edge
    x_eb1 = x_a(ebnd(eb,1),:);
    x_eb2 = x_a(ebnd(eb,2),:);
    elng   = norm(x_eb2-x_eb1);   %length of the tedge
    
    %vector along the edge (tangent to the BD)
    eb_t = [x_eb2(1)-x_eb1(1),x_eb2(2)-x_eb1(2)]/elng;
    
    %vector normal to the BD
    eb_n = [-eb_t(2),eb_t(1)];
    eb_z = cross([eb_t 0],[eb_n 0]);
    if eb_z(3)<0
      eb_n = [eb_t(2),-eb_t(1)];
      eb_z = cross([eb_t 0],[eb_n 0]);
      if eb_z(3)<0
        disp(cross([eb_t 0],[eb_n 0]))
        error('Bad ordering : cross product in the -z direction')
      end
    end
    n_bd = n_bd + eb_n;
  end
  
  n_bd    = n_bd/norm(n_bd);
  
  if abs(n_bd(1))>1e-6
    n_bd(1) = n_bd(1)/abs(n_bd(1));
  end
  if abs(n_bd(2))>1e-6
    n_bd(2) = n_bd(2)/abs(n_bd(2));
  end
  
%   figure(20);clf
%   plot(x_a(:,1),x_a(:,2),'r.')
%   hold on
%   plot(x_a(ibnd(ib),1),x_a(ibnd(ib),2),'k*','markersize',10,'linewidth',2)
%   quiver(x_a(ibnd(ib),1),x_a(ibnd(ib),2),n_bd(1),n_bd(2))
%   hold off
%   axis equal
%   pause
  x_b(ibnd(ib),:) = x_a(ibnd(ib),:) + n_bd*shrink_fct;
end
