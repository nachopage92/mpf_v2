%%
% This function plots the prior functions based in R-functions
% for a given node set, as well as the gradients and the approximation
% for discontinuous nodal data
%

clear all
close all
format long

% ========================================================================================
%% Input parameters

N    = 100;

x_n  = [-1  10;-1 -10;
        1 -10; 1  10];  
edges= [1 2;3 4];

x_s  = [linspace(-1,1,201)',zeros(201,1)];

figure(1)
plot(x_n(:,1),x_n(:,2),'ro','markersize',12,'linewidth',2)
hold on
plot(x_s(:,1),x_s(:,2),'bx')
for ed=1:size(edges,1)
  plot(x_n(edges(ed,:),1),x_n(edges(ed,:),2),'r*-','markersize',12,'linewidth',2)
end
hold off
axis equal

%% RFunction - ------------------------------------------------------

%  Input parameters: CME options
% spow: w^spow, where w is the approximation to the distance function (R-function)
% mpow: distance function derivative maximum degree of the approximation
spow1 = 1; 
spow2 = 2;
spow3 = 3;
mpow1 = 1; 
mpow2 = 2;
mpow3 = 3;

[ddist11,gdist11] = Rfunction_equiv_(edges,x_n,x_s,spow1,mpow1);
[ddist12,gdist12] = Rfunction_equiv_(edges,x_n,x_s,spow1,mpow2);
[ddist13,gdist13] = Rfunction_equiv_(edges,x_n,x_s,spow1,mpow3);
[ddist14,gdist14] = Rfunction_equiv_(edges,x_n,x_s,spow1,10);

[ddist21,gdist21] = Rfunction_equiv_(edges,x_n,x_s,spow2,mpow1);
[ddist22,gdist22] = Rfunction_equiv_(edges,x_n,x_s,spow2,mpow2);
[ddist23,gdist23] = Rfunction_equiv_(edges,x_n,x_s,spow2,mpow3);
[ddist24,gdist24] = Rfunction_equiv_(edges,x_n,x_s,spow2,10);

[ddist31,gdist31] = Rfunction_equiv_(edges,x_n,x_s,spow3,mpow1);
[ddist32,gdist32] = Rfunction_equiv_(edges,x_n,x_s,spow3,mpow2);
[ddist33,gdist33] = Rfunction_equiv_(edges,x_n,x_s,spow3,mpow3);
[ddist34,gdist34] = Rfunction_equiv_(edges,x_n,x_s,spow3,10);

[ddist41,gdist41] = Rfunction_equiv_(edges,x_n,x_s,10,mpow1);
[ddist42,gdist42] = Rfunction_equiv_(edges,x_n,x_s,10,mpow2);
[ddist43,gdist43] = Rfunction_equiv_(edges,x_n,x_s,10,mpow3);
[ddist44,gdist44] = Rfunction_equiv_(edges,x_n,x_s,10,10);


%% Plot basis functions at a given node
createfigure_dist1D(x_s(:,1),(1-abs(x_s(:,1))).^1,ddist11,ddist12,ddist13,ddist14)
title('s^1')

createfigure_dist1D(x_s(:,1),(1-abs(x_s(:,1))).^2,ddist21,ddist22,ddist23,ddist24)
title('s^2')

createfigure_dist1D(x_s(:,1),(1-abs(x_s(:,1))).^3,ddist31,ddist32,ddist33,ddist34)
title('s^3')

createfigure_dist1D(x_s(:,1),(1-abs(x_s(:,1))).^10,ddist41,ddist42,ddist43,ddist44)
title('s^4')