clear all
close all


% Node set, and some pre-processing
data = load(strcat('circle24_mesh.mat'));

x_n  = data.x_n;

dx   = sqrt(sum(x_n.*x_n,2));

x_n  = x_n(abs(dx-1)>0.1,:);
dx   = sqrt(sum(x_n.*x_n,2));

x_n  = x_n/max(dx);
nPts = size(x_n,1);
ids  = 1:nPts;

tri  = delaunay(x_n);
ebnd = freeBoundary(triangulation(tri,x_n));
ibnd = unique(ebnd(:));

id1 = [7;18];
id2 = [6 12;12 25];
x_m = zeros(2,2);
x_m(1,:) = mean(x_n(id2(1,:),:));
x_m(2,:) = mean(x_n(id2(2,:),:));

x_n(id1,:) = x_m;

%We define a non-convex domain
figure(2);clf
plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','b','Markersize',10 )
hold on
%We define a non-convex domain
for k = 1:size(ebnd,1)
  plot(x_n(ebnd(k,:),1),x_n(ebnd(k,:),2),'r-','LineWidth',4)
end
plot(x_n(ibnd,1),x_n(ibnd,2),'ko','MarkerFaceColor','c','Markersize',8 )
plot(x_m(:,1),x_m(:,2),'ks','MarkerFaceColor','g','Markersize',12 )
%Triangulation constrained by the boundary edges
triplot(tri,x_n(:,1),x_n(:,2),'m-','LineWidth',1);
hold off
axis equal

save('prior_mesh.mat','x_n','tri','ibnd');