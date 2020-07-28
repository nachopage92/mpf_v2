function plot_domain(x_bd,x_n,ebnd)

figure
axis off
plot(x_bd(:,1),x_bd(:,2),'k.','LineWidth',1.5)
hold on
plot(x_n(:,1),x_n(:,2),'ko','MarkerFaceColor','k','Markersize',4)
axis equal

if nargin>2
  for k = 1:size(ebnd,1)
    plot(x_n(ebnd(k,:),1),x_n(ebnd(k,:),2),'k-','LineWidth',2)
  end
end