function plotCME2D(x_n,u_n, nSX, nSY, Lx, Ly, centre, options)
% function PlotCME2D(x_n,u_n, nSX, nSY, Lx, Ly, centre, options)
% this functios plots the LME's shape functions in 2D
% x_n: nodes coordinates
% u_n: scalar field for each node to be interpolated
% nSX: number of samples in X
% nSY: number of samples in Y
% Lx, Ly are the size of the rectangular region and center its center.

x_s = UniformGrid2D(nSX, nSY, Lx, Ly, centre);
hs    = Lx/(nSX-1);
x_samp= (0:hs:Lx)-Lx/2+centre(1);
y_samp= (0:hs:Ly)-Ly/2+centre(2);
sPts  = length(x_s);

U_mat = zeros(nSX,nSY);
id_smp= (1:sPts);

outCME = wrapper_cme(x_n,x_s,options);
p_samp = outCME.p_samp;
s_near = outCME.s_near;

u_s= zeros(sPts,1);
for k=1:sPts
  u_s(k) = u_n(s_near{k})'*p_samp{k};
end
k=1;
for i=1:nSX
  for j=1:nSY
    ij  = id_smp(k);
    U_mat(i,j)= u_s(ij);
    k=k+1;
  end
end

axis off
plot3(x_n(:,1),x_n(:,2),u_n,'ko','LineWidth',2,...
  'MarkerFaceColor','g','MarkerSize',6);
hold on
surf(x_samp,y_samp,U_mat, ... %'LineWidth',0.5,...
  'FaceColor','red','FaceLighting','phong')
view([10 40])
hold off
%material shiny
%lighting gouraud
lightangle(80, -40)
lightangle(-90, 60)

end