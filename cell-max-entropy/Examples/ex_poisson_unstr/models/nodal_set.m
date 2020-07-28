clear all
close all


%N_ = [81 289 595 1135 2147 4303 8532 16868 33382 66546];
N_ = [81 120 283 581 1121 2132 4319 8510 16834];

for n = 1:length(N_)
    
    filename = sprintf('nodes_N%05d.txt',N_(n));
    x_n = load(filename);
    
    size(x_n)
    
    p_src = [0.51 0.52;0.31 0.34; 0.73 0.71; 0.28 0.72];
    
    figure(1)
    plot(x_n(:,1),x_n(:,2),'ro',p_src(:,1),p_src(:,2),'k*')
    axis equal
    pause
end