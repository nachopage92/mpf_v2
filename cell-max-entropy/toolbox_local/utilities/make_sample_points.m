function [xsample,x_samp,y_samp]=make_sample_points(x_a,blk_spacing,samp_frac,tri,inside)

x_samp=[min(x_a(:,1)):blk_spacing/samp_frac:max(x_a(:,1))]+blk_spacing/samp_frac/2;
y_samp=[min(x_a(:,2)):blk_spacing/samp_frac:max(x_a(:,2))]+blk_spacing/samp_frac/2;
n_samp=0;
for i=1:(length(x_samp))
    for j=1:(length(y_samp))
        x_trial=[x_samp(i),y_samp(j)];
        el1 = pointLocation(tri, x_trial);
        if ~isnan(el1)
        if inside(el1) == 1
        %[k_tr,area_tr]=convhull([boundary(:,1);x_trial(1)], ...
        %    [boundary(:,2);x_trial(2)]);
        %if (area_tr<=area_bd)
            n_samp=n_samp+1;
            xsample(n_samp,1:2)=x_trial;
        end
        end
    end
end
