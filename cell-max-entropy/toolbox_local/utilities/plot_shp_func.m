function [plot_shp]=plot_shp_func(ishape,nPts,x_samples,js_nears,p_s)

plot_shp=zeros(size(x_samples));
ids = 1:nPts;
for ig=1:length(x_samples)
    if any(js_nears{ig}==ishape)
        ids_i=ids(js_nears{ig}==ishape);
        plot_shp(ig)=p_s{ig}(ids_i);
    end
end
