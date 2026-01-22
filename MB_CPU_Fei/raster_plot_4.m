function [] = raster_plot_4(sptimes_PN, sptimes_KC, sptimes_EN)
%RASTER_PLOT Produces a raster plot of the input.

[no_PN, time] = size(sptimes_PN);
[no_KC] = size(sptimes_KC, 1);
[no_EN] = size(sptimes_EN, 1);
raster_PN = repmat([1:no_PN]', 1, time).*sptimes_PN;
raster_KC = repmat([no_PN+1:no_PN+no_KC]', 1, time).*sptimes_KC;
raster_EN = repmat([no_PN+no_KC+1:no_PN+no_KC+no_EN]', 1, time).*sptimes_EN;
plot(1:time, raster_EN, 'x','MarkerEdgeColor', 'r', 'MarkerSize', 5);
hold on;
plot(1:time, raster_PN, 'g.', 'MarkerSize', 0.1);
plot(1:time, raster_KC, 'b.', 'MarkerSize', 0.1);
hold off;
end

