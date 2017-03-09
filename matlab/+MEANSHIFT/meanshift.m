function [modes0,modes] = meanshift(u,bandwidth,min_cardinality,min_dist)

h = meanshift('profile','epanechnikov', ...
	         'manifold','s2','bandwidth',bandwidth);

[modes0,modes] = meanshift('est_modes',single(u),single(bandwidth), ...
                           int32(min_cardinality), ...
                           single(min_dist));
modes0 = double(modes0);
modes = double(modes);