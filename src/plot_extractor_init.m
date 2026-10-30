function ext = plot_extractor_init(nDoppler)
%PLOT_EXTRACTOR_INIT  Empty plot extractor (see plot_extractor.m).
ext.groups = struct('beam', {}, 'rIdx', {}, 'dIdx', {}, 'hits', {}, 'miss', {}, 'best', {});
ext.nDop   = nDoppler;
end
