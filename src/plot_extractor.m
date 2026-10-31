function [ext, plots] = plot_extractor(ext, reports, opt)
%PLOT_EXTRACTOR  Azimuth peak picking over successive (sliding) CPIs.
%   [ext, plots] = PLOT_EXTRACTOR(ext, reports, opt) groups the detections of
%   the same target found in successive CPIs while the beam sweeps over it
%   (same beam +-1, range cell +-2, Doppler bin +-1). When a group receives no
%   detection for more than opt.maxMiss CPIs it is closed and its strongest
%   detection is returned as a plot (if the group has at least opt.minHits
%   detections; plot.hits = number of detections). Set opt.flush = true to
%   close all open groups (reports = []).
%   Create an empty extractor with plot_extractor_init.

g = ext.groups;
updated = false(1, numel(g));
% beam, range cell and Doppler bin of every open group (kept up to date below)
gb = [g.beam];  gr = [g.rIdx];  gd = [g.dIdx];
for i = 1:numel(reports)
    rp = reports(i);
    dDop = abs(gd - rp.dIdx);
    dDop = min(dDop, ext.nDop - dDop);                 % circular Doppler distance
    j = find(abs(gb - rp.beam) <= 1 & abs(gr - rp.rIdx) <= 2 & dDop <= 1, 1);   % first matching group
    if isempty(j)
        g(end+1) = struct('beam', rp.beam, 'rIdx', rp.rIdx, 'dIdx', rp.dIdx, ...
                          'hits', 1, 'miss', 0, 'best', rp); %#ok<AGROW>
        updated(end+1) = true; %#ok<AGROW>
        j = numel(g);
    else
        g(j).beam = rp.beam;  g(j).rIdx = rp.rIdx;  g(j).dIdx = rp.dIdx;
        g(j).hits = g(j).hits + 1;
        g(j).miss = 0;
        if rp.power > g(j).best.power
            g(j).best = rp;
        end
        updated(j) = true;
    end
    gb(j) = rp.beam;  gr(j) = rp.rIdx;  gd(j) = rp.dIdx;
end

% Groups without a new detection
for j = find(~updated)
    g(j).miss = g(j).miss + 1;
end
closed = [g.miss] > opt.maxMiss;
if opt.flush
    closed(:) = true;
end
plots = [];
for j = find(closed & [g.hits] >= opt.minHits)
    p = g(j).best;
    p.hits = g(j).hits;                                % number of CPIs with a detection
    plots = [plots, p]; %#ok<AGROW>
end
ext.groups = g(~closed);
end
