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
for i = 1:numel(reports)
    rp = reports(i);
    j = find_group(g, rp, ext.nDop);
    if isempty(j)
        g(end+1) = struct('beam', rp.beam, 'rIdx', rp.rIdx, 'dIdx', rp.dIdx, ...
                          'hits', 1, 'miss', 0, 'best', rp); %#ok<AGROW>
        updated(end+1) = true; %#ok<AGROW>
    else
        g(j).beam = rp.beam;  g(j).rIdx = rp.rIdx;  g(j).dIdx = rp.dIdx;
        g(j).hits = g(j).hits + 1;
        g(j).miss = 0;
        if rp.power > g(j).best.power
            g(j).best = rp;
        end
        updated(j) = true;
    end
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

function j = find_group(g, rp, nDop)
% Index of the open group that matches report rp (empty if none).
j = [];
for k = 1:numel(g)
    dDop = abs(g(k).dIdx - rp.dIdx);
    dDop = min(dDop, nDop - dDop);                     % circular Doppler distance
    if abs(g(k).beam - rp.beam) <= 1 && abs(g(k).rIdx - rp.rIdx) <= 2 && dDop <= 1
        j = k;
        return
    end
end
end
