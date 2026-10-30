function targets = draw_target_rcs(targets, nScans, swerling)
%DRAW_TARGET_RCS  RCS of every target for every antenna pass over it.
%   Swerling 1: the RCS is exponentially distributed around the mean
%   tgt.rcs_m2, constant during a pass (dwell) and independent from pass to
%   pass (scan to scan). swerling = 0: constant RCS. The values are stored
%   in tgt.rcsScan; nScans values are drawn (see target_look.m).
for i = 1:numel(targets)
    if swerling == 1
        targets(i).rcsScan = targets(i).rcs_m2 * -log(rand(1, nScans));
    else
        targets(i).rcsScan = targets(i).rcs_m2 * ones(1, nScans);
    end
end
end
