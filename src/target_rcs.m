function rcs = target_rcs(tgt, look)
%TARGET_RCS  RCS of a target during one antenna pass over it.
%   With Swerling 1 (see draw_target_rcs.m) every pass (look, see
%   target_look.m) has its own RCS, constant during the dwell; otherwise the
%   mean RCS tgt.rcs_m2 is used.
if isfield(tgt, 'rcsScan') && ~isempty(tgt.rcsScan)
    rcs = tgt.rcsScan(min(look, numel(tgt.rcsScan)));
else
    rcs = tgt.rcs_m2;
end
end
