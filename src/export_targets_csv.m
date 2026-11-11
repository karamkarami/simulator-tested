function export_targets_csv(targets, file)
%EXPORT_TARGETS_CSV  Write targets (make_target.m) for radar_generate_ddc --targets FILE.
%   EXPORT_TARGETS_CSV(targets, file) writes one line per target:
%   name,rcs_m2,x,y,z,vx,vy,vz followed by its maneuver rows
%   (tStart_s,turnRate_dps,accelAlong_mps2,accelUp_mps2 each): the position and
%   velocity at t = 0 (x East, y North, z up, from the ground below the radar).
%   The C generator then simulates exactly these targets (e.g. a random
%   scenario drawn in MATLAB).
fid = fopen(file, 'w');
if fid < 0, error('export_targets_csv: cannot write %s', file); end
c = onCleanup(@() fclose(fid));
fprintf(fid, '# name,rcs_m2,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps[,tStart_s,turnRate_dps,accelAlong_mps2,accelUp_mps2]...\n');
for i = 1:numel(targets)
    t = targets(i);
    v = [t.rcs_m2, t.pos0_m(:).', t.vel_mps(:).'];
    if isfield(t, 'maneuvers') && ~isempty(t.maneuvers)
        m = t.maneuvers.';
        v = [v, m(:).']; %#ok<AGROW>
    end
    fprintf(fid, '%s', t.name);
    fprintf(fid, ',%.17g', v);
    fprintf(fid, '\n');
end
end
