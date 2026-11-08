function ok = test_mti_fir()
%TEST_MTI_FIR  Clutter filter option cfg.mti.type: canceller or FIR filter from a .mat file.
%   The default FIR file (make_input_files.m) is read; its response has a
%   deeper and wider notch around zero Doppler than the three-pulse
%   canceller; the noise gain of every Doppler bin follows from the
%   coefficients (flat noise floor after the division); one row of
%   coefficients is used for every stagger position; a file of your own
%   (any length, one row per stagger position) is used as given; invalid
%   files are refused; export_receiver_config writes the coefficients and
%   mtiType.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
cfg = radar_config();
cfg.clutter.enable = false;
ok = true;
sysC = setup_radar(cfg);
cfg.mti.type = 'fir';
sysF = setup_radar(cfg);
fir = load(cfg.files.mtiFir);
ok = report(isequal(sysF.mtiCoeffs, fir.coeffs) && strcmp(sysF.mti.type, 'fir') && strcmp(sysC.mti.type, 'canceller') ...
            && isequal(sysC.mtiCoeffs, [1 -2 1]), 'FIR coefficients of cfg.files.mtiFir instead of [1 -2 1]') && ok;
ok = report(sysF.cpi.rowLen == sysC.cpi.rowLen - 2, 'the CPI loses L - 1 = 4 pulses to the 5-tap filter') && ok;
% response (relative to white noise) at 2% and 5% of the PRF (moving clutter) and at PRF/2
H = @(h, f) abs(sum(h .* exp(-2j*pi*f*(0:numel(h)-1))))^2 / sum(h.^2);
ok = report(H(fir.coeffs, 0.02) < H([1 -2 1], 0.02) / 1e3 && H(fir.coeffs, 0.05) < H([1 -2 1], 0.05) / 1e3 ...
            && H(fir.coeffs, 0.5) >= H([1 -2 1], 0.5), ...
            sprintf('notch: %.1f dB at 0.02 PRF (canceller %.1f dB)', 10*log10(H(fir.coeffs, 0.02)), ...
            10*log10(H([1 -2 1], 0.02)))) && ok;
% white noise through MTI + Doppler: flat noise floor after the noise gain
rng(1);
n = sysF.cpi.length;
x = complex(randn(1, 2000, n), randn(1, 2000, n));
P = doppler_process(mti_filter(x, sysF.mtiCoeffs, mod(0:n-1, sysF.nStagger) + 1), sysF);
lvl = squeeze(mean(P, 2));
tested = abs((1:sysF.dop.nfft) - (floor(sysF.dop.nfft/2) + 1)) > sysF.dop.notchBins;
ok = report(max(abs(10*log10(lvl(tested) / 2))) < 0.5, 'flat noise floor after the noise gain of every bin') && ok;
% staggered PRF: the one row is used for every stagger position
cfg.prf.mode = 'staggered';
sysS = setup_radar(cfg);
ok = report(isequal(sysS.mtiCoeffs, repmat(fir.coeffs, sysS.nStagger, 1)), 'staggered PRF: the row for every stagger position') && ok;
% a file of your own: one row per stagger position; invalid files
own = [tempname() '.mat'];
cfg.files.mtiFir = own;
coeffs = [1 -1; 1 -1; 2 -2; 1 -1];
save(own, 'coeffs');
sysO = setup_radar(cfg);
ok = report(isequal(sysO.mtiCoeffs, coeffs) && sysO.cpi.length == sysS.cpi.length, ...
            'own file: 2 taps, one row per stagger position') && ok;
bad = {[1 2; 3 4; 5 6], 'rows', [1 1j], 'real', [0 0 0], 'zero'};
for i = 1:2:numel(bad)
    coeffs = bad{i};
    save(own, 'coeffs');
    try
        setup_radar(cfg);
        refused = false;
    catch e
        refused = ~isempty(strfind(lower(e.message), bad{i+1}));
    end
    ok = report(refused, sprintf('refused: coefficients (%s)', bad{i+1})) && ok;
end
delete(own);
% export
cfg = radar_config();
cfg.mti.type = 'fir';
file = [tempname() '.bin'];
export_receiver_config(cfg, file);
fid = fopen(file, 'r', 'ieee-le');  raw = fread(fid, inf, 'uint8=>uint8').';  fclose(fid);
delete(file);
ok = report(~isempty(strfind(char(raw), 'mtiType')) && ~isempty(strfind(char(raw), 'mtiCoeffs')), ...
            'export_receiver_config writes mtiCoeffs and mtiType') && ok;
fprintf('Clutter filter option (canceller / FIR): %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function ok = report(cond, msg)
ok = logical(cond);
fprintf('  %-62s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
