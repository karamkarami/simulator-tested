function ok = test_generate_ddc()
%TEST_GENERATE_DDC  DDC files from generate_ddc_files.m (parfor over the PRIs).
%   1. Without thermal and ADC noise every PRI is deterministic: the files of
%      generate_ddc_files must be bit-identical to those that run_radar_sim
%      records (cfg.io.saveDdc), header included (48 elements, clutter, one
%      target, 2 MHz and 4 MHz DDC output).
%   2. With noise, a second run over a part of the sector must give the same
%      files for the pulses it shares with the first run: the random numbers
%      of a pulse depend only on its number, not on the order of the pulses
%      or the number of workers.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
out = fullfile(root, 'output', 'test_generate_ddc');
ok = true;

cfg = radar_config();
cfg.display.enable  = false;
cfg.sim.verbose     = false;
cfg.sim.sectors_deg = [10 10.25];
cfg.targets         = make_target('A', 10, 20e3, 10.1, 1000, 0, 0);
cfg.io.saveRaw      = false;
cfg.io.outFolder    = out;

for rate = [2e6 4e6]
    c = cfg;
    c.ddc.outputRate_Hz = rate;
    c.rx.addThermalNoise = false;
    c.adc.addAdcNoise = false;
    c.io.saveDdc = true;
    c.io.ddcFolder = fullfile(out, 'sim');
    clean(c.io.ddcFolder);
    res = run_radar_sim(c);
    c.io.saveDdc = false;
    c.io.ddcFolder = fullfile(out, 'gen');
    clean(c.io.ddcFolder);
    info = generate_ddc_files(c, 0);
    [same, n] = same_files(fullfile(out, 'sim'), c.io.ddcFolder, true);
    good = same && n == res.nProcessed && n == info.nFiles && n > 0;
    fprintf('  %g MHz, no noise: %d files, identical to run_radar_sim: %s\n', rate/1e6, n, ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end

c = cfg;                                              % with noise: a pulse does not depend on the others
c.io.ddcFolder = fullfile(out, 'all');
clean(c.io.ddcFolder);
infoAll = generate_ddc_files(c, 0);
c.sim.sectors_deg = [10.12 10.25];
c.io.ddcFolder = fullfile(out, 'part');
clean(c.io.ddcFolder);
infoPart = generate_ddc_files(c, 0);
[same, n] = same_files(fullfile(out, 'all'), c.io.ddcFolder, false);
d = load(fullfile(c.io.ddcFolder, sprintf('ddc_%07d.mat', infoPart.pulses(1))));
e = load(fullfile(out, 'all', sprintf('ddc_%07d.mat', infoAll.pulses(1))));
noisy = std(double(d.I(:))) > 1 && ~isequal(d.I, e.I);   % noise present, pulses differ
good = same && n == infoPart.nFiles && n > 0 && n < infoAll.nFiles && noisy;
fprintf('  with noise: %d of %d pulses generated again, identical: %s\n', n, infoAll.nFiles, ternary(good, 'ok', 'FAIL'));
ok = ok && good;

rmdir(out, 's');
fprintf('DDC file generation (parfor): %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function clean(folder)
if exist(folder, 'dir'), rmdir(folder, 's'); end
end

function [same, n] = same_files(refFolder, folder, withHeader)
% Every file of folder must exist in refFolder with the same I, Q (and header).
files = dir(fullfile(folder, 'ddc_*.mat'));
n = numel(files);
same = true;
for i = 1:n
    a = load(fullfile(refFolder, files(i).name));
    b = load(fullfile(folder, files(i).name));
    same = same && isequal(a.I, b.I) && isequal(a.Q, b.Q) && (~withHeader || isequal(a.hdr, b.hdr));
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
