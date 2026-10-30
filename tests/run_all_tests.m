% RUN_ALL_TESTS  Run every test of the radar simulator and print a summary.
%   The six-target tests take about 35 minutes each in Octave 8.4 on a
%   4-vCPU VM. Each writes about 18 GB of raw data (output/test_<mode>/raw_pulses),
%   checks it and deletes it again (call test_six_targets(mode, true) to keep
%   the files). When the C receiver is built (radar_receiver/build/radar_receiver,
%   Linux or macOS), test_c_receiver compares it with the MATLAB processing
%   (about 40 minutes per PRF mode). All tests together take about 3 hours.
%   The PPI tests need graphics: in Octave without a display run
%   xvfb-run octave --no-gui --eval run_all_tests (octave-cli has no
%   graphics toolkit: the interactive PPI test is skipped, the other tests
%   run without the PPI).

addpath(fileparts(mfilename('fullpath')));
tests = cell(0, 2);                                  % {name, test function} in the order they run

for mode = {'fixed', 'staggered'}
    for pulse = {'short', 'long'}
        tests(end+1, :) = {sprintf('SNR budget (%s PRF, %s pulse)', mode{1}, pulse{1}), ...
                           @() test_snr_budget(mode{1}, pulse{1})}; %#ok<SAGROW>
    end
end
tests(end+1, :) = {'SNR budget (fixed PRF, long pulse, P4 phase code)', ...
                   @() test_snr_budget('fixed', 'long', 'phasecode')};
tests(end+1, :) = {'SNR budget (fixed PRF, long pulse, P4 phase code, mismatch filter)', ...
                   @() test_snr_budget('fixed', 'long', 'phasecode', 48, 2e6, 'mismatch')};
tests(end+1, :) = {'SNR budget (single element, no beamforming)', @() test_snr_budget('fixed', 'short', 'lfm', 1)};
for pulse = {'short', 'long'}
    tests(end+1, :) = {sprintf('SNR budget (4 MHz DDC output, %s pulse)', pulse{1}), ...
                       @() test_snr_budget('fixed', pulse{1}, 'lfm', 48, 4e6)}; %#ok<SAGROW>
end

tests(end+1, :) = {'UDP streaming of the DDC output', @test_udp_stream};
tests(end+1, :) = {'DDC recording and replay over UDP with the PRI timing', @test_ddc_replay};
tests(end+1, :) = {'DDC file generation with parfor (generate_ddc_files)', @test_generate_ddc};
tests(end+1, :) = {'Phase codes: matched and mismatch filters (and the C core test)', @test_phase_code};
tests(end+1, :) = {'Configuration edge cases (sectors, invalid settings, handover)', @test_setup};
tests(end+1, :) = {'CFAR methods (2-D, range, Doppler), edges, range interval', @test_cfar};
tests(end+1, :) = {'EKF tracker 3-D and 2-D (and the C tracker when built)', @test_tracker};
tests(end+1, :) = {'Interactive PPI (data tips, zoom, pan, filters)', @test_ppi_interactive};
tests(end+1, :) = {'No ghost of the short-pulse echo (same IF)', @test_short_pulse_ghost};
for mode = {'fixed', 'staggered'}
    tests(end+1, :) = {sprintf('Six targets with clutter (%s PRF)', mode{1}), ...
                       @() test_six_targets(mode{1}, true, false)}; %#ok<SAGROW>
end

receiverBin = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'radar_receiver', 'build', 'radar_receiver');
if isunix && exist(receiverBin, 'file')
    for mode = {'fixed', 'staggered'}
        tests(end+1, :) = {sprintf('C receiver vs MATLAB (%s PRF)', mode{1}), ...
                           @() test_c_receiver(mode{1}, receiverBin)}; %#ok<SAGROW>
    end
else
    fprintf('C receiver not built (%s): test_c_receiver skipped\n', receiverBin);
end

% A test that stops with an error counts as failed; the others still run.
passed = false(1, size(tests, 1));
for i = 1:size(tests, 1)
    try
        passed(i) = tests{i, 2}();
    catch err
        fprintf('%s: ERROR %s\n', tests{i, 1}, err.message);
    end
end

fprintf('\n===== Test summary =====\n');
for i = 1:size(tests, 1)
    if passed(i), s = 'PASS'; else, s = 'FAIL'; end
    fprintf('%-66s %s\n', tests{i, 1}, s);
end
