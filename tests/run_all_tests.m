% RUN_ALL_TESTS  Run every test of the radar simulator and print a summary.
%   The six-target tests take about 35 minutes each in Octave 8.4 on a
%   4-vCPU VM. Each writes about 18 GB of raw data (output/test_<mode>/raw_pulses),
%   checks it and deletes it again (call test_six_targets(mode, true) to keep
%   the files). When the C receiver is built (radar_receiver_<os>/build/radar_receiver,
%   Linux or macOS), test_c_receiver compares it with the MATLAB processing
%   (about 20 minutes per run: fixed and staggered PRF, 1-D range and Doppler
%   CFAR, the GO, SO and OS detectors and the PMBM tracker), and test_multi_scan sends three scans through it (about 40 minutes).
%   test_clutter takes about 5 minutes, test_tracker_pmbm about 3 minutes, each
%   test_scenarios run (3 scans, 48 elements) 2 - 3 hours (1 element: about
%   20 minutes). All tests together take about 12 hours.
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
tests(end+1, :) = {'CFAR detectors (CA, GO, SO, OS) x windows (2-D, range, Doppler), edges, range interval', @test_cfar};
tests(end+1, :) = {'EKF tracker 3-D and 2-D (and the C tracker when built)', @test_tracker};
tests(end+1, :) = {'PMBM tracker with IMM-UKF 3-D and 2-D: maneuvers, acceleration, clutter, close targets', @test_tracker_pmbm};
tests(end+1, :) = {'TOMHT tracker (k-best, N-scan) with IMM-UKF: the PMBM cases (and the C tracker)', @test_tracker_tomht};
tests(end+1, :) = {'Options: sub-pulses on/off and carriers, CPI mode, MTI range, MF/Doppler recordings, video floor', ...
                   @test_new_options};
tests(end+1, :) = {'Peak interpolation of the detections (Doppler, elevation)', @test_peak_interpolation};
tests(end+1, :) = {'Clutter filter option: canceller or FIR filter from a .mat file', @test_mti_fir};
tests(end+1, :) = {'Interactive PPI (data tips, zoom, pan, filters)', @test_ppi_interactive};
tests(end+1, :) = {'No ghost of the short-pulse echo (same IF)', @test_short_pulse_ghost};
tests(end+1, :) = {'170 MHz challenges (overload, eclipsed echoes, blind speeds, resolution, load)', @test_vhf_challenges};
tests(end+1, :) = {'Clutter: CNR, terrain types (mountain, hill, desert, sea), shadowing, overload, one CPI per terrain', @test_clutter};
for mode = {'fixed', 'staggered'}
    tests(end+1, :) = {sprintf('Six targets with clutter (%s PRF)', mode{1}), ...
                       @() test_six_targets(mode{1}, true, false)}; %#ok<SAGROW>
end
tests(end+1, :) = {'Thirteen challenging targets in the terrain clutter, 3 scans: plots and tracks (fixed PRF)', ...
                   @() test_scenarios('challenge')};
tests(end+1, :) = {'Thirteen challenging targets, 3 scans, PMBM tracker (IMM-UKF), 48 elements (3-D)', ...
                   @() test_scenarios('challenge', 'method', 'pmbm')};
tests(end+1, :) = {'Thirteen challenging targets, 3 scans, PMBM tracker (IMM-UKF), 1 element (2-D)', ...
                   @() test_scenarios('challenge', 'method', 'pmbm', 'nElements', 1)};
tests(end+1, :) = {'Ten random targets in the terrain clutter, 3 scans: plots and tracks (staggered PRF)', ...
                   @() test_scenarios('random', 'prf', 'staggered')};

receiverBin = fullfile(receiver_dir(), 'build', 'radar_receiver');
if isunix && exist(receiverBin, 'file')
    for mode = {'fixed', 'staggered'}
        tests(end+1, :) = {sprintf('C receiver vs MATLAB (%s PRF)', mode{1}), ...
                           @() test_c_receiver(mode{1}, receiverBin)}; %#ok<SAGROW>
    end
    for method = {'range', 'doppler'}
        tests(end+1, :) = {sprintf('C receiver vs MATLAB (fixed PRF, %s CFAR)', method{1}), ...
                           @() test_c_receiver('fixed', receiverBin, [], [], false, method{1})}; %#ok<SAGROW>
    end
    % other CFAR detectors and the PMBM tracker in both programs
    combos = {'fixed', '2d', 'os', 'pmbm'; 'fixed', 'range', 'go', 'pmbm'; 'fixed', 'doppler', 'so', 'ekf'; ...
              'staggered', '2d', 'ca', 'pmbm'; 'fixed', '2d', 'ca', 'tomht'};
    for i = 1:size(combos, 1)
        c = combos(i, :);
        tests(end+1, :) = {sprintf('C receiver vs MATLAB (%s PRF, CFAR %s %s, %s tracker)', c{1}, upper(c{3}), c{2}, ...
                           upper(c{4})), @() test_c_receiver(c{1}, receiverBin, [], [], false, c{2}, c{3}, c{4})}; %#ok<SAGROW>
    end
    tests(end+1, :) = {'Three scans through the C receiver: tracks on the targets only (fixed PRF)', ...
                       @() test_multi_scan('fixed', 3)};
    tests(end+1, :) = {'C DDC file generator radar_generate_ddc against generate_ddc_files (and the C receiver)', ...
                       @test_generate_ddc_c};
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
