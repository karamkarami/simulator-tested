# Rotating search radar simulator (MATLAB)

A simple, modular MATLAB simulation of a 3 GHz rotating search radar. The
receive antenna is a column of 48 elements (selectable) with digital
beamforming (DBF) in elevation, and every PRI holds a short and a long
sub-pulse. Everything from the 120 MHz ADC samples on the 30 MHz IF down to
the PPI display is simulated. Each processing block is a short function in
`src/`. All parameters are in one file, `radar_config.m`.

The code needs no MATLAB toolbox (UDP streaming uses `udpport` or Java). It
was tested with GNU Octave 8.4. The folder `radar_receiver/` holds a C
version of the whole receiver processing (CMake project for CLion). It
processes the DDC output that the simulator streams over UDP and has an
optional Qt PPI and serial output (see `radar_receiver/README.md`).

```
 PRI:  short pulse 25 us at 0 | receive (near range) | long pulse 150 us at 200 us | receive (far range)
 per target: position from its motion -> power from the radar equation (Swerling 1)
 echoes + clutter + Tx leakage (40 dB isolation) + thermal noise, 30 MHz IF, 48 elements
 x receiver gain (60 dB) -> ADC 120 MHz, 16 bit, ENOB 12.7 -> int16
        |                       \--> optional raw data: output/raw_pulses/pulse_<n>.mat
        v
 DDC: x*cos(30 MHz) | x*(-sin(30 MHz)) -> FIR h1, /10 -> FIR h2, /6 (2 MHz) or /3 (4 MHz)
      (int16 coefficients with 13 fractional bits, time-domain sliding valid convolutions,
       every stage rounded and saturated to int16)                     -> int16 [nDdc x 48]
        |                       \--> optional DDC files: output/ddc_pulses/ddc_<n>.mat
        |                             (replay_ddc_udp.m sends them later with the PRI timing)
        |                       \--> optional UDP: int16 DDC output + pulse/scan/azimuth header
        |                             (-> radar_receiver/, the same processing in C)
        v  split in time: short-pulse samples | long-pulse samples
 DBF (21 x 48 matrix, off for 1 element)                   -> [21 beams x nSamples] per sub-pulse
 Pulse compression (valid convolution): matched filter (Taylor for LFM) or, for a phase code,
         the mismatch filter of inputs/phase_code.mat      -> [21 beams x nRange]   (range/delay axes)
 CPI buffer (circular, the pulses of one CPI)              -> [21 x nRange x 114], short and long
 MTI  fixed PRF: 3-pulse canceller | staggered: weights from file
 Doppler fixed: Hamming + FFT 128 (next power of 2 of 114) | staggered: de-interleave, Hamming, FFT 32, NCI
 CA-CFAR on every range-Doppler map: 2-D, 1-D range or 1-D Doppler window (cfg.cfar.method),
         in a selectable range interval per sub-pulse, range window shifted inwards at the edges
 Pruning (local peaks + sidelobe blanking) -> plot extraction (azimuth peak picking)
 Tracking (optional): EKF, 2-D (range, azimuth) or 3-D (+ elevation), constant velocity or
          constant acceleration, global nearest neighbour (GNN) association
 PPI display (updated while the antenna turns): plots / detections, tracks, block, console and display times
```

## Quick start

```matlab
results = run_radar_sim;            % default configuration (full 360 deg scan)

cfg = radar_config;                 % or change parameters first
cfg.prf.mode = 'staggered';
cfg.clutter.enable = false;
cfg.ant.nElements = 1;              % a single antenna: no beamforming
cfg.sim.sectors_deg = [20 45];      % only simulate this azimuth sector
cfg.io.saveRaw = false;             % recordings per PRI: raw ADC data and/or DDC output
cfg.io.saveDdc = true;
results = run_radar_sim(cfg);
replay_ddc_udp;                     % send the DDC files over UDP with the PRI timing

cfg = radar_config;                 % only the DDC files, in parallel (parfor, Parallel Computing Toolbox)
cfg.io.saveRaw = false;
cfg.ddc.outputRate_Hz = 4e6;        % DDC output at 2 MHz (default) or 4 MHz
info = generate_ddc_files(cfg, 8);  % 8 workers -> output/ddc_pulses/ddc_<n>.mat

cd tests
run_all_tests                       % SNR budget, UDP, ghost and six-target tests
```

C receiver (see `radar_receiver/README.md`):

```matlab
cfg = radar_config;  cfg.udp.enable = true;
export_receiver_config(cfg);        % configuration file for the C receiver
% start radar_receiver/build/radar_receiver --ppi, then:
results = run_radar_sim(cfg);
```

`make_input_files.m` creates the default input files in `inputs/`. The
simulator calls it automatically when a file is missing, or when a default
file was built for other values (for example after you change
`cfg.ant.nElements` or `cfg.ant.tilt_deg`). Replace a file with your own data
(same variable names) to simulate another antenna, filter set or stagger.

## Files

| File | Purpose |
|---|---|
| `radar_config.m` | all user parameters |
| `export_receiver_config.m` | writes the configuration of the C receiver (`radar_receiver/config/radar_config.bin`) |
| `run_radar_sim.m` | main loop: generation, recording, DDC/streaming, processing, display |
| `replay_ddc_udp.m` | sends recorded DDC files (`ddc_<n>.mat`) over UDP, unchanged, with the original PRI timing |
| `generate_ddc_files.m` | simulates every PRI only up to the DDC output and saves the DDC files, the PRIs in parallel (`parfor`) |
| `make_input_files.m` | default input .mat files |
| `src/setup_radar.m` | derives the PRI timeline, NCO, DDC segments, pulse compression filters (matched / mismatch), axes, CPI length, FFT length, CFAR factor ... |
| `src/mismatch_filter.m`, `ddc_pulse.m` | mismatch filter with the smallest peak sidelobe for a code or a sampled pulse; the pulse after the DDC filters |
| `src/generate_rx_signal.m`, `clutter_init.m`, `clutter_signal.m`, `adc_sample.m` | signal generation |
| `src/target_rcs.m`, `draw_target_rcs.m` | Swerling 1 target RCS (a new value every scan) |
| `src/save_pulse_file.m`, `save_ddc_file.m`, `pulse_header.m` | raw ADC data / DDC output of a PRI to a .mat file, its header |
| `src/ddc.m`, `fir_decimate_valid.m`, `ddc_decimation.m` | two-stage DDC with I/Q mixer, decimation and filter file of the 2 / 4 MHz output |
| `src/udp_sender_open.m`, `udp_make_packets.m`, `udp_send_pulse.m`, `udp_send_packets.m`, `udp_send_datagram.m`, `udp_sender_close.m`, `udp_receiver_open.m`, `udp_receive_pulse.m` | UDP streaming of the DDC output (sender and reference receiver) |
| `src/digital_beamforming.m` | DBF |
| `src/matched_filter.m`, `plot_range_profile.m` | pulse compression, range/delay/index plots |
| `src/cpi_buffer_*.m` | circular CPI buffers |
| `src/mti_filter.m`, `doppler_process.m` | clutter filter, Doppler processing |
| `src/cfar_detect.m`, `cfar_alpha.m` | CA-CFAR (2-D, 1-D range, 1-D Doppler window) |
| `src/track_setup.m`, `track_init.m`, `track_update.m`, `track_report.m` | EKF tracker (2-D / 3-D, CV / CA model, GNN association) |
| `src/prune_peaks.m`, `detection_reports.m`, `plot_extractor*.m` | pruning and plot extraction |
| `src/ppi_init.m`, `ppi_update.m`, `ppi_callbacks.m` | interactive PPI display |
| `src/detection_log_init.m`, `detection_log_add.m` | CFAR detections kept for the PPI |
| `src/cfar_pfa.m` | false alarm probability of a CFAR threshold (inverse of `cfar_alpha.m`) |
| `src/snr_budget.m`, `print_snr_budget.m` | expected SNR after every block |
| `src/six_targets_scenario.m`, `make_target.m`, `target_geometry.m`, `illumination_time.m`, `target_sectors.m` | targets |
| `tests/` | tests (see below) |
| `radar_receiver/` | the receiver processing in C (CMake/CLion, Qt PPI, serial output) |

## Input files (`inputs/`)

| File | Variables |
|---|---|
| `antenna_pattern.mat` | `elementSpacing_m` (0.66 wavelength), `az_deg`, `azGain_dB` (one-way, 4.1 deg), `elRel_deg`, `elElemGain_dB` (element pattern vs. angle from broadside, 70 deg), `txElGain_dB`, `rxChannelGain_dBi`, `txGain_dBi` |
| `ddc_filters.mat` | DDC output 2 MHz: `h1` (stage 1 at 120 MHz, 70 taps), `h2` (stage 2 at 12 MHz, 66 taps, pass +-0.5 MHz, stop above 1.5 MHz): int16 with 13 fractional bits (coefficient = value / 8192), DC gain about 1 |
| `ddc_filters_4MHz.mat` | DDC output 4 MHz (`cfg.ddc.outputRate_Hz = 4e6`): the same `h1`, `h2` with 66 taps passing +-1 MHz and stopping above 2 MHz (the Nyquist frequency of 4 MHz) |
| `dbf_matrix.mat` | `W` [21 x nElements] (beam = W * x), `beamEl_deg` |
| `stagger_prf.mat` | `priVector_s` [K x 1], `mtiCoeffs` [K x L] (row k = weights of an MTI output whose first pulse has stagger position k; one row = fixed weights) |
| `phase_code.mat` | phase codes and mismatch filters of both sub-pulses (used when `cfg.wave.<pulse>.type = 'phasecode'`): `shortCode`, `longCode` [nChips x 1] complex chip values, first chip sent first (default P4, 25 and 150 chips of 1 us); `shortMismatch`, `longMismatch`: mismatch filter as a convolution filter (output = conv(code, filter)), by default one coefficient per chip (default: 10 coefficients longer than the code, smallest peak sidelobe); optional `shortMismatchRate_Hz`, `longMismatchRate_Hz`: sample rate of the filter, when it is not one coefficient per chip (e.g. the DDC output rate) |

The default files also hold `madeFor`, the configuration values they were
built for. User files simply leave it out.

## PRI timeline and range coverage

| Time in the PRI | Event |
|---|---|
| 0 - 25 us | short pulse transmitted (`cfg.wave.short.tTx_s`); Tx leakage in the receiver |
| 25 - 200 us | short-pulse receive interval |
| 200 - 350 us | long pulse transmitted (`cfg.wave.long.tTx_s`); Tx leakage in the receiver |
| 350 - 896 us | long-pulse receive interval, up to `cfg.rx.maxRange_m` (80 km) |

The transmit times are parameters. Both sub-pulses use LFM (1 MHz) or a
phase code (1 us chips, P4 by default, from `inputs/phase_code.mat`) on the
same 30 MHz IF and are separated only in
time. The ADC records one continuous block from the start of the PRI to the
end of the long window, Tx leakage included. The DDC processes the whole
block, and its output is then split: a DDC sample belongs to a sub-pulse
when its whole filter input lies inside that sub-pulse's receive interval.
The leakage samples are therefore never processed. The short pulse covers
4.2 - 23.7 km and the long pulse 22.9 - 80.6 km (the long pulse is blind up
to c x 150 us / 2 = 22.5 km). The CFAR tests every range cell of a window
(its range window is shifted inwards at the edges, see CFAR below), so the
short pulse is tested up to the first long-pulse cell (handover at 22.9 km);
its cells beyond that are only CFAR training cells. `cfg.cfar.short` and
`cfg.cfar.long` (`minRange_m`, `maxRange_m`) restrict the tested ranges of
each sub-pulse further. The latest long-pulse start, and with it the
instrumented range, is limited by the shortest PRI; `setup_radar` stops with
an error when a receive interval is too short for the filters and the CFAR
window.

Because both sub-pulses share the IF, the short-pulse echo of a target
beyond 52.5 km also falls into the long-pulse window. The long-pulse matched
filter does not match it, so it spreads over about 10 km of range, 30 - 52 km
in front of the target, at up to about 24 dB below the target. The 2-D CFAR
adapts to it like a local rise of the noise floor, so it causes no false
plots (see `test_short_pulse_ghost`). Its only effect is lower sensitivity
there, at the same beam and Doppler, behind very strong far targets.

## Model summary

* **Geometry:** antenna 6 m above flat ground, tilted 11 deg up, turning
  clockwise at 6 rpm. Targets fly with a constant velocity (`make_target.m`).
  Azimuth is measured clockwise from North.
* **Received signal** (`generate_rx_signal.m`, `adc_sample.m`), built in
  this order for every PRI:
  1. For every target, its position from its motion and the received power
     from the radar equation (20 kW, antenna patterns, 3 dB losses, RCS of
     the current scan). The echo of both sub-pulses on the 30 MHz IF has the
     element phases of the tilted column and the Doppler inside each pulse
     (stop-and-hop between pulses).
  2. Clutter.
  3. Transmitter leakage while a pulse is transmitted: Pt - 40 dB
     isolation, in every channel.
  4. Thermal noise, k T0 F over fs/2 (NF 4 dB).
  5. The receiver gain of 60 dB.
  6. The ADC: full scale +10 dBm (a 2 Vpp sine into 50 ohm), quantization
     to 16 bits (rounding and clipping to int16), then the ADC noise for an
     ENOB of 12.7, after which the samples are rounded and clipped to int16
     again (the raw data stays signed 16-bit).

  Signals are in sqrt(W) at the receiver input. Thermal noise reaches the
  ADC at 180 LSB rms (-42 dBFS), 36 dB above the ADC noise (2.8 LSB rms).
  The leakage is +93 dBm at the ADC and saturates it while a pulse is
  transmitted. It never reaches the processing (see the PRI timeline). Very
  strong near clutter can clip the ADC as well, because there is no STC.
* **Targets:** Swerling 1. Every target gets a new exponentially
  distributed RCS (mean `rcs_m2`) on every pass of the antenna, constant
  during the dwell (`cfg.sim.swerling = 0` for constant RCS). The pass
  number changes when the antenna points away from the target
  (`target_look.m`), so a target near North does not change its RCS in the
  middle of a dwell.
* **Clutter:** constant sigma0 (-30 dB) Rayleigh cells (0.25 deg x 15 m) from
  1 km to the radar horizon, with internal motion. Scan modulation
  comes from the two-way azimuth pattern. Two Doppler spectra of the internal
  motion (`cfg.clutter.spectrum`):
  * `'gaussian'` (default): every cell moves with a random radial velocity
    of std `icmStd_mps` (0.1 m/s), a Gaussian spectrum;
  * `'exponential'`: Billingsley's model of wind-blown land clutter. Every
    cell has a stationary part with r/(1+r) of its power and a moving part
    with 1/(1+r) whose velocity is drawn from the exponential spectrum
    beta/2 exp(-beta |v|). From the wind speed w (`windSpeed_mps`, 4.5 m/s
    = 10 mph): 1/beta = 0.1 log10(w[mph]) + 0.068 m/s (beta = 5.9 s/m) and
    r = 489.8 w[mph]^-1.55 f[GHz]^-1.21 (3.6 at 3 GHz). `expBeta_spm` and
    `dcToAc` set beta and r directly. The echo is convolved with both
  transmitted IF pulses, so partly received (eclipsed) echoes are included.
* **DDC (fixed point):** the real line is `x .* cos(2*pi*30MHz*t)` and the
  imaginary line `x .* (-sin(2*pi*30MHz*t))` (NCO values with 15
  fractional bits). Each line then goes through h1 with decimation by 10
  and h2 with decimation by 6 (2 MHz output, `cfg.files.ddc`) or by 3
  (4 MHz, `cfg.ddc.outputRate_Hz = 4e6`, filters `cfg.files.ddc4MHz`). With
  4 MHz the range cells are 37.5 m instead of 75 m: twice the matched
  filter taps and range cells (the CFAR windows count cells, so they cover
  half the distance), and about twice the processing; the SNR after the
  matched filter is the same. Both are valid convolutions in the time
  domain: the filter slides over the signal and is evaluated only at the
  kept output samples (`fir_decimate_valid.m`, equal to
  `conv(x,h,'valid')` followed by `y(1:D:end)`; no FFT). It is computed in
  polyphase form: one matrix product applies blocks of D filter taps to
  blocks of D input samples, which is several times faster than a loop
  over the taps and gives exactly the same integers. The filter
  coefficients are integers with 13 fractional bits. Each stage shifts its
  accumulator right by 13 bits (`cfg.ddc.shift`, unit gain), rounds and
  saturates it to int16. The DDC output is therefore signed 16-bit, with
  the thermal noise at about 16 LSB rms per I/Q component.
* **Elements and DBF:** the antenna has a fixed 4.1 deg beam in azimuth
  and a column of `cfg.ant.nElements` elements (default 48) spaced 0.66
  wavelength in elevation. The default DBF matrix forms 21 beams equally
  spaced in sine space with a Taylor taper. With 1 element the DBF block is
  skipped, the single channel goes straight to the matched filter, and
  plots carry no elevation (NaN).
* **Pulse compression** (`cfg.wave.<pulse>.compression`): LFM uses its
  matched filter with a Taylor window. A phase code uses its matched filter
  (`'matched'`, the conjugate reversed code) or the **mismatch filter** of
  `cfg.files.phaseCode` (`'mismatch'`): a filter longer than the code that
  lowers the range sidelobes for a small SNR loss. Each coefficient is held
  for the DDC samples of one chip (by default one coefficient per chip:
  chip-matched sum, then the chip-rate filter), so one file serves the 2
  and the 4 MHz DDC output. The filter peaks earlier than the matched
  filter (`wv.mfShift` samples, found from the pulse after the DDC filters),
  which the range axis takes into account; its tail after the echo extends
  the receive window, and the range of each sub-pulse starts and ends a
  little farther inwards (default filters: 0.75 km). With the default P4
  codes and filters (`tests/test_phase_code.m`, target through the whole
  chain): peak sidelobe -21.1 -> -24.5 dB (short pulse) and -28.6 -> -33.1
  dB (long) at 2 MHz, about 4 dB better at 4 MHz, for 0.8 and 0.5 dB SNR
  loss. The 2 MHz DDC (pass band +-0.5 MHz) rounds the 1 us chips, which
  limits a filter designed on ideal chips. A filter designed at the DDC rate
  on the pulse after the DDC, over several sampling phases, does better
  (long pulse: -43 dB):

  ```matlab
  cfg = radar_config;  cfg.wave.long.type = 'phasecode';  cfg.wave.long.compression = 'mismatch';
  sys = setup_radar(cfg);
  f = rmfield(load(cfg.files.phaseCode), 'madeFor');          % a user file: no 'madeFor'
  f.longMismatch = mismatch_filter(ddc_pulse(sys, sys.wave.long, 0:0.25:0.75), 10, 3);
  f.longMismatchRate_Hz = sys.fsOut;                           % one coefficient per DDC sample
  save('inputs/my_phase_code.mat', '-struct', 'f');
  cfg.files.phaseCode = fullfile(pwd, 'inputs', 'my_phase_code.mat');
  ```
  Such a filter belongs to one DDC output rate. The C receiver gets the
  filter taps and range axes from `export_receiver_config.m` and needs no
  other setting.
* **Axes:** the delay and range of every pulse-compression output cell are
  computed in `setup_radar.m` (`wv.delay_s`, `wv.range_m`). They follow from
  the ADC sample times, the valid convolutions of the DDC and pulse
  compression filters (and the earlier peak of a mismatch filter), and the
  transmit time of the sub-pulse.
  `plot_range_profile(rp, wv, 'range'|'delay'|'index')` plots against any
  of them. Set `cfg.display.rangeProfilePulse = <pulse number>` to plot it
  during a run.
* **CPI and buffers:** the CPI is the time a target stays inside the 3 dB
  azimuth beam: 4.1 deg / 36 deg/s = 113.9 ms. Each buffer holds the pulses
  that arrive in that time, `round(CPI / PRI)`: 114 with the fixed PRF. With
  a staggered PRF the number is reduced to a multiple of the number of PRIs
  in the stagger: 116 with the default 4-PRI stagger (mean PRI 980 us). The
  buffers are circular, and MTI, Doppler and CFAR run on every new pulse
  (`cfg.proc.updateEvery`).
* **Doppler:** fixed PRF: 3-pulse canceller (112 outputs), Hamming window
  and one FFT whose length is the next power of 2 of the number of pulses
  (114 -> 128). Staggered PRF: the stagger MTI weights from the file, then
  de-interleaving into 4 rows of 28 samples. Each row is sampled every
  sum(PRI); it gets a Hamming window and a 32-point FFT (the next power of 2
  of 28), and |X|^2 is averaged over the rows.
* **CFAR:** cell-averaging CFAR with one of three windows (`cfg.cfar.method`),
  each with its own guard and training cells and Pfa (or a fixed threshold
  `threshold_dB`):

  | Method | Window (each side) | Training cells | Threshold (Pfa 1e-8) |
  |---|---|---|---|
  | `'2d'` (default) | 3 + 8 range, 3 + 4 Doppler | 296 | 12.79 dB |
  | `'range'` | 3 + 16 range cells | 32 | 13.96 dB |
  | `'doppler'` | 2 + 8 Doppler bins | 16 | 15.39 dB |

  The MTI colours the noise ([1 -2 1] gives +4 dB
  near PRF/2), so every Doppler bin is first divided by its noise gain and
  the CFAR sees a flat noise floor. The threshold factor uses the effective
  number of independent rows: 1 with a fixed PRF, 1.5 with the default
  stagger, because the de-interleaved rows are built from overlapping MTI
  outputs. The bins within one Doppler window mainlobe of zero (+-2 bins: the
  MTI notch, where the clutter residue is) are not tested.
  *Edges:* the Doppler axis is circular (the spectrum is periodic), so the
  Doppler window is always complete. Range is not periodic: wrapping it
  around would mix the near-range clutter with far-range noise. Instead,
  at the first and last range cells of a sub-pulse window the range window is
  shifted inwards: it keeps its length and so its training cells, the cell
  under test moves away from its middle, and its guard cells are cut at the
  edge (so it has a few more training cells there; the noise is their mean).
  Before, the outermost 11 cells of each window (0.8 km at each end) were not
  tested at all. The noise estimate costs the same: the sums of the shifted
  windows are sums already computed for the nearest interior cell.
  Only the cells in the tested range interval (and the cells their windows
  need) are processed, so a smaller interval is also faster.
  The 2-D window adapts best to clutter residue (its Doppler extent sees the
  residue around the notch). A 1-D range window at the start of the
  short-pulse window, where the near-range clutter falls off steeply, only
  sees weaker clutter farther out: in the six-target test the `'range'`
  method gave one false plot at 4.2 km (clutter residue at -2.3 m/s) that
  the `'2d'` and `'doppler'` methods do not give. Raise
  `cfg.cfar.short.minRange_m` (e.g. to 5 km) or the threshold there.
  The 1-D Doppler window has no range extent. The long pulse of a strong
  target inside its blind range (22.5 km) is still received in part
  (eclipsed) at the start of the long-pulse window, where its compressed
  response spreads over several km in the target's Doppler bin: with the
  `'doppler'` method the near target T5 (7.7 km, SNR 62 dB) gave 4 false
  long-pulse plots at 23 - 26 km (25 with phase codes and mismatch
  filters), which the `'2d'` method does not give.
* **Tracking (optional, `cfg.track`):** an extended Kalman filter on the
  plots. State relative to the antenna (East, North, up), `cfg.track.model`:
  `'cv'` (default) position and velocity, nearly constant velocity with
  white-acceleration process noise (`sigmaAccel_mps2` = 5 m/s^2 horizontal,
  `sigmaAccelZ_mps2` = 1 m/s^2 vertical); `'ca'` position, velocity and
  acceleration, nearly constant acceleration with white-jerk process noise
  (`sigmaJerk_mps3` = 0.5, `sigmaJerkZ_mps3` = 0.1 m/s^3) and the initial
  acceleration std `maxAccel_mps2`/sqrt(3) (5 m/s^2: with three unknowns per
  axis, three plots always fit a track, so a wide start lets false plots
  confirm tracks). CA follows accelerating targets (4 m/s^2 in
  `test_tracker`: 136 m instead of 394 m position error) and is less precise
  than CV for targets flying straight. `cfg.track.dim`: `'3d'` measures range, azimuth and elevation
  (needs the elevation of the DBF), `'2d'` ground range and azimuth in the
  horizontal plane (slant range when the elevation is unknown: 1 element),
  `'auto'` (default) takes 3-D with DBF. The measurement h(x) is linearised
  at the predicted state. Plot errors: 50 m, 0.3 deg, 0.5 deg.
  Association: global nearest neighbour (GNN). The plots are collected
  until the antenna has turned `gnnSector_deg` (5 deg) without a new plot,
  at most `gnnMaxSector_deg` (30 deg) after the first one. Then every
  plot-track pair inside the chi-square gate (`gateProb` 0.999) costs
  d^2 + ln|S|, and the assignment (each track at most one plot of the batch,
  each plot at most one track) has the most pairs and, among those, the
  smallest total cost (Hungarian algorithm). Unlike a nearest-neighbour
  association plot by plot, two plots never update the same track, and
  close targets keep their own tracks (`test_tracker`: two targets 200 m
  apart). A plot without a track starts a tentative track (zero velocity
  and acceleration, velocity std `maxSpeed_mps`/sqrt(3)). A track is confirmed after
  3 updates and deleted after more than 3 missed scans (tentative: 1). The
  tracks are in `results.tracker` and on the PPI.
* **Pruning and plots:** only local maxima in beam x range x Doppler are
  kept. Range, elevation and Doppler sidelobes more than 20 dB below a
  stronger peak are removed. Detections of the same target in successive
  (sliding) CPIs are grouped, and the strongest one becomes a plot. A plot
  needs detections over at least 2 deg of antenna rotation, half a beamwidth
  (`cfg.extract.minExtent_deg`). Noise false alarms last about 1 deg, while
  a target is detected over about 6-8 deg. Elevation is interpolated between adjacent
  beams. The plot SNR is measured against the noise floor of the beam. The
  radial velocity is ambiguous (+-25 m/s with the fixed PRF, +-6.4 m/s with
  the default stagger).

## Interactive PPI

The PPI (`ppi_init.m`) is interactive, also while the simulation runs:

| Action | Effect |
|---|---|
| mouse wheel | zoom in/out around the cursor (up to x500) |
| left drag | pan |
| double click, key `r` | reset the view (whole instrumented range) |
| hover over a plot, a detection, a track or a true target | data tip: range, azimuth, elevation, radial velocity, SNR, time, sub-pulse, CPI hits, beam/cell/bin (track: number, state, range, azimuth, height, climb, speed, heading, updates; true target: range, azimuth, elevation, height, speed, RCS) |
| click on a plot or detection | pin its data tip (click again to remove it) |
| keys `+` `-`, arrows | zoom, pan |
| key `c` | remove the pinned data tips |
| key `d` | next display mode (plots / detections / both) |
| key space | pause/resume the display (the simulation continues) |

The panel on the right has:
* **Show**: *Plots (after extractor)*, the plots of the plot extractor as
  before; *Detections (before extractor)*, the pruned CFAR detections of
  every CPI update, i.e. what goes into the plot extractor; or both. Every
  detection is drawn as an arc over the azimuth step between CPI updates
  (rotation rate x PRI x `cfg.proc.updateEvery`, 0.036 deg by default) at its
  ground range, so a target appears as an arc about one beamwidth long (one
  per detected range cell and Doppler bin), and a noise false alarm as a
  short dash. Start value: `cfg.display.showMode`. The newest
  `cfg.display.maxDetections` (100 000) detections are kept.
* show/hide of the short- and long-pulse data, the true targets, range
  rings, plot labels (range and velocity next to each plot), the antenna and
  the **tracks**; a minimum SNR and a choice of what is shown (all, last
  scan, last 2 scans). The filters apply to plots and detections.
  A track is a hollow triangle at its last update, so its plot (a filled
  circle) lies inside it: magenta when confirmed (with its number and a
  velocity leader of `cfg.display.trackLeader_s` = 30 s of flight), gray
  while tentative. Start value of the checkbox: `cfg.display.showTracks`.
* **CFAR threshold [dB]**: the threshold factor 10 log10(alpha) of the
  selected CFAR method (`cfg.cfar.method`). *Apply* hands the value to the
  running simulation at its next display update (`ppi_update` returns `sys`
  with the new `sys.cfar.alpha`), from then on every CFAR uses it. *Default*
  goes back to the value of the method's Pfa. The line below shows the
  method, the threshold in use and its false alarm probability (`cfar_pfa.m`,
  the inverse of `cfar_alpha.m`), e.g. `2-D: 12.79 dB (Pfa 1.0e-08)`.
* Reset view, Clear display, Remove data tips, Pause, and **Save screenshot
  (PNG)**, which saves the window as it is now (file dialog, default name
  `ppi_<date>_<time>.png`).
* **Mean time per PRI [us]** (bottom right, above the help): the mean
  time per PRI of every block since the start: DBF, matched filter, MTI
  (clutter filter, with the CPI buffer), Doppler (FFT and non-coherent
  integration), CFAR of the selected method with the pruning, plot
  extractor and tracking, then the console output (progress lines) and the
  display (the PPI updates, every `cfg.display.updateEvery` PRIs, averaged
  over all PRIs), and their total (`results.timing_us`).

The cursor position (range, azimuth, east/north) is shown under the buttons. In MATLAB the zoom, pan and data tip tools of the axes
toolbar work as well; their data tips also list the plot values.
`ppi_callbacks(action, fig, ...)` does the same from a script (see
`tests/test_ppi_interactive.m`).

## Raw data, DDC files and UDP format

Two recordings can be switched on and off separately in `radar_config.m`:

| Option | Files (one per PRI) | Content |
|---|---|---|
| `cfg.io.saveRaw` (default `true`) | `cfg.io.rawFolder/pulse_<n>.mat` | `adcData`: raw ADC samples, int16 [nSamples x nElements] |
| `cfg.io.saveDdc` (default `false`) | `cfg.io.ddcFolder/ddc_<n>.mat` | `I`, `Q`: DDC output, int16 [nDdc x nElements], both sub-pulses |

Both files also hold `hdr`: pulse number, scan number, time, antenna
azimuth, PRI, stagger index, fs, IF, carrier, record start, pulse timing,
DDC rate and the DDC samples of each sub-pulse. At 80 km instrumented range
one PRI is about 9.5 MB of raw data (a full scan of 10 000 PRIs is about
95 GB) and 342 kB of DDC output (3.4 GB per scan; twice that with the
4 MHz DDC output). Use `cfg.sim.sectors_deg` or a smaller
`cfg.rx.maxRange_m` for long runs.

### Generating the DDC files in parallel

`generate_ddc_files(cfg, nWorkers)` runs only the first three steps of
every PRI (signal generation, ADC, DDC) and saves the same `ddc_<n>.mat`
files (and, with `cfg.io.saveRaw`, the raw data), with no processing after
the DDC. The PRIs are independent of each other, so they run in a `parfor`
loop over the workers of the Parallel Computing Toolbox (`nWorkers`
workers; `0` = this MATLAB session; without the toolbox and in Octave the
loop runs serially). The target RCS of every scan and the clutter map are
drawn first, as in `run_radar_sim`; the noise of every pulse comes from its
own random stream (substream = pulse number), so the files do not depend on
the number of workers or on the order of the pulses. Without noise the
files are bit-identical to those of `run_radar_sim` (`test_generate_ddc`).

```matlab
cfg = radar_config;  cfg.io.saveRaw = false;  cfg.sim.nScans = 3;
info = generate_ddc_files(cfg, 8);                    % -> cfg.io.ddcFolder
radar_receiver/build/radar_stream_send output/ddc_pulses --port 50000   % (shell) to the C receiver
```

### Replaying the DDC files over UDP

`replay_ddc_udp.m` reads the DDC files and sends them to a receiver
without any processing: the stored samples and header go out in the same
datagrams as the live stream below, so the receiver (for example
`radar_receiver/`) cannot tell a replay from the simulator. The PRIs are
sent with the original timing:

```matlab
cfg = radar_config;
cfg.io.saveDdc = true;  cfg.io.saveRaw = false;
cfg.sim.sectors_deg = [20 45];
run_radar_sim(cfg);                                   % -> output/ddc_pulses/ddc_<n>.mat

stats = replay_ddc_udp;                               % to cfg.udp.host:cfg.udp.port
stats = replay_ddc_udp('output/ddc_pulses', 'host', '192.168.1.20', 'port', 50000, ...
                       'payloadBytes', 65400, 'spread', 0.8);
```

* The send time of every PRI comes from the PRI start times in the headers
  (`hdr.time_s`), on one clock, so no error accumulates. The loop sleeps
  until 20 ms before a send time and then busy-waits; on an idle PC a PRI
  leaves within a few microseconds of its time.
* Files are read and split into datagrams ahead of time: up to
  `'preloadMB'` (2048 MB) before the start, the rest only while the next
  datagram is not yet due, so reading never delays a PRI.
* `'spread'` sends the datagrams of a PRI evenly over this fraction of the
  PRI instead of as one burst (smoother data rate for small receive buffers).
* A PRI that still goes out late (the operating system stalled the
  sender) does not cause a burst: the following PRIs keep the PRI spacing.
  `'catchUp', true` instead sends them at once until the original schedule
  is reached again. `'keepGaps', false` shortens a skipped sector to one PRI;
  `'speed'` replays faster or slower.
* `stats.sent_s` (send time of every PRI) and `stats.late_s` (its delay)
  show whether the timing was kept; a summary is printed at the end.
* Data rate: one PRI of the default radar is 342 kB, 2.7 Gbit/s at a 1 ms
  PRI (localhost or 10 Gbit Ethernet). MATLAB needs tens of microseconds
  per datagram, so use few datagrams per PRI: `'payloadBytes', 65400` gives 6
  datagrams per PRI (the network stack fragments them; the C receiver accepts
  datagrams up to 64 kB) instead of 245 with 1400 bytes.

For exact real-time timing use the C sender `radar_stream_send` of
`radar_receiver/` (see its README). It reads the same folder (files saved
with `-v6`, or `-v7` when built with zlib), sends the same datagrams, and
keeps the PRI timing within microseconds at 2.7 Gbit/s:

```sh
radar_receiver/build/radar_stream_send output/ddc_pulses --host 192.168.1.20 --port 50000
```

### UDP stream

With `cfg.udp.enable = true`, the **DDC output** of every PRI (signed
16-bit I/Q, [nDdc x nElements], both sub-pulses) is sent to
`cfg.udp.host:port`. Each datagram has a 64-byte little-endian header
followed by int16 (I, Q) pairs, channel after channel. The default payload
of 1400 bytes keeps each datagram within a standard 1500-byte Ethernet MTU.
One PRI of the default radar is 48 x 1781 x 4 bytes = 342 kB in 245
datagrams.

| Bytes | Type | Field |
|---|---|---|
| 0-3 | uint32 | magic 0x52444443 ('RDDC') |
| 4-7, 8-11 | uint32 | pulse number, scan number |
| 12-15 | single | antenna azimuth [deg] |
| 16-17, 18-19 | uint16 | packet index (0-based), packet count |
| 20-21, 22-23 | uint16 | channels, stagger index |
| 24-27, 28-31 | uint32 | samples per channel, first sample in this packet (0-based) |
| 32-39 | double | PRI start time [s] |
| 40-43, 44-47 | single | DDC sample rate [Hz], time of the first DDC sample in the PRI [s] |
| 48-63 | uint32 | short pulse first/last DDC sample, long pulse first/last DDC sample (1-based) |

`src/udp_receive_pulse.m` is a reference receiver that reassembles one PRI.

## SNR budget

`print_snr_budget` (called at the start of every run) prints the expected
SNR of every target after the ADC, DDC, DBF, matched filter, MTI and Doppler
blocks. Each target is evaluated at the moment the beam crosses it, with its
mean RCS. The table also gives the SNR needed for Pd = 0.9 and the margin:
for a Swerling 1 target (exact for one look, Shnidman's approximation for
more) or, with `cfg.sim.swerling = 0`, for a steady target (Albersheim).
Where a closed form would be approximate, the budget runs the real
processing functions on a noise-free target:
* the matched-filter loss: the pulse, with its Doppler shift, goes through
  the DDC filters and the pulse compression filter. This covers the Taylor
  window, the loss of a mismatch filter, the DDC filtering of phase-code
  chips, and the Doppler mismatch.
* the slow-time gain: the MTI and Doppler functions run on the target and
  on white noise. This covers the MTI response, the Hamming loss, Doppler
  straddle and the azimuth beam shape over the CPI.

`tests/test_snr_budget.m` measures the SNR after every block of the
simulated chain (signal-only and noise-only runs) and compares it with the
budget.

## Tests

| Test | What it checks |
|---|---|
| `test_snr_budget(mode, pulse, type, nElements, ddcRate, compression)` | measured vs. expected SNR after each block (LFM or P4 with its matched or mismatch filter, 1 or 48 elements, 2 or 4 MHz DDC output); the noise power of every block is the mean over the pulses of the CPI (after the Doppler processing: over the Doppler bins outside the MTI notch) |
| `test_phase_code` | phase codes from the file, matched and mismatch filters: a target through the whole chain lands on its range cell with both (2 and 4 MHz), the mismatch filter has lower sidelobes and less than 1.5 dB loss; a user filter at the DDC rate; `run_radar_sim` with mismatch filters gives the right plot; the C core test passes with the exported mismatch configuration |
| `test_udp_stream` | the DDC output of a PRI sent over UDP to localhost is received bit-exact, with its header |
| `test_ddc_replay` | a short run saves only the DDC files (no raw data); `replay_ddc_udp` sends them to localhost bit-exact, with the PRI timing (median delay below 10 % of the PRI); when built, the C sender `radar_stream_send` sends them (saved with `-v6` and `-v7`) bit-exact as well |
| `test_generate_ddc` | `generate_ddc_files`: without noise its DDC files equal those of `run_radar_sim` bit by bit (2 and 4 MHz, header included); with noise a pulse gives the same file whatever other pulses are generated (order and number of workers do not matter) |
| `test_short_pulse_ghost` | a strong target at 65 km is detected, and its short-pulse echo in the long window gives no ghost plot |
| `test_six_targets(mode, clutterOn, keepRaw)` | 6 targets (0.5 - 10 m^2 mean RCS, Swerling 1, 8 - 48 km, up to 290 m/s, around the radar) are detected with the right range, azimuth, elevation and ambiguous velocity, with at most 2 false plots; one raw file per PRI, named by the pulse number |
| `test_ppi_interactive` | PPI interaction: data tips, zoom around the cursor, pan, keys, filters, pinned tips, pause, clear, PNG; detection display (modes, arcs, data tips, filters); CFAR threshold field (Apply, hand-over through `ppi_update`, invalid input, Default, Pfa); tracks (triangles, data tip, checkbox); block processing times |
| `test_setup` | configuration edge cases without a simulation: the simulated sectors (`[0 360]` = whole circle, across North, several sectors), settings that cannot work stop with a message naming them (MTI notch over all Doppler bins, CPI too short, unknown CFAR method, DDC rate, PRF mode), a short pulse that hands over within its first CFAR window keeps one whole window (and the C receiver accepts the exported configuration) |
| `test_cfar` | the three CFAR methods on exponential noise: false alarm rate 1e-3 in the interior and at the range edges (shifted window), targets in the first and last range cell detected, only the tested range interval is tested, the noise estimate equals the direct sums over the (shifted) window, also when the whole tested interval lies within the edge cells |
| `test_tracker` | EKF tracker on synthetic plots of the six targets over 5 scans with false plots, 3-D and 2-D, CV and CA: one confirmed track per target, close to it, no false track; GNN with two targets 200 m apart (nearest neighbour would merge them); CA against CV on a target accelerating at 4 m/s^2; when built, the C tracker gives the same track updates |
| `test_c_receiver(mode, ..., cfarMethod)` | the six-target scenario streamed over UDP to the C receiver (`radar_receiver/`): its plots and track updates must match MATLAB; the recorded stream is replayed to check repeatability and speed |

Latest results (Octave 8.4 on a 4-vCPU VM, every test of `run_all_tests`
passed): the SNR budget checks agree within 0.45 dB (ADC, DDC and DBF within
0.25 dB). The six-target tests find all 6 targets with no false plot, for
both the fixed and the staggered PRF. In `test_c_receiver` the C receiver
gives the same 6 plots as MATLAB (same cells, range, azimuth, elevation,
velocity, SNR and number of detections) and the same track updates for both
PRF modes. It also matches MATLAB exactly (two or three of the targets each)
with unusual configurations: staggered PRF with the 4 MHz DDC output, the
1-D range CFAR and the CA tracker; phase codes with mismatch filters, the
1-D Doppler CFAR and 2-D tracking; a single element with Doppler/CFAR every
4 pulses and the exponential clutter spectrum; restricted CFAR range
intervals; the default 80 km receive window. A six-target run takes about
0.85 s per PRI (2400 PRIs in 35 minutes, raw files included). The polyphase
DDC halved the run time (one target with clutter, 397 PRIs: 1.65 s per PRI
before, 0.82 s after, same plot). After the CFAR changes (methods, shifted
range window, no extra CFAR cells in the receive window): the 2-D CFAR still
finds all 6 targets with no false plot (fixed PRF), and the C receiver gives
exactly the same plots and track updates as MATLAB with the 2-D (six
targets, fixed PRF; two targets, staggered PRF), the 1-D range and the 1-D
Doppler CFAR (two targets). With 3 scans of two targets both tracks are
confirmed: speed within 4 m/s, heading within 0.3 deg, height within 25 m of
the truth.

`tests/run_all_tests.m` runs all of them for the fixed and the staggered
PRF (about 3 hours in Octave on a 4-vCPU VM). The six-target test simulates
only the azimuth sectors around the targets, about 2400 PRIs. Its raw files
are about 18 GB per PRF mode; `run_all_tests` deletes them after checking
them. Octave without a display (`octave-cli` has no graphics toolkit) runs
the simulator without the PPI and skips the interactive PPI test; for the
PPI use a virtual display: `xvfb-run octave --no-gui --eval run_all_tests`.

## Assumptions and limits

* Flat earth for target geometry. The clutter extends to the 4/3-earth
  radar horizon, and all clutter arrives from one elevation (0 deg).
* The transmitter uses the wide element pattern in elevation. Every
  receive element sees the full azimuth pattern.
* Swerling 1 targets (constant RCS within a scan), no multipath, no echoes
  from the previous PRI.
* LFM range-Doppler coupling is not corrected (about 130 m at 290 m/s with the
  150 us pulse), because the measured velocity is ambiguous.
* Signal generation uses single precision. That is enough for a 16-bit ADC
  and about twice as fast. The DDC computes with integers (exact in
  double). The processing after the DDC uses single precision, as the C
  receiver does.
