# Radar receiver in C (CLion / CMake project, Linux)

The receiver is kept as two CLion projects with the same sources:
`radar_receiver_linux` (this one: Ubuntu with GCC and the Qt of the
distribution, optionally the Ubuntu real-time kernel, see "Linux and real
time") and `../radar_receiver_windows`. Their `src`, `tests`, `tools` and
`config` folders are identical (test `same_sources`), so every change goes
into both; only `CMakeLists.txt` (the default Qt location) and the README
differ.

This program does all the receiver processing of the MATLAB simulator in C.
Its input is the signed 16-bit DDC output of every PRI, which the simulator
(or real hardware using the same packet format) sends over Ethernet/UDP. It
writes the plots to the console and a CSV file. Optionally it also shows them
on a Qt PPI and sends the CFAR detections (or the plots) over a serial port
at 921600 baud.

```
 UDP datagrams (int16 I/Q + header)  -- or a recorded file (--replay)
        |  input thread: reassembly of each PRI (udp_input.c) -> queue of PRIs (input_thread.c)
        v  split in time: short-pulse samples | long-pulse samples
 DBF (21 x 48 matrix, off for 1 element)
 Pulse compression: matched filter, or the mismatch filter of a phase code
         (valid convolution; FFT fast convolution when it is faster)
 CPI ring buffer (the pulses of one CPI) -> MTI ring (MTI output of each pulse computed once)
 Doppler: fixed PRF Hamming + FFT 128 | staggered: de-interleave, Hamming, FFT 32, NCI;
          division by the noise gain of every bin
 CFAR (cfar.c): CA, GO, SO or OS detector with the 2-D, 1-D range or 1-D Doppler window
   (tested range interval, shifted range window at the edges)
   -> pruning (local peaks, sidelobe blanking) -> reports -> plot extraction
   -> tracking (optional, tracker.c): EKF + GNN (tracker_ekf.c) or PMBM + IMM-UKF (tracker_pmbm.c)
        |                                   |                         |                |
        |                                   +-> serial (detections)   +-> CSV, console, Qt PPI,
        |                                                                 serial (plots) +-> CSV, Qt PPI
```

The algorithms are the same as in the MATLAB simulator (`../src/*.m`, see
`../README.md`). The C receiver reads every number that MATLAB derives
(DBF matrix, pulse compression filters, range axes, CPI length, MTI
weights, Doppler window and noise gains, CFAR detector, window, threshold
factor and tested range cells, extraction and tracker settings) from a
configuration file that `export_receiver_config.m` writes. Both programs
therefore always use the same values.

**Phase codes and mismatch filters** (`cfg.wave.<pulse>.type =
'phasecode'`, `cfg.wave.<pulse>.compression = 'matched' | 'mismatch'`,
codes and filters in `../inputs/phase_code.mat`): the exported taps are the
matched or the mismatch filter at the DDC rate, and the exported range axis
already holds the earlier peak of the longer mismatch filter. The filter may
have any length: the pulse compression input of each sub-pulse
(`<pulse>_ddcFirst`, `<pulse>_ddcCount`) then starts before and ends after
its receive interval (`<pulse>_rxFirst`, `<pulse>_rxCount`), and those samples
are zero, so the range cells are the same as with the matched filter. The
receiver filters with whatever taps the file has (direct or FFT, as for
LFM), so it needs no setting of its own; it prints the waveform and the
filter of each sub-pulse at the start, e.g. `long pulse: phase code of 150
chips, mismatch filter of 320 taps (22.90 ... 80.61 km)`. The file also holds the pulse
itself at the DDC rate (`<pulse>_pulse`), which `test_core` uses as the
echo. `config/radar_config_phasecode.bin` is such a configuration (both
sub-pulses with mismatch filters) for the `core_phasecode` test.

## Build

Requirements: a C11 compiler (GCC, Clang, or MSVC 2019 16.8 or newer), CMake 3.21 or newer,
OpenMP (optional, for multi-threading), Qt 5.12+ or Qt 6 Widgets
(optional, for the PPI) and zlib (optional: compressed `-v7` MAT-files, the
DDC files of `radar_stream_send` and the input files of `radar_export_config`,
are then inflated by zlib; without it a built-in inflate reads them).

**CLion:** *File > Open* the `radar_receiver_linux` folder; CLion loads
`CMakeLists.txt`. Choose the **Release** profile (*Settings > Build,
Execution, Deployment > CMake*, add a Release profile). Without it the code
is not optimised and runs about 10 times slower. Run configurations are
created for `radar_receiver` and the tests. Set the program arguments in the
run configuration, for example `--config ../config/radar_config.bin --ppi`.
If CMake does not find Qt, add `-DCMAKE_PREFIX_PATH=<Qt>/lib/cmake` (for
example `~/Qt/6.7.2/gcc_64/lib/cmake` of the Qt installer) to the CMake
options.

**Command line:**

```sh
cmake -S radar_receiver_linux -B radar_receiver_linux/build -DCMAKE_BUILD_TYPE=Release
cmake --build radar_receiver_linux/build -j
cd radar_receiver_linux/build && ctest --output-on-failure
```

CMake options:

| Option | Default | Meaning |
|---|---|---|
| `RADAR_WITH_QT` | ON | build the PPI (the build continues without Qt when it is not found) |
| `RADAR_OPENMP` | ON | OpenMP threads |
| `RADAR_NATIVE` | ON | optimise for the CPU of the build machine (`-march=native`, MSVC `/arch:AVX2`); switch it off when the program must run on other PCs |
| `RADAR_BUILD_TESTS` | ON | unit tests (`ctest`) |

On Ubuntu the packages of "Linux and real time" are all that is needed;
CMake finds the Qt of the distribution. The serial port is opened as
`/dev/ttyUSB0` and so on (user in the group `dialout`).

## Use with the MATLAB simulator

```matlab
cfg = radar_config;                       % the same cfg for both programs
cfg.udp.enable = true;                    % stream the DDC output
export_receiver_config(cfg);              % -> radar_receiver_linux/config/radar_config.bin
```

```sh
radar_receiver_linux/build/radar_receiver --config radar_receiver_linux/config/radar_config.bin --ppi
```

Then run `run_radar_sim(cfg)` in MATLAB. Export the configuration again
whenever `cfg` changes: the receiver checks that the PRI layout (number of
channels and DDC samples) matches and ignores PRIs that do not. The file
has a format version (2); the receiver refuses files of another version and
needs every value that `export_receiver_config.m` writes.

**Pruning:** as in MATLAB (`prune_peaks.m`): local peaks in beam x range x
Doppler; range, elevation and Doppler sidelobes more than
`sidelobeMargin_dB` below a stronger peak are dropped. The local peak test
is done inside the CFAR: a cell is detected when it is above the threshold
*and* the maximum of its 3 x 3 x 3 neighbourhood (the same cells as the
threshold followed by the peak test). In a row of the map with many cells
above the threshold the peak test is one vectorised pass with the maxima of
the neighbouring beams and Doppler bins, in a row with few each cell is
checked, so its time no longer grows with the number of threshold
crossings. The neighbourhood is 3 x 3 (x 3 with DBF) also for the 1-D range
and 1-D Doppler CFAR: a target spreads over neighbouring range cells and
Doppler bins whatever the CFAR window, and a peak test along one dimension
only would keep 3 - 5 cells of every target (one per Doppler bin of its
mainlobe with the range CFAR). Only the peaks more than the margin above the
weakest peak can blank another one; only those are sorted by power, and
each peak is compared with the ones that much stronger (binary search).
When those comparisons would cost more than about 20 passes over the map (a
low CFAR threshold with many strong peaks, e.g. clutter residue),
`blank_dense()` takes the largest peak power of each sidelobe neighbourhood
over the whole map instead (running maxima, van Herk / Gil-Werman along
range): the same result, in a time bounded by the map size.
`tests/test_c_receiver.m` does all of this automatically and compares the
plots of both programs.

## Configuration without MATLAB (`radar_export_config`)

`radar_export_config` (`tools/export_config.c`) writes the receiver
configuration from `radar_config.m` itself, without MATLAB: the same file,
the same records and values as `export_receiver_config.m` (within rounding,
1e-14; the test `export_config` and the comparison below). In CLion it is a
run configuration of its own, **radar_export_config** (CLion makes one for
every program of `CMakeLists.txt`): choose it and press *Run*. Without
arguments it reads `../radar_config.m` of this project and writes
`config/radar_config.bin` of this project, whatever the working directory,
then start `radar_receiver` as usual. Options (*Program arguments* of the run
configuration, or the command line):

| Option | Meaning |
|---|---|
| `--config-m FILE` | the MATLAB configuration (default `../radar_config.m` of the project) |
| `--out FILE` | the receiver configuration (default `config/radar_config.bin` of the project) |
| `--set STATEMENT` | a MATLAB assignment run after the file, as a change of `cfg` before `export_receiver_config(cfg)` in MATLAB; repeatable, e.g. `--set "cfg.prf.mode = 'staggered';" --set "cfg.cfar.method = 'range';"` |
| `--quiet` | print only errors |

```sh
radar_receiver_linux/build/radar_export_config --set "cfg.mti.type = 'fir';" --set "cfg.report.dopplerInterp = true;"
```

How it works: `tools/mfile.c` evaluates `radar_config.m` as MATLAB would, for
what a configuration file uses: the function line, assignments to `cfg`
fields and variables (`root = fileparts(mfilename('fullpath'))`), comments
and `...` lines, numbers, strings, matrices (`[0 40; 1000 40]`, `[1 -2 1]`),
cell arrays, `struct(...)`, `fullfile`, `Inf`, arithmetic, comparisons and
common functions (`sqrt`, `round`, `sind`, `min`, `max`, `linspace`, ...).
A function it does not know (`target_scenario(cfg)`: the simulated targets,
which the receiver does not need) is left unevaluated; it is an error only
if the export needs that value. `if`/`for` blocks and indexed assignments
(`x(2) = 1`) are refused with the line. Then `tools/export_config.c` does
what `setup_radar.m` does for the export: it reads the input files of
`cfg.files` (antenna pattern, DDC filters of the selected output rate, DBF
matrix, stagger PRIs and MTI weights, phase codes and mismatch filters, the
FIR clutter filter of `cfg.mti.type = 'fir'`), and derives the 3 dB
azimuth beamwidth, the record timing and DDC sample times, the receive
intervals, the pulse compression filters (Taylor-windowed LFM, phase code
with the matched or the mismatch filter and its peak shift from the pulse
after the DDC filters), range axes and the handover, the CFAR-tested cells,
the CPI, the Doppler window, the noise gain and the looks of every Doppler
bin (impulses through the MTI, the window and the DFT; the eigenvalues of
the row covariance), the median/mean of the noise, the CFAR detector of
`cfg.cfar.type` (CA, GO, SO, OS: cells per half window, rank, the factor
from the Pfa by the incomplete beta function or the quadrature of
`cfar_pfa.m`, root finding, the noise scale of the statistic), the tracker
settings (EKF gates by the inverse incomplete gamma function; PMBM with
IMM-UKF: densities, pruning, models, process noise, sigma-point weights,
latency) and the plot extraction settings.

A default input file whose `madeFor` does not match the configuration (for
example after changing `cfg.ant.nElements`, `cfg.fc` or `cfg.adc.fs_Hz`) is
rebuilt by MATLAB; `radar_export_config` computes the same default data in
memory (antenna pattern, DDC filters, DBF matrix; a note says so) and leaves
the file as it is. Only the phase codes and mismatch filters are not
recomputed: it then asks to run `make_input_files(cfg, 'phaseCode')` in
MATLAB once, or to give a file without `madeFor`. Missing input files are
an error as well (the repository has all default files). Checked against
MATLAB's export (Octave 8.4) for the default radar, staggered PRF, phase
codes with mismatch filters, one element, 4 MHz DDC output, the three CFAR
windows with and without a fixed threshold, the CA, GO, SO and OS detectors
(also another OS rank), the FIR clutter filter (fixed and staggered), the EKF
with the CA model in 2-D and the speed gate, the PMBM tracker (3-D, 2-D with
two IMM models, one element), a given CPI length and FFT length, and rebuilt
default files (32 elements, other carrier frequency, other ADC rate): every
record within 3e-14 (relative).

## Options

| Option | Meaning |
|---|---|
| `--config FILE` | configuration from `export_receiver_config.m` (default `config/radar_config.bin`) |
| `--port N` | UDP port of the DDC stream (default 50000) |
| `--record FILE` | also store every received datagram in FILE |
| `--replay FILE` | read the datagrams from a recorded FILE instead of the network (as fast as possible) |
| `--ppi` | show the interactive PPI (Qt build, see below) |
| `--ppi-snapshot FILE` | save the PPI as a PNG at the end and close it |
| `--ppi-show MODE` | what the PPI shows at the start: `plots` (default, after the plot extractor), `detections` (CFAR detections before it), `both`, `video` (the raw video: DDC output averaged over the elements), `dbf` (DBF output of the selected beam) or `mf` (pulse compression output of the selected beam) |
| `--ppi-beam N` | the beam after the DBF of `dbf`, `mf` and the range-Doppler tab at the start (default 1; also in the panel) |
| `--ppi-tab TAB` | tab at the start: `ppi` (default) or `rd` (range-Doppler map of the selected beam next to a small PPI) |
| `--cfar-threshold DB` | CFAR threshold factor 10 log10(alpha) in dB of the configured detector and window (default: the configuration, i.e. the method's Pfa); can also be changed in the PPI |
| `--serial DEVICE` | send the CFAR output to a serial port (`/dev/ttyUSB0`, `COM3`, ...) |
| `--baud N` | serial baud rate (default 921600) |
| `--serial-content X` | `detections` (default: the pruned CFAR detections of every CPI update) or `plots` |
| `--threads N` | processing threads (default: `OMP_NUM_THREADS`, else number of cores - 1) |
| `--mf auto\|direct\|fft` | matched filter: direct valid convolution or FFT fast convolution (same result within float rounding); `auto` (default) takes the faster one for the thread count (with one element: the FFT, see Speed) |
| `--update-every N` | run Doppler, CFAR and extraction every N pulses (default from the configuration: 1) |
| `--plots-csv FILE` | plot list (default `plots.csv`; `''` = none) |
| `--detections-csv FILE` | all pruned CFAR detections (default none) |
| `--tracks-csv FILE` | every track update and deletion: id, confirmed, deleted, time, x, y, z, vx, vy, vz (East, North, up relative to the antenna), hits, misses, ax, ay, az (0 with the constant-velocity model) (default none) |
| `--idle-timeout S` | stop when no PRI arrived for S seconds after the first one (default 0 = never) |
| `--quiet` | do not print every plot and track confirmation (default: `cfg.console.plots` of the configuration). With a low CFAR threshold there can be thousands of plots per second; printing them can take longer than the processing |
| `--print-plots` | print them (overrides `cfg.console.plots = false`) |
| `--queue N` | PRIs buffered between the input thread and the processing (default 64) |
| `--no-skip` | never skip Doppler/CFAR updates when the processing falls behind (see Speed) |
| `--rt-priority N` | Linux: SCHED_FIFO priority N (1 ... 98) of the processing thread and its OpenMP workers, N + 1 for the UDP input thread (the PPI stays at normal priority); see "Linux and real time" |
| `--cpus LIST` | Linux: CPUs of the processing threads, e.g. `2-7` or `2,4,6`; without `--threads` one thread per CPU |
| `--input-cpus LIST` | Linux: CPUs of the UDP input thread (default: those of `--cpus`) |
| `--mlock` | Linux: lock all memory of the process (`mlockall`), so the processing gets no page faults |

Ctrl+C (or closing the PPI window) stops the receiver. The open plots are
closed and the processing times are printed.

## Input: UDP packets

Each datagram is a 64-byte little-endian header followed by int16 (I, Q)
pairs, channel after channel (see "Raw data and UDP format" in
`../README.md`). The receiver reassembles every PRI from its packets in any
order. It ignores duplicates and foreign datagrams, drops PRIs with missing
packets (counted as incomplete), and follows the sender when it restarts.
A jump in the pulse numbers (a PRI lost, or a skipped azimuth sector)
restarts the CPI buffers and closes the open plots, as the MATLAB simulator
does between simulated sectors. A PRI whose layout (channels, DDC samples,
receive interval of each sub-pulse, DDC rate) differs from the configuration
is ignored and counted: export the configuration again. The layout depends
only on the pulse timing, the instrumented range, the DDC filters and the
decimation, so the CFAR method and the pulse compression filters of the
receiver may differ from those of the sender.

The datagrams are read in a thread of their own, which puts the complete
PRIs into a queue (`--queue`). The socket is therefore emptied also while a
Doppler/CFAR update runs. With a live stream a PRI that finds the queue full
is dropped and counted; a replay waits instead.

At the full rate of the default radar (48 channels, about 1400 DDC samples
per PRI and 1000 PRIs/s), the stream is 340 MB/s (2.7 Gbit/s). It needs
10 Gbit Ethernet, and the socket receive buffer should be large. The
receiver asks for 32 MB; on Linux, allow it with
`sysctl -w net.core.rmem_max=33554432`.

## Outputs

**Plots** (console, `--plots-csv`, PPI, serial): time, azimuth, range,
elevation (NaN with one element), radial velocity (ambiguous, + = away),
SNR, sub-pulse, number of CPI detections (hits), beam, range cell, Doppler
bin (0-based) and pulse number.
The SNR is power / noise - 1 with the noise floor of the beam (the median
of its map, scaled to the mean) or, with `cfg.report.snrNoise = 'cfar'`,
the CFAR noise estimate of the cell; the start of the console output says
which (`SNR noise: ...`). **Peak interpolation** (a block of its own after
the CFAR and the pruning, `peak_interpolation()` in `processor.c`, as
`peak_interpolation.m`): with `cfg.report.dopplerInterp` the radial velocity
comes from a parabola through the powers [dB] of the Doppler bins d-1, d,
d+1 of the peak (circular), with `cfg.report.elevationInterp` (default on,
as before) the elevation from a parabola through the beams b-1, b, b+1 when
the peak beam has both neighbours. The detections, the plots, the serial
output and the tracker get the refined values; `dIdx` and `beam` stay the
cells of the peak. Both options are optional in the configuration file
(older files: Doppler off, elevation on); the start of the console output
says which are on, and which clutter filter (`cfg.mti.type`: the canceller
or an FIR filter from `inputs/mti_fir.mat`, given to the receiver as
`mtiCoeffs`).

**Serial frames** (little-endian; CRC-16/CCITT-FALSE, polynomial 0x1021,
init 0xFFFF, over the bytes after the sync word up to the CRC):

| Frame | Layout |
|---|---|
| Detections (type 1), 16 + 24 N bytes, N <= 200 per frame | `A5 5A`, u8 type = 1, u8 pulse (0 short, 1 long), u32 pulse number, f32 azimuth [deg], u16 N, then N records: u16 beam, u16 range cell, u16 Doppler bin, u16 reserved, f32 range [m], f32 radial velocity [m/s], f32 elevation [deg], f32 SNR [dB]; u16 CRC |
| Plot (type 2), 42 bytes | `A5 5A`, u8 type = 2, u8 pulse, u32 pulse number, f64 time [s], f32 azimuth, range, elevation, radial velocity, SNR, u16 beam, u16 hits, u16 CRC |

A CPI update with more than 200 detections is sent as several frames with
the same header. The serial output never blocks the processing. Frames go
into a 256 kB buffer; when the port cannot keep up, whole frames are dropped
and counted. 921600 baud carries about 92 kB/s, which is about 3800
detections per second. `tools/decode_serial.py` decodes a port or a capture
file (`decode_serial.py /dev/ttyUSB0`, `decode_serial.py capture.bin --csv`).

**PPI:** range rings, azimuth spokes (0 = North, clockwise) and the antenna
line. Short-pulse plots are cyan and long-pulse plots yellow. They are drawn
at ground range and, by default, removed when the antenna passes them again
in the next scan.

The window opens maximized and has two tabs. **PPI**: the PPI as above.
**Range-Doppler (selected beam)**: the range-Doppler map of the beam chosen
under *Beam* in the panel (`processor_set_rdmap()`: the power maps of both
sub-pulses after the noise-gain equalisation, the map the CFAR tests,
joined in range like the video in range bins of c / (2 fsOut)), colour = dB
above its noise floor (the median of the cells), from -3 dB (black) to the
strongest cell (at least 20 dB), with a colour scale; the CFAR detections of
the beam (white circles) and the MTI notch (dashed). A pixel column shows
the largest of its range bins, so a peak of one range cell stays visible.
Below the map: the Doppler spectrum (dB above noise) of one range bin, the
one under the cursor or else the bin of the strongest detection (or cell),
with the detections circled, to see a peak and the spectrum around it;
hovering shows range, velocity and dB of the cell. The PPI moves into this
tab, small next to the map, so both are seen live. The worker copies a map
to the display at most 10 times a second, and only while this tab is shown
(nothing is computed otherwise); a map is made when both sub-pulses had a
CPI update. The panel serves both tabs; *Save screenshot* (and
`--ppi-snapshot`) saves the PPI in the first tab and the whole tab (map and
PPI) in the second.

The panel on the right of the window has:
* **Display**: *Plots (after extractor)* as before; *Detections (before
  extractor)*: the pruned CFAR detections of every CPI update, the input of
  the plot extractor, each drawn as an arc over the azimuth step between CPI
  updates (rotation rate x PRI x `--update-every`) at its ground range
  (darker cyan / yellow), so a target appears as an arc about one beamwidth
  long; or both; or *Video (DDC, elements avg.)*: the raw video instead of
  plots and detections. Key `D` switches too (plots, detections, both,
  video). The newest 60 000 detections are kept; they disappear like the
  plots. Video (`video.c`, as `ppi_video_add.m`): while the PPI shows it, the
  worker averages the DDC output of every PRI non-coherently over the
  elements (mean of I^2 + Q^2, before DBF and pulse compression), joins the
  receive intervals of the short pulse (up to the first long-pulse range
  cell) and of the long pulse, puts the result into range bins of
  c / (2 fsOut) (75 m at 2 MHz) and hands it to the PPI (about 25 us per PRI
  with 48 elements, counted as display time; nothing is computed in the other
  modes). The PPI keeps the newest video of every 0.25 deg azimuth bin (the
  mean of the PRIs in it) until the antenna passes again and draws it as an
  image: the colour is the power from 3 dB below the noise floor (black; a
  running 10th percentile of the video, as most of the raw video is clutter)
  over blue, cyan, green, yellow and red to
  50 dB above it (white); the legend shows the scale. A pixel shows the
  largest range bin under it, so an echo stays visible when zoomed out.
  *Clear display* clears the video too.
  *DBF out (selected beam)* and *MF out (selected beam)*: the same display
  for the beam chosen under **Beam after DBF** (default 1, the list shows
  the beam elevations): `processor_set_profile()` makes the processor hand
  over, in every PRI, |x|^2 of the beamformer output of that beam (the pulse
  compression input inside the receive intervals) or |y|^2 of the pulse
  compression output (matched filter, or the mismatch filter of a phase
  code), both sub-pulses joined like the video in the same range bins
  (`ProcessorCallbacks.profile`, a few microseconds per PRI, display time;
  nothing is computed in the other modes). One complex sample per bin has an
  exponential power, so the noise floor is the 10th percentile + 9.77 dB.
  Changing the mode or the beam starts an empty map. Key `D` cycles plots,
  detections, both, video, DBF out, MF out.
* **Save screenshot (PNG) ...**: saves the PPI as it is now (without the
  panel); the default file name has the date and time.
* **Tracks (T)**: shows or hides the tracks of the tracker (EKF or PMBM): a hollow
  triangle at the last update, so the plot (filled circle) lies inside it;
  magenta with its number and a 30 s velocity leader when confirmed, gray
  while tentative. Hover over it: number, state, tracker, range, azimuth, height and
  climb (3-D), speed, heading, acceleration (constant-acceleration model, PMBM),
  for the PMBM the existence probability and the IMM model probabilities,
  updates. A deleted track disappears.
* **CFAR threshold [dB]**: shows the detector and window (`cfg.cfar.type`, `cfg.cfar.method`), the
  threshold the processing uses and its false alarm probability
  (`cfar_pfa()`, needs `nLooks` in the configuration file: export it again
  with the current `export_receiver_config.m`). *Apply* sends the value in
  the field to the processing thread, which uses it from the next CPI update
  on and reports it back (also printed on the console); *Default* returns to
  the value of the configuration.
* **Mean time per PRI [us]** (bottom right, above the help): the mean
  time per PRI over the PRIs since the previous panel update (at least 100;
  the last line says how many), so a change such as a lower CFAR threshold
  shows at once instead of creeping into a mean since the start, of every
  block of the processing thread:
  DBF, matched filter, MTI (clutter filter), Doppler (FFT and non-coherent
  integration), CFAR of the configured detector and window (e.g. `CFAR OS 2-D`,
  in the panel `CFAR OS 2D`) with the pruning, peak interpolation (*Peak
  interp.*), plot extractor, tracking (e.g. *Track 3D PMBM*), the console output (printing the plots, track
  confirmations and messages; 0 with `--quiet`) and the display calls
  (handing the antenna, plots, detections, tracks and times to the PPI),
  and their total: the time the processing thread spends per PRI, display
  calls included. Every label has at most 15 characters (longer names are
  shortened, `1-D` -> `1D`), so the values stay inside the panel. Then
  *PPI draw (GUI)*: the time the GUI thread spends
  taking the data and drawing the PPI, per PRI; it runs on another thread
  (core) at the same time as the processing, so it does not lengthen the
  PRI processing and is not in the total. The worker posts its
  accumulated times 4 times a second; the console and display times are
  also printed at the end (`Output per PRI`).

**CFAR detectors** (`cfar.c`, `cfg.cfar.type` and `cfg.cfar.method` in the
configuration): the same CA (mean of the training cells), GO / SO (larger /
smaller mean of the two half windows: before and after the cell in range
with the 2-D and range windows, below and above it in Doppler with the
Doppler window) and OS (k-th smallest training cell) detectors as
`cfar_detect.m`, each with the 2-D, range or Doppler window. CA keeps its
running box sums; GO and SO use running Doppler sums and range prefix sums
in double; OS counts, for the local peaks only (the other cells cannot
survive the pruning), the training cells below P/alpha, and selects the
k-th one only for a detection (its noise estimate). The false alarm
probability of a threshold (`cfar_pfa()`, the PPI threshold field) is the
closed form for CA and a numerical integral (the same as `cfar_pfa.m`) for
GO, SO and OS.

**Tracking** (`cfg.track.enable` in the configuration, `tracker.c`, which
selects `cfg.track.method`): `tracker_ekf.c`, the same EKF as `track_ekf.m`
(2-D or 3-D, constant velocity or constant acceleration, global nearest
neighbour association on batches of plots), or `tracker_pmbm.c`, the same
multi-target Bayesian tracker as `track_pmbm.m` (PMBM filter: Poisson
undetected targets, multi-Bernoulli mixture with Murty's k best
assignments; IMM-UKF state estimation with the CV, CA and coordinated-turn
models; see `../README.md`). The track updates of both are identical to
MATLAB's (`tests/test_tracker.c` with `../tests/test_tracker.m` and
`../tests/test_tracker_pmbm.m`). The tracks CSV (`--tracks-csv`) has, for the
PMBM, the existence probability and the model probabilities (`pCV, pCA,
pCT`) as its last columns. The radial velocity of the plots is a measurement as
well (`cfg.track.useRadialVelocity`): unambiguous when no target can be
faster than half the span of the Doppler axis (fixed PRF at 170 MHz:
+-441 m/s), otherwise unfolded with the range rate the track predicts. Once
per PRI it gets the plots of that PRI, associates a complete batch and
deletes old tracks (and, with `trackMinSpeed` > 0 from
`cfg.track.minSpeedEnable`, the tracks with 2+ updates slower than
`cfg.track.minSpeed_mps`), as `run_radar_sim.m` does; a gap in the pulse numbers
and the end of the data associate the open batch. Confirmations are printed
on the console (with their radial velocity; `cfg.console.plots`,
`--quiet`), and the number of updates at the end.

The PPI is interactive:

| Action | Effect |
|---|---|
| mouse wheel | zoom in/out around the cursor (up to x500); rings adapt to the visible range |
| left drag | pan |
| double click, `R`, `Home` | reset the view |
| hover over a plot or detection | data tip: range, azimuth, elevation, radial velocity, SNR, time, sub-pulse, CPI hits (plots), beam/cell/bin, pulse number |
| click on a plot or detection | pin its data tip on the display (click again to remove it) |
| right click | options: show plots / detections / both, short/long pulse, plot labels, range rings, antenna trail, plots shown (until the antenna passes again / two scans / all until cleared), minimum SNR, reset view, clear display, remove data tips, pause, save a screenshot |
| `+` `-`, arrows | zoom, pan |
| `C` | remove the pinned data tips |
| `D` | next display mode (plots / detections / both) |
| space | pause/resume the display (the processing continues; the plots are shown on resume) |
| `S` | save a screenshot (PNG) |
| `T` | show/hide the tracks |

The status line shows the zoom and the range, azimuth and east/north
position of the cursor.

## Speed

The processing is written for x86 CPUs with AVX2:
* The data is float and stored as separate real and imaginary arrays, so
  the inner loops of every stage vectorise.
* DBF: three beams are formed per pass over the channel data.
* Matched filter: FFT fast convolution (the valid convolution, within float
  rounding) or direct convolution, whichever is faster for the filter
  length and the number of threads. The FFTs of up to 24 beams run as the
  lanes of one batched FFT. A single element (no DBF) is split into the 8
  interleaved sub-sequences of a four-step FFT (8 FFTs of length N/8 as the
  lanes of one vector, twiddles, FFTs of length 8 over the lanes), so its
  fast convolution is vectorised too: 15 us instead of 55 us for the long
  pulse of the 80 km configuration (300 taps, 771 range cells).
* Only the matched filter outputs that the MTI needs are kept, and the MTI
  output of every pulse is computed once and kept in a ring, stored by
  blocks of 32 cells (`[block][slot][32]`), so that the Doppler processing
  of a block reads contiguous memory.
* Doppler: radix-2^2 FFTs of 32 cells at once (a block of 128 x 32 complex
  values fits the L1 cache).
* CFAR: beam by beam, so the data stays in the cache. The range sums of a
  Doppler row come from partial sums of 1, 2, 4, ... cells, added in one
  pass; only the rows d-oD ... d+1+oD are kept (a ring of 16 rows instead
  of the range sums of the whole map); the Doppler sums are running sums
  (double). The local peak test is part of the CFAR test (see Pruning).
  With fewer beams than threads (one element) the tested range cells are
  split into chunks of at least 64 cells for the threads; each chunk sums
  the cells its windows need, so the sums and decisions are those of the
  whole beam.
* Noise floor (the SNR of the reports): the exact median of the beam's map
  (as MATLAB's `median`), selected among the values of a bracket around the
  median of the previous update (about 1 % of them; one vectorised pass
  counts the values below the bracket and copies those inside it); a sample
  of 2048 values gives the bracket when that one misses. It took more time
  than the CFAR itself with the histogram of the previous version. With
  `cfg.report.snrNoise = 'cfar'` (record `snrNoise` of the configuration;
  older files without it: median) the SNR uses the CFAR noise estimate of
  the cell, which the CFAR keeps for its detections, and no median is
  computed (about 0.03 ms less per PRI with one element; the SNR of a
  target near a strong one or near clutter residue is too low).
* Plot extractor: a detection joins the oldest open group near it, found in
  a hash grid of the cells (beam, range cell, Doppler bin / 8) instead of a
  search through all open groups (with a low threshold there are thousands).
* All stages are split over threads with OpenMP, and all buffers are
  allocated once. A loop never gets more threads than it has work items
  (beams, range chunks, cell blocks). Otherwise every idle thread has to be
  woken and reach the closing barrier; with the PPI and the UDP input on the
  other cores that wait made the single-beam matched filter take more than
  1 ms per PRI instead of about 0.1 ms.
* `radar_bench` measures the time of every block on synthetic PRIs (noise,
  moving targets, optional clutter) for a list of CFAR thresholds, without
  the network, the display and the console: `radar_bench --config
  config/radar_config_1ch.bin --threads 11`. The default thread count is
  all logical CPUs but one: on a CPU with performance and efficiency cores
  (e.g. a Core i5-13420H: 4 + 4 cores, 12 threads, so 11 threads) that
  leaves one logical CPU for the UDP input and the PPI; compare
  `--threads 4`, `6`, `8` and `11` with `radar_bench` and start the
  receiver with the fastest one.

Processing time per PRI with `radar_bench` (synthetic PRIs with moving
targets, Doppler/CFAR on every pulse, mean over 1000 PRIs (300 with 48
elements) and two runs; the 80 km default radar: 262 + 771 range cells;
4-vCPU Intel Xeon VM, built with `-march=native`), before and after the
changes to the CFAR, pruning, noise floor, plot extractor and single-beam
matched filter listed above (the reports, plots and tracks are the same):

| Configuration | Threads | CFAR threshold | before | now |
|---|---|---|---|---|
| 1 element (`config/radar_config_1ch.bin`) | 3 | 12.8 dB (Pfa 1e-8, default) | 1.44 ms | 0.35 ms |
| | 3 | 8 dB (Pfa 2e-3) | 1.53 ms | 0.44 ms |
| | 3 | 6 dB (Pfa 2e-2) | 2.33 ms | 0.87 ms |
| | 1 | 12.8 dB | 1.56 ms | 0.64 ms |
| 48 elements, 21 beams (`config/radar_config.bin`) | 3 | 12.8 dB | 7.9 ms | 6.0 ms |
| | 3 | 8 dB | 24.7 ms | 9.2 ms |

With one element and the default threshold the CFAR block (CFAR, pruning
and reports, the "CFAR" row of the PPI panel) took 1.22 ms of the 1.44 ms,
most of it the median of the noise floor (computed for every sub-pulse with
a detection); now 0.21 ms. The matched filter went from 0.07 to 0.02 ms,
the Doppler processing from 0.14 to 0.11 ms. At 6 dB the plot extractor
took 0.31 ms (now 0.07 ms), with 48 elements at 8 dB 8.2 ms (now 0.65 ms).
With a threshold at which most cells are false alarms (4 dB, Pfa 8e-2:
1100 detections per update, 14 000 plots in a second) the EKF tracker
takes 7 - 8 ms per PRI: it gets thousands of false plots. The replays of
`tests/test_c_receiver.m` (fixed PRF, 50 km: 262 + 371 range cells, 3
threads) now take 0.18 ms per PRI with one element and 2.7 ms with 48
elements.

Measured with the previous version (replay of the six-target streams of
`tests/test_c_receiver.m`: 48 channels, 21 beams, both sub-pulses, 273 +
382 range cells; 4-vCPU Intel Xeon VM at 2.1 GHz, AVX2):

| PRF | Threads | Doppler/CFAR every pulse | every 8 pulses |
|---|---|---|---|
| fixed (CPI 114, FFT 128) | 1 | 8.6 ms | 1.6 ms |
| fixed | 3 | 3.6 ms | 0.80 ms |
| staggered (4 x 28, FFT 32) | 1 | 4.4 ms | 0.91 ms |
| staggered | 3 | 1.7 ms | 0.63 ms |

DBF, matched filter and MTI together take about 0.45 ms per PRI with 3
threads. One Doppler/CFAR update of both sub-pulses takes about 4 ms (fixed
PRF), so with `--update-every 8` the receiver keeps up with 1000 PRIs/s on
this small machine. In a real-time test (`radar_stream_send`, 2.7 Gbit/s
over localhost) it dropped no PRIs and plotted every target, with
`--update-every` 1, 4 and 8. With every-pulse updates, about three quarters
of the updates were skipped automatically. A desktop CPU with 8 or more cores is
correspondingly faster. The 80 km default configuration has about twice as
many long-pulse range cells. Another 4-vCPU VM (Xeon at 2.8 GHz, AVX-512)
was about 1.6 times slower (fixed PRF, 3 threads: 6.0 ms per PRI, 1.35 ms
every 8 pulses; the same with the receiver before the latest fixes). In the
real-time test with `--threads 2` it still dropped no PRIs and plotted every
target: with `--update-every 8` about half of the updates were skipped (fixed
PRF; a third with the staggered PRF), with every-pulse updates nine tenths.

**CFAR detectors: processing time** (`radar_bench`, 4-vCPU VM, 3 threads,
Doppler/CFAR on every pulse, each detector at its threshold for Pfa 1e-8,
2-D window; "CFAR" = CFAR with the pruning and the reports, per PRI):

| Elements | CA | GO | SO | OS (rank 0.9) | total per PRI with CA / OS |
|---|---|---|---|---|---|
| 48 (21 beams) | 3.1 ms | 7.7 ms | 7.9 ms | 5.0 ms | 6.0 / 8.2 ms |
| 1 | 0.22 ms | 0.43 ms | 0.42 ms | 0.45 ms | 0.36 / 0.62 ms |

CA keeps its float box sums; GO and SO need the two half sums of every
cell (double prefix sums, a branch-free loop for the interior cells); OS
counts the training cells of the local peaks only and stops as soon as too
many lie above P/alpha. The PMBM tracker costs about as much as the EKF at
the design threshold (about 1 us per PRI with a few plots per scan).

**Real time.** Doppler processing and CFAR on every new pulse, as in MATLAB,
take longer than a 1 ms PRI on a small PC. The receiver handles this in two
ways:
* `--update-every N` (for example 4 or 8) runs Doppler, CFAR and plot
  extraction every N pulses: 0.14 or 0.29 deg of antenna rotation. The fast
  time processing (DBF, matched filter, MTI) still runs for every pulse.
* Independently of that, when the processing falls behind a live stream
  (more than a quarter of `--queue` waiting), Doppler/CFAR updates are
  skipped until it has caught up. DBF, matched filter and MTI continue, so
  the CPI buffers stay complete. An update after skipped ones counts for all
  the pulses it covers. A plot then still needs detections over the same
  azimuth extent (2 deg), so targets are still plotted as long as some
  updates run while they are in the beam.
  `--no-skip` switches this off.

Without skipped updates the plots are exactly those of MATLAB. The stream of
the default radar is 340 MB/s (2.7 Gbit/s). It needs 10 Gbit Ethernet and a
receive buffer of at least 32 MB (see "Input: UDP packets").

## Linux and real time

The receiver builds and runs on Linux as it is (all tests run on Linux);
nothing has to be changed. For CLion and Qt on Ubuntu 24.04 / Debian 12:

```sh
sudo apt install build-essential cmake qt6-base-dev zlib1g-dev
```

(OpenMP comes with GCC; Qt 5 works too: `qtbase5-dev`.) In CLion open the
`radar_receiver` folder and choose the Release profile. Allow the large UDP
receive buffer (2.7 Gbit/s with 48 elements) and, for `--serial`, the serial
port:

```sh
echo 'net.core.rmem_max=33554432' | sudo tee /etc/sysctl.d/90-radar.conf
sudo sysctl --system
sudo usermod -aG dialout $USER
```

**Ubuntu real-time kernel** (https://ubuntu.com/real-time, PREEMPT_RT).
From Ubuntu 26.04 on it is in the normal archive, free and without Ubuntu
Pro:

```sh
sudo apt install ubuntu-realtime
sudo reboot
uname -v                            # contains PREEMPT_RT
```

Ubuntu 22.04 and 24.04 need Ubuntu Pro (free for personal use on up to 5
machines): `sudo pro attach <token>`, `sudo pro enable realtime-kernel`.
The previous kernel stays installed and can be chosen in the GRUB menu
(Esc or Shift while booting); remove the real-time one with
`sudo apt remove ubuntu-realtime` and its `linux-image-*realtime*`
packages. Before installing it, measure whether it is needed at all: the
average processing time does not get shorter (the kernel bounds the delays,
it does not speed up the processing; throughput can even drop a little),
and with the PRI queue of the input thread delays of a few ms cost no PRI.
`sudo apt install rt-tests; sudo cyclictest -m -p 90 -t -i 200 -D 5m` shows
the largest scheduling delay (`Max`, us) of this computer under its normal
load. The proprietary NVIDIA driver and DKMS modules (VirtualBox, some
Wi-Fi drivers) may not build for the real-time kernel.
(For Intel Core / Atom CPUs Ubuntu also offers `--variant=intel-iotg`, the
same kernel with Intel's time-coordinated computing support; the generic
variant is enough here.) The real-time kernel bounds the time until a
thread of high priority runs, but only threads with a real-time priority
profit from it: the receiver therefore sets it itself (nothing else changes
in the program):

```sh
radar_receiver --config cfg.bin --ppi --rt-priority 40 --cpus 2-7 --input-cpus 1 --mlock
```

* `--rt-priority 40`: the processing thread and its OpenMP workers run with
  SCHED_FIFO 40, the UDP input thread with 41 (it must empty the socket
  even while an update runs). Stay below 50: the interrupt threads of a
  PREEMPT_RT kernel (network card, timer) run at 50 and must not wait for
  the processing. The PPI (Qt GUI thread) keeps the normal priority, so a
  slow redraw never delays the processing.
* `--cpus 2-7`, `--input-cpus 1`: the processing threads only on CPUs 2-7,
  the input thread on CPU 1; the PPI, the network interrupts and the rest of
  the system use the others. The thread count is the number of CPUs given.
  The OpenMP workers inherit priority and CPUs (do not also set
  `OMP_PLACES` / `OMP_PROC_BIND`). On a CPU with performance and efficiency
  cores give the processing the performance cores: `lscpu --extended` shows
  them (the larger `MAXMHZ`; on a Core i5-13420H usually CPUs 0-7 are the
  4 performance cores with their two hardware threads each and CPUs 8-11
  the efficiency cores). Measure the choice with `radar_bench`.
* `--mlock`: all memory is locked (`mlockall`) and freed memory stays in the
  process, so no buffer causes a page fault during the processing.

Without the permission the receiver prints a warning and runs with normal
priority (or unlocked memory). Allow it once for the user (log in again):

```sh
echo "$USER - rtprio 95"          | sudo tee -a /etc/security/limits.d/90-radar.conf
echo "$USER - memlock unlimited"  | sudo tee -a /etc/security/limits.d/90-radar.conf
```

or for the program: `sudo setcap cap_sys_nice,cap_ipc_lock+ep radar_receiver`
(again after every build). Further settings of the system that make the
timing steadier (none is needed to run the receiver):

* keep the processing CPUs free of other work: kernel parameters
  `isolcpus=2-7 nohz_full=2-7 rcu_nocbs=2-7` (in `GRUB_CMDLINE_LINUX`,
  then `sudo update-grub`); the receiver's `--cpus 2-7` still uses them;
* the interrupts of the network card on the CPU of the input thread or next
  to it: `/proc/irq/<n>/smp_affinity_list` (stop `irqbalance`);
* CPU frequency governor `performance` (`cpupower frequency-set -g
  performance`), and in the BIOS/UEFI deep C-states off when the latency
  still varies;
* idle OpenMP workers wait actively for a while (GNU OpenMP: about 0.1 -
  0.3 ms) and with SCHED_FIFO nothing of lower priority runs on their CPU
  meanwhile; on CPUs that are not isolated set `OMP_WAIT_POLICY=passive`
  (they then sleep at once, and waking them costs a few microseconds).

The kernel stops real-time threads for 50 ms in every second when they would
take more than 95 % of a CPU (`kernel.sched_rt_runtime_us`); with the
processing well below the PRI this never happens.

## Sending DDC files or a recording with the PRI timing (`radar_stream_send`)

`radar_stream_send` sends recorded DDC data over UDP at the real PRI rate
(or `--speed X` times faster), so the receiver can be tested in real time
without MATLAB. The source is either

* a folder of DDC files `ddc_<n>.mat`, written by the simulator with
  `cfg.io.saveDdc = true` (one file per PRI: `I`, `Q`, `hdr`), or
* a stream recorded with `radar_receiver --record FILE`.

The data are not processed: every PRI of a DDC file is split into exactly
the datagrams that the simulator sends live (`udp_make_packets.m`), with
`--payload` bytes of samples per datagram; a recording is sent datagram by
datagram as it was received.

```sh
radar_stream_send ../output/ddc_pulses --host 192.168.1.20 --port 50000
radar_stream_send ../output/ddc_pulses --payload 8192 --spread 0.8 --timing-csv timing.csv
radar_receiver --config cfg.bin --port 50000 --record stream.rec   # while MATLAB streams
radar_stream_send stream.rec --port 50001 --loops 10
```

| Option | Default | Meaning |
|---|---|---|
| `--host H`, `--port N` | 127.0.0.1, 50000 | destination (IPv4) |
| `--speed X` | 1 | replay speed (0 = as fast as possible) |
| `--payload BYTES` | 1400 | DDC bytes per datagram (DDC files only; 1400 fits a 1500-byte MTU) |
| `--spread F` | 0 | send the datagrams of a PRI evenly over the fraction F of the PRI instead of as one burst |
| `--no-gaps` | | a skipped sector (gap in the pulse numbers) lasts one PRI instead of its real duration |
| `--catch-up` | | after a late PRI, send the next ones at once until the original schedule is reached (default: keep the PRI spacing, no burst) |
| `--loops N` | 1 | send the data N times; pulse numbers and times continue |
| `--queue N` | 256 | PRIs read ahead (filled before the first PRI is sent) |
| `--timing-csv FILE` | | scheduled and real send time and delay of every PRI |
| `--quiet` | | only the summary |

Timing:
* The send time of every PRI is its PRI start time from the header
  (`hdr.time_s`) divided by `--speed`, on one monotonic clock: errors do not
  accumulate. The sender sleeps until 2 ms before the send time (20 ms on
  Windows, where `Sleep` can wake a timer tick late) and busy-waits for the
  rest. On this 4-vCPU VM, with nothing else running, the PRIs left with a
  median delay of 0.1 us and the PRI intervals were within 13 ... 34 us
  (99 %) of the 1 ms PRI, 1000 PRIs/s at 2.2 Gbit/s with 1400-byte payloads.
  Late PRIs come from the operating system (the VM stalled the process for
  up to a few ms now and then).
* A reader thread reads and splits the files ahead into a queue of
  `--queue` PRIs, so reading never delays a PRI (1 ms PRI: 256 PRIs =
  0.26 s of data, 88 MB for 48 channels). Reading a 48-channel file takes
  less than 0.5 ms.
* At the end a summary gives the number of late PRIs, the largest and the
  median delay and the error of the PRI intervals.
* When the receiver runs on the same PC, leave a core for the sender
  (`--threads`) and use `--update-every` so the receiver keeps up; on this
  4-vCPU VM, with `--threads 2 --update-every 8`, all 397 PRIs of a
  48-channel recording arrived complete and gave the MATLAB plot.

## Tests

| Test | What it checks |
|---|---|
| `tracker` (`tests/test_tracker.c`) | (the tracker of the configuration, `cfg.track.method`) two crossing targets over 6 scans, 3-D and 2-D (plots with radial velocity): two confirmed tracks that stay on their targets; GNN with two targets 200 m apart: each track stays on its own target; speed gate: a 3 m/s track is deleted at its second update with `trackMinSpeed` = 5 m/s, a 100 m/s one kept; with plots, a configuration and the track log from `../tests/test_tracker.m` (CV and CA, 2-D and 3-D, with and without the radial velocity, staggered PRF with folded plot velocities), every C track update equals MATLAB's (same plot, state within 1e-6); with the calls of `../tests/test_tracker_pmbm.m` (16 PMBM cases: maneuvers, constant acceleration, 30 false plots per scan, close and crossing targets, the challenge scenario) every update of the C PMBM tracker equals MATLAB's |
| `core` (`tests/test_core.c`) | batched FFT against a direct DFT; plot extractor, also on random dense detection streams (hundreds of open groups) against the linear search of `plot_extractor.m`: the same plots in the same order; a moving target through the whole chain with the exported configuration (its velocity a quarter of the Doppler axis from zero, for any carrier frequency; one plot at the right range, velocity and beam); direct and fast-convolution matched filters give the same plot; the same target with the CFAR threshold raised to 45 dB at run time gives no detection; `cfar_pfa()` against the closed forms of CA, GO, SO and OS (one look) and the Pfa 1e-8 of the configured threshold; the display taps of the processor (`processor_set_profile`, `processor_set_rdmap`): a DBF and an MF profile every PRI, the MF peak in the target's range bin, the DBF echo above the noise, the range-Doppler map of the target's beam with its peak at the target's range bin and Doppler bin and only that beam's detections; the raw video of the PPI (`video.c`): mean power over the channels, the short-pulse samples up to the handover, then the long-pulse samples, in the right range bins |
| `core_phasecode` | the same core test with `config/radar_config_phasecode.bin` (both sub-pulses P4 phase codes with mismatch filters): the echo is the real code, the plot lies on the target cell |
| `core_1ch` | the same core test with `config/radar_config_1ch.bin` (the default radar with one receive element, no DBF): four-step FFT and direct matched filter give the same plot |
| `cfar` (`tests/test_cfar.c`) | the CFAR detectors of `cfar.c` (CA, GO, SO, OS) against the statistic computed directly from the training cells combined with the direct local peak test: every decision and noise estimate, for every detector with the 2-D, 1-D range and 1-D Doppler windows and tested range intervals in the middle, at both range edges (also when the whole interval lies within the cells of the shifted windows) and over all cells, with the beams of the configuration and with one beam split into range chunks over 4 threads; the noise floor median equals that of a sorted copy (400 maps: odd and even sizes, ties, zeros, a structured half; with and without a bracket from the previous median); the pruning of tens of thousands of CFAR detections (a strong target with sidelobes, 3000 strong cells) gives exactly the cells of the direct pairwise rules; the peak interpolation: a Gaussian peak gives its exact fractional Doppler bin and beam, the Doppler bin wraps at the edge of the axis, the lowest beam keeps its centre, nothing changes with the options off |
| `cfar_1ch` | the same CFAR test with `config/radar_config_1ch.bin` |
| `udp_input` | reassembly with shuffled, duplicate, missing and foreign datagrams; timeouts; record and replay; fuzzing |
| `ddc_files` (`tests/test_ddc_files.c`) | MAT-file reader with files written as MATLAB writes them (data in smaller storage types, small elements, skipped char field, compressed `-v7`: with zlib, or stored blocks for the built-in inflate without zlib; with zlib the built-in inflate is checked against zlib's dynamic and fixed Huffman blocks); folder listing; the datagrams of a DDC file are accepted by `udp_input.c` and give back the same PRI; corrupted, truncated, `-v7.3` and big-endian files are refused |
| `mfile` (`tests/test_mfile.c`) | the MATLAB subset of `radar_export_config`: function line, comments, block comments, continuation lines, matrices with spaces and unary minus, multi-line cells, `struct`, `fullfile`/`fileparts`/`mfilename`, precedence (`-2^2`), indexing, functions, `--set` statements, unevaluated MATLAB functions, subfunctions; errors with the line |
| `export_config`, `core_exported` | `radar_export_config` writes the configuration of `../radar_config.m`, and the core test runs the whole chain with it |
| `serial` | frame encoding and CRC; real output to a pseudo-terminal (Linux): exact bytes, 200-detection split, non-blocking drops, close |
| `ppi` | the Qt PPI runs off-screen with antenna updates, plots, detections and a CFAR threshold report and saves a snapshot |
| `ppi_interactive` | mouse, wheel and key events: data tips, cursor readout, zoom around the cursor, pan, pin/unpin, reset, filters, persistence options, pause, clear; detections (display modes, arcs drawn, data tips, key D); window panel: display mode box, CFAR threshold Apply/Default and its report, screenshot, mean block times (labels of at most 15 characters, the peak interpolation row); raw video (box and key D, rows of 160 azimuth bins, noise floor, echo drawn, nothing where the antenna did not pass, plots hidden, Clear); DBF out / MF out (beam selector, profile wanted only in these modes and for the selected beam, rows of another kind or beam not shown, noise floor of one-sample rows, echo bright); range-Doppler tab (map wanted only in the tab, the PPI moved next to it and back, map drawn with the peak bright and the noise dark, spectrum at the strongest detection, hover readout, a map of another beam ignored, screenshot of the tab); tracks (triangle drawn around the plot, data tip, PMBM data tip with existence and model probabilities, checkbox and key T, deletion) |
| `../tests/test_c_receiver.m` (Octave/MATLAB) | the six-target scenario streamed over UDP to this receiver: the C plots must match the MATLAB plots; the recording is replayed to check that the result is repeatable and to measure the time |
