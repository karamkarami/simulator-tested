# Radar receiver in C (CLion / CMake project)

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
 CA-CFAR (2-D, 1-D range or 1-D Doppler window, tested range interval, shifted range window at the edges)
   -> pruning (local peaks, sidelobe blanking) -> reports -> plot extraction -> EKF tracking (optional)
        |                                   |                         |                |
        |                                   +-> serial (detections)   +-> CSV, console, Qt PPI,
        |                                                                 serial (plots) +-> CSV, Qt PPI
```

The algorithms are the same as in the MATLAB simulator (`../src/*.m`, see
`../README.md`). The C receiver reads every number that MATLAB derives
(DBF matrix, pulse compression filters, range axes, CPI length, MTI
weights, Doppler window and noise gains, CFAR method, window, threshold
factor and tested range cells, extraction and tracker settings) from a
configuration file that `export_receiver_config.m` writes. Both programs
therefore always use the same values.

**Phase codes and mismatch filters** (`cfg.wave.<pulse>.type =
'phasecode'`, `cfg.wave.<pulse>.compression = 'matched' | 'mismatch'`,
codes and filters in `../inputs/phase_code.mat`): the exported taps are the
matched or the mismatch filter at the DDC rate, and the exported range axis
already holds the earlier peak of the longer mismatch filter. The receiver
filters with whatever taps the file has (direct or FFT, as for LFM), so it
needs no setting of its own; it prints the waveform and the filter of each
sub-pulse at the start, e.g. `long pulse: phase code of 150 chips, mismatch
filter of 320 taps (23.65 ... 80.61 km)`. The file also holds the pulse
itself at the DDC rate (`<pulse>_pulse`), which `test_core` uses as the
echo. `config/radar_config_phasecode.bin` is such a configuration (both
sub-pulses with mismatch filters) for the `core_phasecode` test.

## Build

Requirements: a C11 compiler (GCC, Clang, or MSVC 2019 16.8 or newer), CMake 3.21 or newer,
OpenMP (optional, for multi-threading), Qt 5.12+ or Qt 6 Widgets
(optional, for the PPI) and zlib (optional: `radar_stream_send` then also
reads compressed `-v7` DDC files; `-v6` files need nothing).

**CLion:** *File > Open* the `radar_receiver` folder; CLion loads
`CMakeLists.txt`. Choose the **Release** profile (*Settings > Build,
Execution, Deployment > CMake*, add a Release profile). Without it the code
is not optimised and runs about 10 times slower. Run configurations are
created for `radar_receiver` and the tests. Set the program arguments in the
run configuration, for example `--config ../config/radar_config.bin --ppi`.
If CMake does not find Qt, add `-DCMAKE_PREFIX_PATH=<Qt>/lib/cmake` (for
example `C:/Qt/6.7.2/mingw_64/lib/cmake`) to the CMake options.

**Command line:**

```sh
cmake -S radar_receiver -B radar_receiver/build -DCMAKE_BUILD_TYPE=Release
cmake --build radar_receiver/build -j
cd radar_receiver/build && ctest --output-on-failure
```

CMake options:

| Option | Default | Meaning |
|---|---|---|
| `RADAR_WITH_QT` | ON | build the PPI (the build continues without Qt when it is not found) |
| `RADAR_OPENMP` | ON | OpenMP threads |
| `RADAR_NATIVE` | ON | optimise for the CPU of the build machine (`-march=native`, MSVC `/arch:AVX2`); switch it off when the program must run on other PCs |
| `RADAR_BUILD_TESTS` | ON | unit tests (`ctest`) |

On Windows, the MinGW toolchain bundled with CLion and MSVC both work.
OpenMP and `-march=native` need no extra setup with MinGW. With MSVC,
`/openmp` is used, and `RADAR_NATIVE` means `/arch:AVX2` (every Intel Core
CPU since 2013). When Qt is found, `windeployqt` copies the Qt DLLs next to
the programs after the build, so they also start outside Qt's environment.
The serial port is opened as `COM3` and so on.

## Use with the MATLAB simulator

```matlab
cfg = radar_config;                       % the same cfg for both programs
cfg.udp.enable = true;                    % stream the DDC output
export_receiver_config(cfg);              % -> radar_receiver/config/radar_config.bin
```

```sh
radar_receiver/build/radar_receiver --config radar_receiver/config/radar_config.bin --ppi
```

Then run `run_radar_sim(cfg)` in MATLAB. Export the configuration again
whenever `cfg` changes: the receiver checks that the PRI layout (number of
channels and DDC samples) matches and ignores PRIs that do not.
`tests/test_c_receiver.m` does all of this automatically and compares the
plots of both programs.

## Options

| Option | Meaning |
|---|---|
| `--config FILE` | configuration from `export_receiver_config.m` (default `config/radar_config.bin`) |
| `--port N` | UDP port of the DDC stream (default 50000) |
| `--record FILE` | also store every received datagram in FILE |
| `--replay FILE` | read the datagrams from a recorded FILE instead of the network (as fast as possible) |
| `--ppi` | show the interactive PPI (Qt build, see below) |
| `--ppi-snapshot FILE` | save the PPI as a PNG at the end and close it |
| `--ppi-show MODE` | what the PPI shows at the start: `plots` (default, after the plot extractor), `detections` (CFAR detections before it) or `both` |
| `--cfar-threshold DB` | CA-CFAR threshold factor 10 log10(alpha) in dB of the configured method (default: the configuration, i.e. the method's Pfa); can also be changed in the PPI |
| `--serial DEVICE` | send the CFAR output to a serial port (`/dev/ttyUSB0`, `COM3`, ...) |
| `--baud N` | serial baud rate (default 921600) |
| `--serial-content X` | `detections` (default: the pruned CFAR detections of every CPI update) or `plots` |
| `--threads N` | processing threads (default: `OMP_NUM_THREADS`, else number of cores - 1) |
| `--mf auto\|direct\|fft` | matched filter: direct valid convolution or FFT fast convolution (same result); `auto` (default) takes the faster one for the thread count |
| `--update-every N` | run Doppler, CFAR and extraction every N pulses (default from the configuration: 1) |
| `--plots-csv FILE` | plot list (default `plots.csv`; `''` = none) |
| `--detections-csv FILE` | all pruned CFAR detections (default none) |
| `--tracks-csv FILE` | every track update and deletion: id, confirmed, deleted, time, x, y, z, vx, vy, vz (East, North, up relative to the antenna), hits, misses, ax, ay, az (0 with the constant-velocity model) (default none) |
| `--idle-timeout S` | stop when no PRI arrived for S seconds after the first one (default 0 = never) |
| `--quiet` | do not print every plot |
| `--queue N` | PRIs buffered between the input thread and the processing (default 64) |
| `--no-skip` | never skip Doppler/CFAR updates when the processing falls behind (see Speed) |

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
samples of each sub-pulse, DDC rate) differs from the configuration is
ignored and counted: export the configuration again.

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

The panel on the right of the window has:
* **Display**: *Plots (after extractor)* as before; *Detections (before
  extractor)*: the pruned CFAR detections of every CPI update, the input of
  the plot extractor, each drawn as an arc over the azimuth step between CPI
  updates (rotation rate x PRI x `--update-every`) at its ground range
  (darker cyan / yellow), so a target appears as an arc about one beamwidth
  long; or both. Key `D` switches too. The newest 60 000 detections are kept;
  they disappear like the plots.
* **Save screenshot (PNG) ...**: saves the PPI as it is now (without the
  panel); the default file name has the date and time.
* **Tracks (T)**: shows or hides the tracks of the EKF tracker: a hollow
  triangle at the last update, so the plot (filled circle) lies inside it;
  magenta with its number and a 30 s velocity leader when confirmed, gray
  while tentative. Hover over it: number, state, range, azimuth, height and
  climb (3-D), speed, heading, acceleration (constant-acceleration model),
  updates. A deleted track disappears.
* **CFAR threshold [dB]**: shows the method (`cfg.cfar.method`), the
  threshold the processing uses and its false alarm probability
  (`cfar_pfa()`, needs `nLooks` in the configuration file: export it again
  with the current `export_receiver_config.m`). *Apply* sends the value in
  the field to the processing thread, which uses it from the next CPI update
  on and reports it back (also printed on the console); *Default* returns to
  the value of the configuration.
* **Mean time per PRI [us]** (bottom right, above the help): the mean
  time per PRI since the start of every block of the processing thread:
  DBF, matched filter, MTI (clutter filter), Doppler (FFT and non-coherent
  integration), CFAR of the configured method with the pruning, plot
  extractor, tracking, the console output (printing the plots, track
  confirmations and messages; 0 with `--quiet`) and the display calls
  (handing the antenna, plots, detections, tracks and times to the PPI),
  and their total. Then *PPI drawing (GUI)*: the time the GUI thread spends
  taking the data and drawing the PPI, per PRI; it runs on another core at
  the same time, so it is not in the total. The worker posts its
  accumulated times 4 times a second; the console and display times are
  also printed at the end (`Output per PRI`).

**Tracking** (`cfg.track.enable` in the configuration, `tracker.c`): the
same EKF as `track_update.m` (2-D or 3-D, constant velocity or constant
acceleration, global nearest neighbour association on batches of plots,
see `../README.md`). Once per PRI it gets the plots of that PRI, associates
a complete batch and deletes old tracks, as `run_radar_sim.m` does; a gap
in the pulse numbers and the end of the data associate the open batch.
Confirmations are printed on the console (unless `--quiet`), and the number
of updates at the end.

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
* Matched filter: FFT fast convolution (same result as the valid
  convolution) or direct convolution, whichever is faster for the filter
  length and the number of threads.
* Only the matched filter outputs that the MTI needs are kept, and the MTI
  output of every pulse is computed once and kept in a ring.
* Doppler: radix-2^2 FFTs of many cells at once, in cache-sized blocks.
* CFAR: beam by beam, so the data stays in the cache. Range sums are built
  from partial sums, and Doppler sums are running sums. The noise floor
  median uses a histogram.
* All stages are split over threads with OpenMP, and all buffers are
  allocated once. A loop never gets more threads than it has work items
  (beams, cell blocks): with a single element (no DBF) the loops over beams
  (matched filter, CFAR) run in the processing thread. Otherwise every idle
  thread has to be woken and reach the closing barrier; with the PPI and the
  UDP input on the other cores that wait made the single-beam matched filter
  take more than 1 ms per PRI instead of about 0.1 ms.

Measured processing time per PRI (replay of the six-target streams of
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
| `tracker` (`tests/test_tracker.c`) | two crossing targets over 6 scans, 3-D and 2-D: two confirmed tracks that stay on their targets; GNN with two targets 200 m apart: each track stays on its own target; with plots, a configuration and the track log from `../tests/test_tracker.m` (CV and CA, 2-D and 3-D), every C track update equals MATLAB's (same plot, state within 1e-6) |
| `core` (`tests/test_core.c`) | batched FFT against a direct DFT; plot extractor; a moving target through the whole chain with the exported configuration (one plot at the right range, velocity and beam); direct and fast-convolution matched filters give the same plot; the same target with the CFAR threshold raised to 45 dB at run time gives no detection; `cfar_pfa()` against the closed form (one look) and the Pfa 1e-8 of the configured threshold |
| `core_phasecode` | the same core test with `config/radar_config_phasecode.bin` (both sub-pulses P4 phase codes with mismatch filters): the echo is the real code, the plot lies on the target cell |
| `cfar` (`tests/test_cfar.c`) | the CA-CFAR of `processor.c` against direct sums over the training cells: every decision, for the 2-D, 1-D range and 1-D Doppler windows and tested range intervals in the middle, at both range edges (also when the whole interval lies within the cells of the shifted windows) and over all cells |
| `udp_input` | reassembly with shuffled, duplicate, missing and foreign datagrams; timeouts; record and replay; fuzzing |
| `ddc_files` (`tests/test_ddc_files.c`) | MAT-file reader with files written as MATLAB writes them (data in smaller storage types, small elements, skipped char field, compressed `-v7` with zlib); folder listing; the datagrams of a DDC file are accepted by `udp_input.c` and give back the same PRI; corrupted, truncated, `-v7.3` and big-endian files are refused |
| `serial` | frame encoding and CRC; real output to a pseudo-terminal (Linux): exact bytes, 200-detection split, non-blocking drops, close |
| `ppi` | the Qt PPI runs off-screen with antenna updates, plots, detections and a CFAR threshold report and saves a snapshot |
| `ppi_interactive` | mouse, wheel and key events: data tips, cursor readout, zoom around the cursor, pan, pin/unpin, reset, filters, persistence options, pause, clear; detections (display modes, arcs drawn, data tips, key D); window panel: display mode box, CFAR threshold Apply/Default and its report, screenshot, mean block times; tracks (triangle drawn around the plot, data tip, checkbox and key T, deletion) |
| `../tests/test_c_receiver.m` (Octave/MATLAB) | the six-target scenario streamed over UDP to this receiver: the C plots must match the MATLAB plots; the recording is replayed to check that the result is repeatable and to measure the time |
