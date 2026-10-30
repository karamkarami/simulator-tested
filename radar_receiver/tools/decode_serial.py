#!/usr/bin/env python3
"""decode_serial.py - decode the serial output of the radar receiver.

Reads the frames written by serial_out.c from a serial port or from a binary
file (a capture of the port), finds them by the sync bytes A5 5A and the
CRC-16/CCITT-FALSE, and prints one line per detection (frame type 1) or plot
(frame type 2). Frames with a wrong CRC are skipped (the decoder resyncs on
the next sync bytes). A summary goes to stderr at the end (or on Ctrl-C).

Usage:
  decode_serial.py /dev/ttyUSB0             serial port, 921600 baud
  decode_serial.py COM3 --baud 921600       Windows (needs pyserial)
  decode_serial.py capture.bin              binary file
  decode_serial.py - < capture.bin          standard input
  decode_serial.py capture.bin --csv        CSV output

Serial ports are read with pyserial when it is installed; on Linux/macOS
the port is also read without it (termios, raw 8N1).

Frame format (little-endian):
  Detections (type 1), 16 + 24*N bytes:
     0: A5 5A   2: u8 type = 1   3: u8 pulse (0 short, 1 long)   4: u32 pulse number
     8: f32 azimuth [deg]   12: u16 N
    14: N x (u16 beam, u16 range cell, u16 Doppler bin, u16 reserved,
             f32 range [m], f32 radial velocity [m/s], f32 elevation [deg], f32 SNR [dB])
    14+24N: u16 CRC
  Plot (type 2), 42 bytes:
     0: A5 5A   2: u8 type = 2   3: u8 pulse   4: u32 pulse number   8: f64 time [s]
    16: f32 azimuth, range, elevation, radial velocity, SNR   36: u16 beam   38: u16 hits
    40: u16 CRC
  CRC over the bytes after the sync up to the byte before the CRC.
"""

import argparse
import binascii
import os
import re
import signal
import stat
import struct
import sys

SYNC = b"\xa5\x5a"
TYPE_DETECTIONS = 1
TYPE_PLOT = 2
DET_HEADER = struct.Struct("<2sBBIfH")         # 14 bytes: sync, type, pulse, pulse number, az, N
DET_RECORD = struct.Struct("<HHHHffff")        # 24 bytes
PLOT_FRAME = struct.Struct("<2sBBIdfffffHHH")  # 42 bytes, CRC included
CRC_LEN = 2
PULSE_NAMES = {0: "short", 1: "long"}

CSV_HEADER = ("kind,pulse_number,pulse,time_s,az_deg,beam,range_cell,doppler_bin,"
              "range_m,vr_mps,el_deg,snr_db,hits")


def crc16_ccitt(data):
    """CRC-16/CCITT-FALSE (poly 0x1021, init 0xFFFF)."""
    return binascii.crc_hqx(data, 0xFFFF)


def decode_frame(frame):
    """Decodes a frame with a valid CRC into a dict."""
    if frame[2] == TYPE_PLOT:
        (_, _, pulse, pn, time_s, az, rng, el, vr, snr, beam, hits, _) = PLOT_FRAME.unpack(frame)
        return {"type": TYPE_PLOT, "pulse": pulse, "pulse_number": pn, "time_s": time_s,
                "az_deg": az, "range_m": rng, "el_deg": el, "vr_mps": vr, "snr_db": snr,
                "beam": beam, "hits": hits}
    _, _, pulse, pn, az, n = DET_HEADER.unpack_from(frame, 0)
    dets = []
    for k in range(n):
        beam, rcell, dbin, _, rng, vr, el, snr = DET_RECORD.unpack_from(
            frame, DET_HEADER.size + k * DET_RECORD.size)
        dets.append({"beam": beam, "range_cell": rcell, "doppler_bin": dbin, "range_m": rng,
                     "vr_mps": vr, "el_deg": el, "snr_db": snr})
    return {"type": TYPE_DETECTIONS, "pulse": pulse, "pulse_number": pn, "az_deg": az,
            "detections": dets}


class FrameDecoder:
    """Finds frames in a byte stream (sync + length + CRC)."""

    def __init__(self, max_dets=200):
        self.buf = bytearray()
        self.max_dets = max_dets
        self.det_frames = 0
        self.plot_frames = 0
        self.detections = 0
        self.rejected = 0      # sync bytes that did not start a valid frame
        self.skipped = 0       # bytes outside valid frames
        self.incomplete = 0    # bytes still waiting for the end of a frame when the input ended

    def feed(self, data):
        """Adds received bytes and returns the list of complete decoded frames."""
        self.buf += data
        frames = []
        while True:
            frame = self._next_frame()
            if frame is None:
                return frames
            frames.append(frame)

    def finish(self):
        """Call at the end of the input. The bytes still waiting for the end of a
        frame may start with a false sync that hides complete frames behind it:
        they are scanned again byte by byte. Returns the frames found."""
        self.incomplete += len(self.buf)
        frames = []
        while self.buf:
            frames.extend(self.feed(b""))
            if self.buf:                  # starts with an incomplete frame: skip its first byte
                self.skipped += 1
                del self.buf[:1]
        return frames

    def _frame_length(self):
        """Length of the frame at the start of buf: 0 = need more bytes, -1 = not a frame."""
        if len(self.buf) < 4:
            return 0
        if self.buf[2] == TYPE_PLOT:
            return PLOT_FRAME.size
        if self.buf[2] == TYPE_DETECTIONS:
            if len(self.buf) < DET_HEADER.size:
                return 0
            n = self.buf[12] | (self.buf[13] << 8)
            if n > self.max_dets:
                return -1
            return DET_HEADER.size + n * DET_RECORD.size + CRC_LEN
        return -1

    def _next_frame(self):
        while True:
            start = self.buf.find(SYNC)
            if start < 0:
                keep = 1 if self.buf.endswith(SYNC[:1]) else 0   # maybe the first sync byte
                self.skipped += len(self.buf) - keep
                del self.buf[:len(self.buf) - keep]
                return None
            if start > 0:
                self.skipped += start
                del self.buf[:start]

            length = self._frame_length()
            if length == 0 or len(self.buf) < length:
                return None                                       # wait for more bytes
            if length > 0:
                frame = bytes(self.buf[:length])
                crc = frame[-2] | (frame[-1] << 8)
                if crc == crc16_ccitt(frame[2:-CRC_LEN]):
                    del self.buf[:length]
                    return self._count(decode_frame(frame))
            self.rejected += 1                                    # resync after this byte
            self.skipped += 1
            del self.buf[:1]

    def _count(self, frame):
        if frame["type"] == TYPE_PLOT:
            self.plot_frames += 1
        else:
            self.det_frames += 1
            self.detections += len(frame["detections"])
        return frame

    def summary(self):
        return ("frames: %d detections frames (%d detections), %d plots; "
                "rejected syncs: %d, skipped bytes: %d, incomplete bytes at the end: %d"
                % (self.det_frames, self.detections, self.plot_frames, self.rejected,
                   self.skipped, self.incomplete + len(self.buf)))


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def pulse_name(pulse):
    return PULSE_NAMES.get(pulse, str(pulse))


def format_frame(frame, csv):
    """Returns the output lines of one frame."""
    if frame["type"] == TYPE_PLOT:
        f = frame
        if csv:
            return ["plot,%d,%s,%.6f,%.3f,%d,,,%.2f,%.3f,%.3f,%.2f,%d" % (
                f["pulse_number"], pulse_name(f["pulse"]), f["time_s"], f["az_deg"], f["beam"],
                f["range_m"], f["vr_mps"], f["el_deg"], f["snr_db"], f["hits"])]
        return ["PLOT pn=%-10d t=%12.6f s %-5s az=%7.2f range=%9.1f m el=%6.2f vr=%7.2f m/s "
                "snr=%5.1f dB beam=%2d hits=%d" % (
                    f["pulse_number"], f["time_s"], pulse_name(f["pulse"]), f["az_deg"],
                    f["range_m"], f["el_deg"], f["vr_mps"], f["snr_db"], f["beam"], f["hits"])]
    lines = []
    for d in frame["detections"]:
        if csv:
            lines.append("det,%d,%s,,%.3f,%d,%d,%d,%.2f,%.3f,%.3f,%.2f," % (
                frame["pulse_number"], pulse_name(frame["pulse"]), frame["az_deg"], d["beam"],
                d["range_cell"], d["doppler_bin"], d["range_m"], d["vr_mps"], d["el_deg"],
                d["snr_db"]))
        else:
            lines.append("DET  pn=%-10d %-5s az=%7.2f beam=%2d rcell=%5d dbin=%4d range=%9.1f m "
                         "el=%6.2f vr=%7.2f m/s snr=%5.1f dB" % (
                             frame["pulse_number"], pulse_name(frame["pulse"]), frame["az_deg"],
                             d["beam"], d["range_cell"], d["doppler_bin"], d["range_m"],
                             d["el_deg"], d["vr_mps"], d["snr_db"]))
    return lines


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def is_serial_port(path):
    if re.match(r"^(\\\\\.\\)?COM\d+$", path, re.IGNORECASE):
        return True
    try:
        return stat.S_ISCHR(os.stat(path).st_mode)
    except OSError:
        return False


def read_pyserial(path, baud):
    import serial  # optional dependency
    port = serial.Serial(path, baudrate=baud, bytesize=8, parity="N", stopbits=1, timeout=0.2)
    try:
        while True:
            data = port.read(max(1, port.in_waiting))
            if data:
                yield data
    finally:
        port.close()


def read_termios(path, baud):
    import termios
    speed = getattr(termios, "B%d" % baud, None)
    if speed is None:
        raise SystemExit("baud rate %d is not supported without pyserial" % baud)
    fd = os.open(path, os.O_RDONLY | os.O_NOCTTY)
    try:
        iflag, oflag, cflag, lflag, _, _, cc = termios.tcgetattr(fd)
        iflag &= ~(termios.IGNBRK | termios.BRKINT | termios.PARMRK | termios.ISTRIP |
                   termios.INLCR | termios.IGNCR | termios.ICRNL | termios.IXON | termios.IXOFF)
        oflag &= ~termios.OPOST
        lflag &= ~(termios.ECHO | termios.ECHONL | termios.ICANON | termios.ISIG | termios.IEXTEN)
        cflag &= ~(termios.CSIZE | termios.PARENB | termios.CSTOPB | getattr(termios, "CRTSCTS", 0))
        cflag |= termios.CS8 | termios.CREAD | termios.CLOCAL
        cc[termios.VMIN] = 1
        cc[termios.VTIME] = 0
        termios.tcsetattr(fd, termios.TCSANOW, [iflag, oflag, cflag, lflag, speed, speed, cc])
        while True:
            try:
                data = os.read(fd, 65536)
            except OSError:      # EIO: the other side of a pseudo-terminal was closed
                return
            if not data:
                return
            yield data
    finally:
        os.close(fd)


def read_file(stream):
    while True:
        data = stream.read(65536)
        if not data:
            return
        yield data


def open_source(path, baud):
    """Returns (generator of byte chunks, live) for a port, a file or '-'."""
    if path == "-":
        return read_file(sys.stdin.buffer), False
    if not is_serial_port(path):
        return read_file(open(path, "rb")), False
    try:
        import serial  # noqa: F401
        return read_pyserial(path, baud), True
    except ImportError:
        if os.name == "nt":
            raise SystemExit("reading %s needs pyserial: pip install pyserial" % path)
        return read_termios(path, baud), True


def stop(signum, frame):
    raise KeyboardInterrupt   # SIGTERM: print the summary like Ctrl-C


def main():
    parser = argparse.ArgumentParser(description="Decode the serial output of the radar receiver.")
    parser.add_argument("source", help="serial port (/dev/ttyUSB0, COM3), binary file or - (stdin)")
    parser.add_argument("--baud", type=int, default=921600, help="baud rate (default 921600)")
    parser.add_argument("--max-dets", type=int, default=200,
                        help="largest number of detections accepted in one frame (default 200)")
    parser.add_argument("--csv", action="store_true", help="print CSV lines")
    args = parser.parse_args()

    signal.signal(signal.SIGTERM, stop)
    chunks, live = open_source(args.source, args.baud)
    decoder = FrameDecoder(args.max_dets)
    if args.csv:
        print(CSV_HEADER)
    def emit(frames):
        lines = []
        for frame in frames:
            lines.extend(format_frame(frame, args.csv))
        if lines:
            print("\n".join(lines), flush=live)

    try:
        try:
            for data in chunks:
                emit(decoder.feed(data))
        except KeyboardInterrupt:
            pass
        emit(decoder.finish())
    except BrokenPipeError:
        sys.stdout = open(os.devnull, "w")   # output closed (e.g. piped into head)
        return 0
    print(decoder.summary(), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
