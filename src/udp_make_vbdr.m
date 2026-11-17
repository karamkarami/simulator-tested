function packets = udp_make_vbdr(hdr, I, Q, vbdr, loss)
%UDP_MAKE_VBDR  The DDC output of one PRI as the stream of the VBDR board.
%   packets = UDP_MAKE_VBDR(hdr, I, Q, vbdr, loss) returns the nPackets UDP packets (cell
%   array of uint8 row vectors) of the PRI described by hdr (pulse_header.m), I, Q int16
%   [nSamples x nChannels] (the first 4 channels are sent, missing ones are zero).
%   vbdr: payloadBytes, nPackets (setup_radar.m: sys.vbdr; nPackets [] = computed);
%   loss: packets the sender lost so far (label field "Loss=").
%
%   Stream (little-endian): two halves, half 1 = CH1 & CH2, half 2 = CH3 & CH4. A half is
%     32 uint16 system parameters:
%        1 PRF number (pulse number mod 65536)   2 azimuth (0 ... 65535 = 0 ... 360 deg)
%        3 N/A (test, 0)                         4 0x5555
%        5 previous PRT [us]                     6 current PRT [us]
%        7, 8 N/A (test, 0)                      9 ... 12 Code1[0:63] (0)   13 ... 32 Code2[0:319] (0)
%     followed by int16 CHa_I, CHa_Q, CHb_I, CHb_Q of every DDC sample (from the start of the
%     DDC record). Each half is cut into nPackets/2 packets of payloadBytes bytes (the last one
%     shorter); every packet starts with the 28-byte label:
%        uint32 stream number (the pulse number), uint32 packet number (1 ... nPackets),
%        uint32 packets per stream, 16 bytes ASCII "Loss=<n>" (zero padded).
nS = size(I, 1);
nCh = size(I, 2);
I4 = zeros(nS, 4, 'int16');  Q4 = I4;
m = min(nCh, 4);
I4(:, 1:m) = int16(I(:, 1:m));
Q4(:, 1:m) = int16(Q(:, 1:m));

prevPri = hdr.pri_s;
if isfield(hdr, 'prevPri_s') && ~isempty(hdr.prevPri_s), prevPri = hdr.prevPri_s; end
par = zeros(1, 32, 'uint16');
par(1) = mod(hdr.pulseNumber, 65536);
par(2) = mod(round(mod(hdr.antennaAz_deg, 360) / 360 * 65536), 65536);
par(4) = hex2dec('5555');
par(5) = min(65535, round(prevPri * 1e6));
par(6) = min(65535, round(hdr.pri_s * 1e6));
parBytes = typecast(par, 'uint8');

payload = vbdr.payloadBytes;
halfBytes = 64 + 8*nS;
nHalf = ceil(halfBytes / payload);
if isfield(vbdr, 'nPackets') && ~isempty(vbdr.nPackets) && vbdr.nPackets ~= 2*nHalf
    error('udp_make_vbdr: %d DDC samples need %d packets of %d bytes, not nPackets = %d.', nS, 2*nHalf, ...
          payload, vbdr.nPackets);
end
lossTxt = zeros(1, 16, 'uint8');
t = uint8(sprintf('Loss=%d', loss));
lossTxt(1:min(16, numel(t))) = t(1:min(16, numel(t)));

packets = cell(1, 2*nHalf);
for h = 1:2
    ch = 2*h - 1;
    d = [I4(:, ch).'; Q4(:, ch).'; I4(:, ch+1).'; Q4(:, ch+1).'];    % CHa_I, CHa_Q, CHb_I, CHb_Q per sample
    bytes = [parBytes, typecast(d(:).', 'uint8')];
    for j = 1:nHalf
        k = (h-1)*nHalf + j;
        part = bytes((j-1)*payload + 1 : min(j*payload, halfBytes));
        packets{k} = [typecast(uint32([hdr.pulseNumber, k, 2*nHalf]), 'uint8'), lossTxt, part];
    end
end
end
