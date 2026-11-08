function rd = ppi_rd_map(V, sys, slices, dets, beam, t_s, az_deg)
%PPI_RD_MAP  Range-Doppler map of one beam for the second tab of the PPI.
%   rd = PPI_RD_MAP(V, sys, slices, dets, beam, t_s, az_deg) joins the power
%   maps of the beam of both sub-pulses (slices{p} = [nRange x nfft] of
%   doppler_process.m, after the noise-gain equalisation, p = 1 short, 2 long;
%   [] if none yet) in the range bins of the PPI map (ppi_video_init.m): the
%   short-pulse cells nearer than the first long-pulse cell, then the
%   long-pulse cells. rd.P [nBins x nfft] (NaN = no cell), rd.vel_mps (the
%   Doppler axis), rd.noise (the median of the cells: the noise floor),
%   rd.dets (the reports of the beam: [range_m dopplerBin power]),
%   rd.notchBins (MTI notch). The C receiver does the same (processor.c,
%   rdmap_add).
nfft = sys.dop.nfft;
rd.P = NaN(V.nR, nfft);
for p = 1:2
    if isempty(slices{p}), continue; end
    L = V.mf{p};
    rd.P(L.bin, :) = double(slices{p}(L.idx, :));
end
rd.vel_mps = -sys.dop.f_Hz * sys.lambda / 2;
c = rd.P(~isnan(rd.P));
if isempty(c), rd.noise = NaN; else, rd.noise = median(c); end
rd.dets = dets;
rd.notchBins = sys.dop.notchBins;
rd.beam = beam;
rd.t_s = t_s;
rd.az_deg = az_deg;
rd.dR_m = V.dR_m;
end
