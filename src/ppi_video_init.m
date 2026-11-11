function V = ppi_video_init(sys, rMax_m)
%PPI_VIDEO_INIT  Raw data map of the PPI (display modes Video, DBF out, MF out): map and pixel lookup.
%   V = PPI_VIDEO_INIT(sys, rMax_m) prepares the video map of the PPI out to
%   the range rMax_m: azimuth bins of V.azRes_deg x range bins of V.dR_m.
%   ppi_video_add.m writes the video of every PRI into the bins of its
%   azimuth, ppi_video_image.m draws the map as an image (colour = power).
%   The same map holds the DBF output or the pulse compression output of the
%   selected beam (V.kind 'dbf', 'mf'; ppi_profile_row.m, ppi_video_add_row.m);
%   V.dbf{p}, V.mf{p} (p = 1 short, 2 long) give the samples and range bins of
%   both outputs, joined like the video; ppi_rd_map.m puts the range-Doppler
%   map of the selected beam into the same range bins.
%
%   Video of a PRI: the DDC output of all receive elements averaged
%   non-coherently (mean of I^2 + Q^2 over the elements, before the
%   beamforming and the pulse compression). The receive intervals of the short
%   and the long sub-pulse are cut out of it and joined: the short-pulse
%   samples nearer than the start of the long-pulse range (the handover),
%   then the long-pulse samples. Sample k of an interval is the range of an
%   echo whose leading edge arrives at it: range_m(1) + (k-1)*c/(2 fsOut).
c = sys.cfg.const.c;
dR = c / (2 * sys.fsOut);                                     % range per DDC sample [m]
ws = sys.wave.short;  wl = sys.wave.long;
handover = sys.handover_m;                                    % (Inf without the long pulse)
rS = first_range(ws) + (0:numel(ws.rxIdx)-1) * dR;
keepS = rS < handover;
rL = first_range(wl) + (0:numel(wl.rxIdx)-1) * dR;
idx = [ws.rxIdx(keepS), wl.rxIdx];                            % DDC samples of the joined video
r   = [rS(keepS), rL];
V.dR_m = dR;
V.nR = max(1, ceil(rMax_m / dR));                             % range bins (at least one sample each)
V.kind = 'video';                                             % what the map holds: 'video', 'dbf', 'mf'
V.beam = 1;                                                   % beam of 'dbf' / 'mf'
% DBF output: the pulse compression input of a sub-pulse (wv.ddcIdx, ddc_segment.m) inside its
% receive interval; MF output: its range cells. Short pulse: nearer than the long pulse's first cell.
waves = {ws, wl};                                             % (own variable names: r, ok follow below)
for p = 1:2
    wv = waves{p};
    if ~wv.enable                                             % sub-pulse switched off: nothing
        V.dbf{p} = struct('idx', zeros(1, 0), 'bin', zeros(1, 0));
        V.mf{p} = V.dbf{p};
        continue
    end
    jj = find(ismember(wv.ddcIdx, wv.rxIdx));                 % segment samples in the receive interval
    rr = wv.range_m(1) + (wv.ddcIdx(jj) - wv.rxIdx(1)) * dR;
    keep = rr < handover | p == 2;
    bb = floor(rr(keep) / dR) + 1;  inside = bb >= 1 & bb <= V.nR;
    jj = jj(keep);
    V.dbf{p} = struct('idx', jj(inside), 'bin', bb(inside));
    rr = wv.range_m(:).';
    keep = rr < handover | p == 2;
    bb = floor(rr(keep) / dR) + 1;  inside = bb >= 1 & bb <= V.nR;
    kk = find(keep);
    V.mf{p} = struct('idx', kk(inside), 'bin', bb(inside));
end
bin = floor(r / dR) + 1;
ok = bin >= 1 & bin <= V.nR;
V.sampleIdx = idx(ok);
V.sampleBin = bin(ok);
V.binCount = accumarray(V.sampleBin(:), 1, [V.nR 1]).';       % samples per range bin (0: none)
V.azRes_deg = 0.5;
V.nAz = round(360 / V.azRes_deg);
V.map = zeros(V.nAz, V.nR, 'single');                         % mean power [LSB^2] (0 = no data)
V.time = -Inf(V.nAz, 1);                                      % time of the newest PRI in every bin
V.sum = zeros(1, V.nR);  V.n = 0;  V.azBin = 0;                % PRIs of the current azimuth bin
V.noise_dB = NaN;                                             % noise floor (running 10th percentile of the video)
V.azNoise = nan(V.nAz, 1);                                    % its estimate in every azimuth bin (ppi_video_add_row.m)
x = gammaincinv(0.1, sys.nCh);                                % its offset: the 10th percentile of the mean of
V.p10Video_dB = 10*log10(sys.nCh / x);                        % nCh exponential powers is below the mean (1: 9.77
V.p10Sample_dB = 10*log10(1 / gammaincinv(0.1, 1));           % dB, 48: 0.86 dB); one sample (DBF, MF out): 9.77 dB
V.span_dB = 50;                                               % colour scale: noise - 3 dB ... + span
% Image of the PPI: nPix x nPix pixels over [-rMax, rMax] (km in the axes), lookup into the map
V.nPix = 601;
V.rMax_km = rMax_m / 1e3;
g = linspace(-rMax_m, rMax_m, V.nPix);
[X, Y] = meshgrid(g, g);                                      % row 1 = South (image YData from -rMax)
R = hypot(X, Y);
A = mod(atan2(X, Y) * 180/pi, 360);
iAz = min(floor(A / V.azRes_deg) + 1, V.nAz);
% A pixel is larger than a range bin: it shows the largest of the V.binsPerPix range bins
% under it (a small echo does not fall between the pixels)
V.binsPerPix = max(1, ceil(2*rMax_m / (V.nPix - 1) / dR));
V.nRPix = ceil(V.nR / V.binsPerPix);
iR  = floor(R / (dR * V.binsPerPix)) + 1;
V.pixInside = R <= rMax_m & iR <= V.nRPix;
V.pixAz = iAz(V.pixInside);
V.pixLin = iAz(V.pixInside) + V.nAz * (iR(V.pixInside) - 1);  % linear index into the reduced map
V.cmap = ppi_video_colormap(256);
end

function r = first_range(wv)
% Range of the first receive sample of a sub-pulse (none when it is switched off).
r = 0;
if ~isempty(wv.range_m), r = wv.range_m(1); end
end
