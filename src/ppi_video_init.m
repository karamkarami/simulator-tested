function V = ppi_video_init(sys, rMax_m)
%PPI_VIDEO_INIT  Raw video of the PPI ('Video' display mode): map and pixel lookup.
%   V = PPI_VIDEO_INIT(sys, rMax_m) prepares the video map of the PPI out to
%   the range rMax_m: azimuth bins of V.azRes_deg x range bins of V.dR_m.
%   ppi_video_add.m writes the video of every PRI into the bins of its
%   azimuth, ppi_video_image.m draws the map as an image (colour = power).
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
rS = ws.range_m(1) + (0:numel(ws.rxIdx)-1) * dR;
keepS = rS < wl.range_m(1);
rL = wl.range_m(1) + (0:numel(wl.rxIdx)-1) * dR;
idx = [ws.rxIdx(keepS), wl.rxIdx];                            % DDC samples of the joined video
r   = [rS(keepS), rL];
V.dR_m = dR;
V.nR = max(1, ceil(rMax_m / dR));                             % range bins (at least one sample each)
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
