function dets = eclipse_blanking(dets, short, sys)
%ECLIPSE_BLANKING  Drop long-pulse detections that come from the eclipsed long-pulse echo of a nearer target.
%   dets = ECLIPSE_BLANKING(dets, short, sys) prunes the long-pulse detections
%   dets of a CPI update (prune_peaks.m) with the short-pulse detections short
%   of the same update. A target nearer than the blind range of the long
%   pulse (range R) sends back the end of the long pulse while the long-pulse
%   receive window has already started: the pulse compression spreads this
%   partial echo over the long-pulse range cells up to R + sys.eclipseExtent_m
%   (the pulse length plus the earlier peak of a mismatch filter and the DDC
%   filters), in the Doppler bin of the target, in its beam and, through the
%   elevation sidelobes of the DBF, in the other beams. With a strong target
%   (e.g. at 170 MHz) its peaks reach the CFAR threshold. A long-pulse
%   detection is dropped when a short-pulse detection
%     - is in the same Doppler bin (+-1, circular; any beam),
%     - is at most sys.eclipseExtent_m nearer,
%     - is more than sys.cfg.prune.sidelobeMargin_dB stronger,
%   the rules of the sidelobe blanking in prune_peaks.m. short = [] (no
%   short-pulse update) keeps dets as they are.

if isempty(short) || isempty(short.beam) || isempty(dets.beam), return; end
nD = sys.dop.nfft;
ratio = 10^(sys.cfg.prune.sidelobeMargin_dB/10);
rS = sys.wave.short.range_m(short.rIdx);  rS = rS(:);
rL = sys.wave.long.range_m(dets.rIdx);
pS = double(short.power(:));
keep = true(size(dets.beam));
for i = 1:numel(dets.beam)
    dDop = abs(short.dIdx(:) - dets.dIdx(i));
    hit = min(dDop, nD - dDop) <= 1 & rL(i) <= rS + sys.eclipseExtent_m & pS > ratio * double(dets.power(i));
    keep(i) = ~any(hit);
end
dets.beam  = dets.beam(keep);
dets.rIdx  = dets.rIdx(keep);
dets.dIdx  = dets.dIdx(keep);
dets.power = dets.power(keep);
dets.noise = dets.noise(keep);
end
