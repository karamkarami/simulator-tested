function xc = clutter_signal(clut, sys, az_deg, t_s)
%CLUTTER_SIGNAL  Complex IF clutter echo of one PRI, one column per elevation group.
%   xc = CLUTTER_SIGNAL(clut, sys, az_deg, t_s) returns the clutter echo of
%   both sub-pulses during the receive window [nAdc x G] for the antenna
%   azimuth az_deg at time t_s, column g holding the clutter of elevation
%   group g (clutter_terrain.m: flat ground, and every hill or mountain
%   region at its own elevation). The two-way azimuth pattern weights the
%   clutter cells (this also creates the scanning modulation); every cell
%   rotates with its internal-motion Doppler (plus a stationary part with
%   the exponential spectrum model, see clutter_init.m). Element i receives
%   real(xc * sys.clutterSteer(i, :).') (generate_rx_signal.m).
%   With the sub-pulses on different carriers (cfg.wave.<pulse>.fc_Hz) xc is
%   [nAdc x G x 2]: page p the echo of sub-pulse p (its IF, internal-motion
%   Doppler and power), for the steering vectors sys.wave.<pulse>.clutterSteer.

nAz = numel(clut.az_deg);
azRel = wrap180(clut.az_deg - az_deg);
sel   = find(abs(azRel) <= sys.ant.azSpan_deg);
gRow  = zeros(nAz, 1);
gRow(sel) = pattern_gain(sys.ant.az_deg, sys.ant.azGain_dB, azRel(sel));   % one-way power = two-way amplitude
% cells of the selected rows (one or two contiguous blocks of rows)
brk = [0; find(diff(sel) > 1); numel(sel)];
idx = cell(numel(brk) - 1, 1);
for b = 1:numel(brk) - 1
    a1 = sel(brk(b) + 1);  a2 = sel(brk(b + 1));
    idx{b} = (clut.rowStart(a1):clut.rowStart(a2 + 1) - 1).';
end
idx = vertcat(idx{:});
G = clut.nGroups;
nR = numel(clut.range_m);
same = ~isfield(clut, 'sameCarrier') || clut.sameCarrier;
nPages = 1 + ~same;
xc = zeros(sys.nAdc, G, nPages);
if isempty(idx), return; end
names = {'short', 'long'};
refl = cell(1, 2);
for i = 1:2
    p = clut.(names{i});
    fdScale = 1;  ampScale = 1;
    if isfield(p, 'fdScale'), fdScale = p.fdScale;  ampScale = p.ampScale; end
    if i == 2 && same && fdScale == clut.short.fdScale, refl{2} = refl{1}; continue; end
    w = gRow(clut.cellRow(idx)) .* clut.rho(idx) .* exp(1j*2*pi*fdScale*clut.fd_Hz(idx)*t_s);
    if ~isempty(clut.dc)                                 % stationary part (exponential model)
        w = w + gRow(clut.cellRow(idx)) .* clut.dc(idx);
    end
    refl{i} = ampScale * full(sparse(clut.cellRange(idx), clut.cellGroup(idx), w, nR, G));   % [nRange x G]
end

for g = 1:G
    for i = 1:2
        p = clut.(names{i});
        if isfield(p, 'enable') && ~p.enable, continue; end
        if ~any(refl{i}(:, g)), continue; end
        % Echo of each range cell starts at ADC sample n0 (0-based, may be negative)
        n0  = p.start0 + clut.delaySamples;
        use = n0 > -p.nTx & n0 < sys.nAdc;
        % Impulses shifted by nTx-1 so that partially received echoes are kept
        imp = zeros(p.nFft, 1);
        imp(n0(use) + p.nTx) = refl{i}(use, g);
        y = ifft(fft(imp) .* p.txFft);                     % convolution with the IF pulse
        pg = min(i, nPages);
        xc(:, g, pg) = xc(:, g, pg) + y(p.nTx : p.nTx + sys.nAdc - 1);
    end
end
end
