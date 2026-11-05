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
xc = zeros(sys.nAdc, G);
if isempty(idx), return; end
w = gRow(clut.cellRow(idx)) .* clut.rho(idx) .* exp(1j*2*pi*clut.fd_Hz(idx)*t_s);
if ~isempty(clut.dc)                                     % stationary part (exponential model)
    w = w + gRow(clut.cellRow(idx)) .* clut.dc(idx);
end
refl = full(sparse(clut.cellRange(idx), clut.cellGroup(idx), w, nR, G));   % [nRange x G]

names = {'short', 'long'};
for g = 1:G
    if ~any(refl(:, g)), continue; end
    for i = 1:2
        p = clut.(names{i});
        % Echo of each range cell starts at ADC sample n0 (0-based, may be negative)
        n0  = p.start0 + clut.delaySamples;
        use = n0 > -p.nTx & n0 < sys.nAdc;
        % Impulses shifted by nTx-1 so that partially received echoes are kept
        imp = zeros(p.nFft, 1);
        imp(n0(use) + p.nTx) = refl(use, g);
        y = ifft(fft(imp) .* p.txFft);                     % convolution with the IF pulse
        xc(:, g) = xc(:, g) + y(p.nTx : p.nTx + sys.nAdc - 1);
    end
end
end
