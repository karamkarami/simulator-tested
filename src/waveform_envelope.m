function u = waveform_envelope(wv, t)
%WAVEFORM_ENVELOPE  Complex envelope of a sub-pulse.
%   u = WAVEFORM_ENVELOPE(wv, t) evaluates the envelope at the times t [s]
%   (t = 0 is the start of the pulse). The envelope is zero outside the pulse.
%   wv.type = 'lfm'       : linear FM chirp from -B/2 to +B/2
%   wv.type = 'phasecode' : the chips wv.code (from cfg.files.phaseCode, set by
%                           setup_radar.m), each wv.chip_s seconds long

switch lower(wv.type)
    case 'lfm'
        k = wv.bandwidth_Hz / wv.width_s;               % chirp rate [Hz/s]
        u = exp(1j*pi*k*(t - wv.width_s/2).^2);
    case 'phasecode'
        nChips = numel(wv.code);
        m = min(max(floor(t / wv.chip_s), 0), nChips-1); % chip index (0-based)
        u = reshape(wv.code(m + 1), size(t));
    otherwise
        error('Unknown waveform type "%s"', wv.type);
end
u(t < 0 | t >= wv.width_s) = 0;
end
