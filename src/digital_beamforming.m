function beams = digital_beamforming(bb, W)
%DIGITAL_BEAMFORMING  Elevation beamforming.
%   beams = DIGITAL_BEAMFORMING(bb, W) combines the element signals
%   bb [nSamples x nElements] with the DBF matrix W [nBeams x nElements] and
%   returns beams [nBeams x nSamples].
beams = W * bb.';
end
