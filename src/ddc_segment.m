function seg = ddc_segment(bb, wv)
%DDC_SEGMENT  Input of the pulse compression of a sub-pulse.
%   seg = DDC_SEGMENT(bb, wv) returns the DDC samples wv.ddcIdx of bb
%   [nDdc x nChannels] (seg [numel(wv.ddcIdx) x nChannels]). The samples
%   outside the receive interval wv.rxIdx of the sub-pulse (transmission,
%   the other sub-pulse, after the end of the record) are zero: a
%   compression filter longer than the pulse (mismatch filter) then gives
%   the same range cells as the matched filter (setup_radar.m).
in = wv.ddcIdx >= wv.rxIdx(1) & wv.ddcIdx <= wv.rxIdx(end);
if all(in)
    seg = bb(wv.ddcIdx, :);
else
    seg = zeros(numel(wv.ddcIdx), size(bb, 2), class(bb));
    seg(in, :) = bb(wv.ddcIdx(in), :);
end
end
