function [v, u] = steering_vector(sys, azRel_deg, el_deg)
%STEERING_VECTOR  Phase of a plane wave on the elevation elements of the tilted array.
%   [v, u] = STEERING_VECTOR(sys, azRel_deg, el_deg) returns the nElements x 1
%   steering vector for a source at azimuth azRel_deg (relative to the antenna
%   boresight) and elevation el_deg (above the horizon). u is the sine of the
%   angle between the source and the array broadside, measured along the
%   element column (u = sin(el - tilt) when azRel = 0).

tilt = sys.ant.tilt_deg;
u = -sind(tilt)*cosd(el_deg)*cosd(azRel_deg) + cosd(tilt)*sind(el_deg);
n = (0:sys.nCh-1)';
v = exp(1j*2*pi*sys.ant.elementSpacing_m/sys.lambda * n * u);
end
