function a = wrap180(a)
%WRAP180  Wrap angles in degrees to the interval [-180, 180).
a = mod(a + 180, 360) - 180;
end
