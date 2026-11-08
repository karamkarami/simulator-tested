function cmap = ppi_video_colormap(n)
%PPI_VIDEO_COLORMAP  Colours of the PPI video (low to high power): black, blue,
%   cyan, green, yellow, red, white. The C receiver's PPI uses the same.
x = [0 0.15 0.35 0.55 0.75 0.9 1];
c = [0 0 0; 0 0 0.55; 0 0.6 1; 0 1 0.3; 1 1 0; 1 0.25 0; 1 1 1];
t = linspace(0, 1, n).';
cmap = [interp1(x, c(:, 1), t), interp1(x, c(:, 2), t), interp1(x, c(:, 3), t)];
end
