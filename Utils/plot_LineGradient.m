function hline = plot_LineGradient(x, y, colorData, varargin)

p = inputParser();
p.addRequired('x')
p.addRequired('y')
p.addRequired('colorData')
p.addOptional('linestyle', '-', @isstr)
p.addParameter('colorShiftSign', -1)
p.KeepUnmatched = true;
parse(p, x, y, colorData, varargin{:})

offsetsign = p.Results.colorShiftSign;

hline = plot(x, y, p.Results.linestyle, p.Unmatched);
% modified jet-colormap
n = length(x);
cd = [ones(n, 1), linspace(.7, 1, n)', linspace(.7, 1, n)'];
cd = cd .* rgb2hsv(colorData);
cd(:, 1) = mod(cd(:, 1) + offsetsign * linspace(.25, 0, n)', 1);
cd = uint8([hsv2rgb(cd), ones(n, 1)]*255).';
drawnow
set(hline.Edge, 'ColorBinding', 'interpolated', 'ColorData', cd)

hMarkers = hline.MarkerHandle;
set(hMarkers,'FaceColorBinding','interpolated', 'FaceColorData', cd)
set(hMarkers,'EdgeColorBinding','interpolated', 'EdgeColorData', cd)
end