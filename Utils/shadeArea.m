function shadeArea(x, curve1, curve2)

if ~isrow(x)
    x = x';
end

if ~isrow(curve1)
    curve1 = curve1';
end

if ~isrow(curve2)
    curve2 = curve2';
end

mainLine = plot(x, curve1, 'LineWidth', 2);
hold on;
plot(x, curve2, 'LineWidth', 2, 'Color', mainLine.Color);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fillColor = min(1, mainLine.Color + .2*ones(1,3));
fill(x2, inBetween, fillColor);

end

