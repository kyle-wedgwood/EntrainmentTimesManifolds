% Plot points along manifold sequentially
% 29.10.2020
% Kyle Wedgwood

function plotOneManSeq(filename)

  y = load(filename);

  figure;
  hold on;
  g = plot(y(1:2,1), y(1:2,2), 'Linewidth', 2, 'Color', 'b');
  h = plot(y(2:3,1), y(2:3,2), 'Linewidth', 2, 'Color', 'r');

  for i = 2:size(y, 1)
    xdata = [get(g, 'XData'), get(h, 'XData')];
    ydata = [get(g, 'YData'), get(h, 'YData')];
    set(g, 'XData', xdata, 'YData', ydata);
    set(h, 'XData', y(i:i+1,1), 'YData', y(i:i+1,2));
    title(sprintf('i: %d', i));
    drawnow; pause;
  end

end