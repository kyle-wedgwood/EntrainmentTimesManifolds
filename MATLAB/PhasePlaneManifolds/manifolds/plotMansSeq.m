% Plot manifolds to test which is which
% 29.10.2020
% Kyle Wedgwood

close all; clear; clc;

files = dir('*.dat');

fig = figure;
ax = axes(fig);

h = plot(0, 0, 'Linewidth', 2);

set(ax, 'Xlim', [-2, 2]);
set(ax, 'YLim', [-2, 2]);
xlabel(ax, 'A');
ylabel(ax, 'C', 'Rotation', 0);

for i = 1:length(files)
  data = load(files(i).name);
  set(h, 'XData', data(:,1), 'YData', data(:,2));
  pause;
end

