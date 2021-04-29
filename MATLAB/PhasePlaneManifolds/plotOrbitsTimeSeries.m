% Plot manifolds as time series
% 30.10.2020
% Kyle Wedgwood

close all; clear; clc;

N = 12;

p = setDefaultParameters;
p = setDayLengthParameters(p, N);
p.tShift = 0.0;

folder = '~/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds/';
stable_orbit = load([folder 'Orbit_I_50_UZ3.dat']);
saddle_orbit = load([folder 'Orbit_I_50_UZ4.dat']);
unstable_orbit = load([folder 'Orbit_I_50_UZ5.dat']);
t = linspace(0, 24, length(stable_orbit));

stab_max = max(stable_orbit(:,1));
stab_min = min(stable_orbit(:,1));
sad_max = max(saddle_orbit(:,1));
sad_min = min(saddle_orbit(:,1));

[fig,ax] = setupFigure;

ind_light = 1:(length(stable_orbit)-1)/2;
ind_dark = ind_light(end)+1:length(stable_orbit);

h = fill(ax, [12 24 24 12], [-2,-2,2,2], 0.7*ones( 1, 3), ...'
         'FaceAlpha', 0.4, 'EdgeColor', 'none');

plot(ax, t(ind_light), stable_orbit(ind_light,1), ...
     'Linewidth', line_width, 'Color', grey);
plot(ax, t(ind_dark), stable_orbit(ind_dark,1), ...
     'Linewidth', line_width, 'Color', 'black');
plot(ax, t(ind_light), saddle_orbit(ind_light,1), ...
     'Linewidth', line_width, 'Color', grey, 'Linestyle', STYs{2});
plot(ax, t(ind_dark), saddle_orbit(ind_dark,1), ...
     'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{2});
plot(ax, t(ind_light), unstable_orbit(ind_light,1), ...
     'Linewidth', line_width, 'Color', grey, 'Linestyle', STYs{3});
plot(ax, t(ind_dark), unstable_orbit(ind_dark,1), ...
     'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{3});

xlab = xlabel(ax, 'Time (h)');
ylab = ylabel(ax, '$A$', 'Rotation', 0);
ylab.Position(1) = ylab.Position(1) - 0.4;
ylab.Position(2) = ylab.Position(2) + 0.4;
set(ax, 'XLim', [0 24], 'YLim', [-2,2]);
set(ax, 'XTick', 0:4:24, 'YTick', -2:2);

[fig,ax] = postProcessFig(fig, ax, 3, 1, 'time_series_I_50', 'C', ax_pos);

%% Plot as phase plane
[fig,ax] = setupFigure;

plot(ax, stable_orbit(ind_light,1), stable_orbit(ind_light,2), ...
     'Linewidth', line_width, 'Color', grey);
plot(ax, stable_orbit(ind_dark,1), stable_orbit(ind_dark,2), ...
     'Linewidth', line_width, 'Color', 'black');
plot(ax, saddle_orbit(ind_light,1), saddle_orbit(ind_light,2), ...
     'Linewidth', line_width, 'Color', grey, 'Linestyle', STYs{2});
plot(ax, saddle_orbit(ind_dark,1), saddle_orbit(ind_dark,2), ...
     'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{2});
plot(ax, unstable_orbit(ind_light,1), unstable_orbit(ind_light,2), ...
     'Linewidth', line_width, 'Color', grey, 'Linestyle', STYs{3});
plot(ax, unstable_orbit(ind_dark,1), unstable_orbit(ind_dark,2), ...
     'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{3});
   
pt = plot(ax, stable_orbit(1,1), stable_orbit(1,2), ...
         'Marker', 'o', 'Markersize', marker_size, 'Color', 'black', ...
         'MarkerFaceColor', 'black');pause(1);
for sz = 1:-0.3:0.2
  pt = plot(ax, saddle_orbit(1,1), saddle_orbit(1,2), ...
         'Marker', 'o', 'Markersize', sz*marker_size, 'Color', 'black');
end
pt = plot(ax, unstable_orbit(1,1), unstable_orbit(1,2), ...
          'Marker', 'o', 'Markersize', marker_size, 'Color', 'black', ...
          'MarkerFaceColor', 'white');

xlab = xlabel(ax, '$A$');
xlab.Position(1) = xlab.Position(1) + 1.0;
xlab.Position(2) = xlab.Position(2) - 0.6;
ylab = ylabel(ax, '$C$', 'Rotation', 0);
ylab.Position(1) = ylab.Position(1) - 0.8;
ylab.Position(2) = ylab.Position(2) + 0.9;
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2:2, 'YTick', -2:2:2);
axis square
[fig,ax] = postProcessFig(fig, ax, 3, 1, 'orbits_pp_I_50', 'E', ax_pos);

%% Next do situation for I = 1000
p.I = 130;
stable_orbit = load([folder 'Orbit_I130pt00_taux_24.2_UZ2.dat']);
stab_max = max(stable_orbit(:,1));
stab_min = min(stable_orbit(:,1));

[fig,ax] = setupFigure;

h = fill(ax, [12 24 24 12], [-2,-2,2,2], 0.7*ones( 1, 3), ...'
         'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(ax, t(ind_light), stable_orbit(ind_light,1), ...
     'Linewidth', line_width, 'Color', grey);
plot(ax, t(ind_dark), stable_orbit(ind_dark,1), ...
     'Linewidth', line_width, 'Color', 'black');
   
xlab = xlabel(ax, 'Time (h)');
ylab = ylabel(ax, '$A$', 'Rotation', 0);
ylab.Position(1) = ylab.Position(1) - 0.4;
ylab.Position(2) = ylab.Position(2) + 0.4;
set(ax, 'Xlim', [0 24], 'YLim', [-2,2]);
set(ax, 'XTick', 0:4:24, 'YTick', -2:2);

[fig,ax] = postProcessFig(fig, ax, 3, 1, 'time_series_I_130', 'B', ax_pos);

%% Plot as phase plane
[fig,ax] = setupFigure;

plot(ax, stable_orbit(ind_light,1), stable_orbit(ind_light,2), ...
     'Linewidth', line_width, 'Color', grey);
plot(ax, stable_orbit(ind_dark,1), stable_orbit(ind_dark,2), ...
     'Linewidth', line_width, 'Color', 'black');
   
pt = plot(ax, stable_orbit(1,1), stable_orbit(1,2), ...
         'Marker', 'o', 'Markersize', marker_size, 'Color', 'black', ...
         'MarkerFaceColor', 'black');pause(1);

xlab = xlabel(ax, '$A$');
xlab.Position(1) = xlab.Position(1) + 1.0;
xlab.Position(2) = xlab.Position(2) - 0.6;
ylab = ylabel(ax, '$C$', 'Rotation', 0);
ylab.Position(1) = ylab.Position(1) - 0.8;
ylab.Position(2) = ylab.Position(2) + 0.9;
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2:2, 'YTick', -2:2:2);

% [fig,ax] = postProcessFig(fig, ax, 3, 1, 'orbits_pp_I_130');