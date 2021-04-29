% Plot (A,C) solutions

close all; clear; clc;

p = setDefaultParameters();
p = setDayLengthParameters(p, 12);
p.I = 1000;
p.tShift = 0.0;
p.taux = 26.5;

orbit_folder = '~/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds/';

% Prepare figure
[fig,ax] = setupFigure();

stable_orbit_filename = [orbit_folder, 'Orbit_I_1000_UZ1.dat'];
stable_orbit = load(stable_orbit_filename);

y0 = stable_orbit(1,:)';
odeopts = odeset('AbsTol', 1e-12, 'RelTol', 1e-10);
no_cycles = 1e6;

for i = 1:no_cycles
  [t,y] = ode45(@Kronauer2SingleFS, [0,24], y0, odeopts, p);
  y0 = y(end,:)';
  plot(ax, y0(1), y0(2), 'Marker', '.', 'Color', 'black', 'Markersize', ...
       marker_size);
  drawnow;
end

xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);
