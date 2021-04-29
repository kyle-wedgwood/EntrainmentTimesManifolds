% Plot time simulations near the saddle and stable limit cycles
% 30.10.2020
% Kyle Wedgwood

close all; clear; clc;

N = 12;

p = setDefaultParameters;
p = setDayLengthParameters(p, N);
p.tShift = 0.0;

folder = '/Users/kcaw201/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds/';
stable_orbit = load([folder 'Orbit_I_50_UZ3.dat']);
saddle_orbit = load([folder 'Orbit_I_50_UZ4.dat']);
unstable_orbit = load([folder 'Orbit_I_50_UZ5.dat']);

stab_max = max(stable_orbit(:,1));
stab_min = min(stable_orbit(:,1));
sad_max = max(saddle_orbit(:,1));
sad_min = min(saddle_orbit(:,1));

no_cycles = 1;
tfinal = no_cycles*24;

odeopts = odeset('AbsTol', 1e-10, 'RelTol', 1e-8);

[t_stable,y_stable] = ode45(@Kronauer2SingleFS, [0 tfinal], ...
                            stable_orbit(1,:)', odeopts, p);
[t_saddle,y_saddle] = ode45(@Kronauer2SingleFS, [0 tfinal], ...
                            saddle_orbit(1,:)'+0.1*randn(2,1), odeopts, p);
[t_saddle,y_saddle] = ode45(@Kronauer2SingleFS, [0 tfinal], ...
                            stable_orbit(2700,:), odeopts, p);
[t_unstable,y_unstable] = ode45(@Kronauer2SingleFS, [0 tfinal], ...
                            unstable_orbit(1,:)'+0.2*randn(2,1), odeopts, p);

[fig,ax] = setupFigure;

light_ind = logical(size(t_saddle));
for n = 0:no_cycles
  ind = (n*24 < t_saddle) & (t_saddle < n*24 + N);
  light_ind(ind) = true;
end

% plot(ax, t_stable, y_stable, 'Linewidth', line_width);
line(ax, [0, no_cycles], [stab_min,stab_min], 'Linewidth', 0.5*line_width, ...
     'Color', 'black', 'Linestyle', '--');
line(ax, [0, no_cycles], [stab_max,stab_max], 'Linewidth', 0.5*line_width, ...
     'Color', 'black', 'Linestyle', '--');
line(ax, [0, no_cycles], [sad_min,sad_min], 'Linewidth', 0.5*line_width, ...
     'Color', 'black', 'Linestyle', ':');
line(ax, [0, no_cycles], [sad_max,sad_max], 'Linewidth', 0.5*line_width, ...
     'Color', 'black', 'Linestyle', ':');

temp = y_saddle;
temp(~light_ind) = nan;
plot(ax, t_saddle/24, temp(:,1), 'Linewidth', line_width, ...
     'Color', grey);
temp = y_saddle;
temp(light_ind) = nan;
plot(ax, t_saddle/24, temp(:,1), 'Linewidth', line_width, ...
     'Color', 'black');

% plot(ax, t_unstable, y_unstable, 'Linewidth', line_width);

xlab = xlabel(ax, 'Time (d)');
ylab = ylabel(ax, '$A$', 'Rotation', 0);
ylab.Position(1) = ylab.Position(1) - 0.4;
ylab.Position(2) = ylab.Position(2) + 0.4;
set(ax, 'YLim', [-2,2]);
set(ax, 'XTick', 0:2:no_cycles, 'YTick', -2:2);

[fig,ax] = postProcessFig(fig, ax, 3, 1, 'time_series_I_50');

%% Next do situation for I = 1000
p.I = 1000;
stable_orbit = load([folder 'Orbit_I_1000_UZ1.dat']);
stab_max = max(stable_orbit(:,1));
stab_min = min(stable_orbit(:,1));

[t,y] = ode45(@Kronauer2SingleFS, [0 tfinal], ...
              stable_orbit(2700,:), odeopts, p);

[fig,ax] = setupFigure;

light_ind = logical(size(t));
for n = 0:no_cycles
  ind = (n*24 < t) & (t < n*24 + N);
  light_ind(ind) = true;
end

line(ax, [0, no_cycles], [stab_min,stab_min], 'Linewidth', 0.5*line_width, ...
     'Color', 'black', 'Linestyle', '--');
line(ax, [0, no_cycles], [stab_max,stab_max], 'Linewidth', 0.5*line_width, ...
     'Color', 'black', 'Linestyle', '--');

temp = y;
temp(~light_ind) = nan;
plot(ax, t/24, temp(:,1), 'Linewidth', line_width, ...
     'Color', grey);
temp = y;
temp(light_ind) = nan;
plot(ax, t/24, temp(:,1), 'Linewidth', line_width, ...
     'Color', 'black');