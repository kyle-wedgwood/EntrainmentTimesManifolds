% Plot Poincare map of reentrainment trajectories
% 4.11.2020
% Kyle Wedgwood

close all; clear; clc;

p = setDefaultParameters;
p = setDayLengthParameters(p, 12);
p.tShift = 0.0;

orbit_folder = '~/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds/';
orbit_filename = [orbit_folder, 'Orbit_I_50_UZ3.dat'];
man_name = 'Manifold_I_50_UZ*';
man_ind = 5;
files = dir([orbit_folder man_name]);
orbit = load(orbit_filename);
man = load([orbit_folder files(man_ind).name]);

% Now do simulations
odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
no_iterates = 80;
z = zeros(no_iterates+1, 2);
cm = lines(10);

% Find start point
res = zeros(size(orbit, 1), 1);
for i = 1:size(orbit, 1)
  res(i) = norm(man(1,:) - orbit(i,:));
end

[~,ind] = min(res);

% First do orthodromic
fig_filename = 'reentrain_orthodromic_I_50';
y0 = orbit(ind,:)';
z(1,:) = y0;

[fig,ax] = setupFigure();

for n = 1:no_iterates
  [t,y] = ode45(@Kronauer2SingleFS, [0 24], y0, odeopts, p);
  y0 = y(end,:);
  z(n+1,:) = y0;
end

% Plot
light_ind = 1:(length(orbit)-1)/2;
dark_ind  = light_ind(end)+1:length(orbit);
plot(ax, orbit(light_ind,1), orbit(light_ind,2), 'Linewidth', line_width, ...
     'Color', grey, 'HandleVisibility', 'off');
plot(ax, orbit(dark_ind,1), orbit(dark_ind,2), 'Linewidth', line_width, ...
     'Color', 'black', 'HandleVisibility', 'off');
   
plot(ax, z(:,1), z(:,2), 'Marker', '.', 'Markersize', marker_size, ...
       'Color', [0.75 0.75 0]);
     
% Add marker for point of zero phase on stable limit cycle
plot(ax, orbit(1,1), orbit(1,2), 'Marker', '.', 'Markersize', 2*marker_size, ...
     'Color', 'black');
     
xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);
axis square
[fig,ax] = postProcessFig(fig, ax, 3, 1, 2.4,fig_filename, 'B', ax_pos);

% Next do anti-dromic
fig_filename = 'reentrain_antidromic_I_50';
y0 = orbit(ind-50,:)';
z(1,:) = y0;

for n = 1:no_iterates
  [t,y] = ode45(@Kronauer2SingleFS, [0 24], y0, odeopts, p);
  y0 = y(end,:);
  z(n+1,:) = y0;
end

[fig,ax] = setupFigure();

% Plot
light_ind = 1:(length(orbit)-1)/2;
dark_ind  = light_ind(end)+1:length(orbit);
plot(ax, orbit(light_ind,1), orbit(light_ind,2), 'Linewidth', line_width, ...
     'Color', grey, 'HandleVisibility', 'off');
plot(ax, orbit(dark_ind,1), orbit(dark_ind,2), 'Linewidth', line_width, ...
     'Color', 'black', 'HandleVisibility', 'off');

plot(ax, z(:,1), z(:,2), 'Marker', '.', 'Markersize', marker_size, ...
       'Color', [0 0.5 0]);
     
% Add marker for point of zero phase on stable limit cycle
plot(ax, orbit(1,1), orbit(1,2), 'Marker', '.', 'Markersize', 2*marker_size, ...
     'Color', 'black');
     
xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);
axis square
[fig,ax] = postProcessFig(fig, ax, 3, 1, 2.4,fig_filename, 'A', ax_pos);

% Repeat for I = 1000
orbit_filename = [orbit_folder, 'Orbit_I_1000_UZ1.dat'];
man_name = 'Manifold_I_1000_UZ*';
man_ind = 2;
files = dir([orbit_folder man_name]);
orbit = load(orbit_filename);
man = load([orbit_folder files(man_ind).name]);

% Find shortcut
p.I = 1000;
fig_filename = 'reentrain_shortcut_I_1000';
y0 = orbit(ind-215,:)'; % 340
z(1,:) = y0;

for n = 1:no_iterates
  [t,y] = ode45(@Kronauer2SingleFS, [0 24], y0, odeopts, p);
  y0 = y(end,:);
  z(n+1,:) = y0;
end

[fig,ax] = setupFigure();

% Plot
light_ind = 1:(length(orbit)-1)/2;
dark_ind  = light_ind(end)+1:length(orbit);
plot(ax, orbit(light_ind,1), orbit(light_ind,2), 'Linewidth', line_width, ...
     'Color', grey, 'HandleVisibility', 'off');
plot(ax, orbit(dark_ind,1), orbit(dark_ind,2), 'Linewidth', line_width, ...
     'Color', 'black', 'HandleVisibility', 'off');

plot(ax, z(:,1), z(:,2), 'Marker', '.', 'Markersize', marker_size, ...
       'Color', cm(5,:));
     
% Add marker for point of zero phase on stable limit cycle
plot(ax, orbit(1,1), orbit(1,2), 'Marker', '.', 'Markersize', 2*marker_size, ...
     'Color', 'black');

xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);
axis square
[fig,ax] = postProcessFig(fig, ax, 3, 1, 2.4,fig_filename, 'C', ax_pos);