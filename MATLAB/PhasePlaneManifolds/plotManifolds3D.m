% Integrate manifolds at one tShift value to fill out entire phase space
% 22.5.2019
% Kyle Wedgwood

close all; clear; clc;

data_folder = ...
  '/Users/kcaw201/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds';
addpath(data_folder); 

orbit = load('orbits/Orbit_lux_50_UZ2.dat');
saddle_orbit = load('orbits/Orbit_lux_50_UZ3.dat');
unstable_orbit = load('orbits/Orbit_lux_50_UZ4.dat');

% % Load stable point manifolds
% stable_man = [];
% 
% for i = [1,4]
%   stable_man_temp = load(sprintf('orbits/Manifold_lux_50_UZ2_%d.dat', i));
%   dm = diff(stable_man_temp, 1);
%   ind = find(dm(:,1) ~= 0);
%   stable_man_temp = stable_man_temp(ind,:);
%   stable_man = [flipud(stable_man); stable_man_temp];
% end
% 
% plot(stable_man(:,1), stable_man(:,2));
% 
% % Load saddle point manifolds
% saddle_man = load('Manifold_lux_50_UZ3_4.dat');
% saddle_man_temp = load('Manifold_lux_50_UZ3_5.dat');
% dm = diff(saddle_man_temp, 1);
% stop_ind = find(abs(dm(:,1)) > 0.4);
% saddle_man = [flipud(saddle_man); saddle_man_temp(1:stop_ind,:)];
% 
% plot(saddle_man(:,1), saddle_man(:,2));

% Load precomputed manifolds
stable_man = load('manifolds/stable_stable_I_50.dat');
saddle_man = load('manifolds/stable_saddle_I_50.dat');
p = setDefaultParameters();
p = setDayLengthParameters(p, 12);

% Plot options
dark_blue = [0.2549019607843137, 0.4117647058823529, 0.882352941176470];
light_blue = [0.6901960784313725, 0.7686274509803922, 0.8705882352941177];

p.domain_radius = 2.3;
plot_radius = 1.8;

% Interpolate
no_pts = 2000;

% Now integrate manifolds
[T_saddle_forward,A_saddle_forward,C_saddle_forward] = ...
   integrateManifoldForwards(saddle_man, p, no_pts);
[T_stable_forward,A_stable_forward,C_stable_forward] = ...
   integrateManifoldForwards(stable_man, p, no_pts);
[T_saddle_back,A_saddle_back,C_saddle_back] = ...
   integrateManifoldBackwards(saddle_man, p, no_pts);
[T_stable_back,A_stable_back,C_stable_back] = ...
   integrateManifoldBackwards(stable_man, p, no_pts);

%% Draw line around points
no_pts = size(T_saddle_back,1);
saddle_boundary = [];
T_saddle_boundary = [];
for i = 1:no_pts
  normed_dist = sqrt(A_saddle_back(i,:).^2 + C_saddle_back(i,:).^2);
  ind = find((normed_dist(1:end-1) < plot_radius).*(normed_dist(2:end) > plot_radius));
  if ~isempty(ind)
    ind = [ind:ind+1];
    A_bound = interp1(normed_dist(ind), ...
                      A_saddle_back(i,ind), plot_radius);
    C_bound = interp1(normed_dist(ind), ...
                      C_saddle_back(i,ind), plot_radius);                
    saddle_boundary(end+1,:) = [A_bound; C_bound];
    T_saddle_boundary = [T_saddle_boundary; T_saddle_back(i,1)];
  end
end

no_pts = size(T_stable_back,1);
stable_boundary = [];
T_stable_boundary = [];
for i = 1:no_pts
  normed_dist = sqrt(A_stable_back(i,:).^2 + C_stable_back(i,:).^2);
  ind = find((normed_dist(2:end) < plot_radius).*(normed_dist(1:end-1) > plot_radius));
  if ~isempty(ind)
    ind = [ind:ind+1];
    A_bound = interp1(normed_dist(ind), ...
                      A_stable_back(i,ind), plot_radius);
    C_bound = interp1(normed_dist(ind), ...
                      C_stable_back(i,ind), plot_radius);                
    stable_boundary(end+1,:) = [A_bound; C_bound];
    T_stable_boundary = [T_stable_boundary; T_stable_back(i,1)];
  end
end
    
%% Set circular limits
ind_saddle = (A_saddle_back.^2 + C_saddle_back.^2) > plot_radius^2;
A_saddle_back(ind_saddle(:)) = nan;
C_saddle_back(ind_saddle(:)) = nan;

ind_stable = (A_stable_back.^2 + C_stable_back.^2) > plot_radius^2;
A_stable_back(ind_stable(:)) = nan;
C_stable_back(ind_stable(:)) = nan;

%% Plot results
[fig,ax] = setupFigure();
hold(ax, 'on'); grid(ax, 'on');

plot3(ax, T_saddle_boundary(:,1), saddle_boundary(:,1), ...
      saddle_boundary(:,2), 'Linewidth', 2*line_width, 'Color', light_blue);
plot3(ax, T_stable_boundary(:,1), stable_boundary(:,1), ...
      stable_boundary(:,2), 'Linewidth', 2*line_width, 'Color', dark_blue);

surf(T_saddle_back, A_saddle_back, C_saddle_back, ...
     'Facecolor', light_blue, ...
     'Edgecolor', 'none', 'FaceAlpha', 0.5,'AmbientStrength',0.5);
surf(T_stable_back, A_stable_back, C_stable_back, ...
     'Facecolor', dark_blue, ...
     'Edgecolor', 'none', 'FaceAlpha', 0.5,'AmbientStrength',0.5);
% plot edges
 plot3(T_saddle_back(1,:), A_saddle_back(1,:), C_saddle_back(1,:), ...
     'color', light_blue, 'Linewidth', 2*line_width);
plot3(T_stable_back(1,:), A_stable_back(1,:), C_stable_back(1,:), ...
     'color', dark_blue, 'Linewidth', 2*line_width);
plot3(T_saddle_back(end,:), A_saddle_back(end,:), C_saddle_back(end,:), ...
     'color', light_blue, 'Linewidth', 2*line_width);
plot3(T_stable_back(end,:), A_stable_back(end,:), C_stable_back(end,:), ...
     'color', dark_blue, 'Linewidth', 2*line_width);
 
% surf(ax, T_saddle_forward, A_saddle_forward, C_saddle_forward, ...
%   'Facecolor', 'b', 'Edgecolor', 'none', 'FaceAlpha', 0.5);
% surf(ax, T_stable_forward, A_stable_forward, C_stable_forward, ...
%   'Facecolor', 'r', 'Edgecolor', 'none', 'FaceAlpha', 0.5);
 
lighting gouraud;
xlit=light('Position',[1 0 0]);
negxlit=light('Position',[-1 0 0]);

set(ax, 'View', [-60.0000 18.7790]);
grid(ax, 'on');
xlabel(ax, 'Time (h)', 'Position', [11.7206,2.2393,-2.5792]);
ylabel(ax, '$A$', 'Position', [0.8880,0.6100,-2.520]);
zlabel(ax, '$C$', 'Rotation', 0, 'Position', [-0.1188,-2.5257,0.7227]);
%title(ax, '$I=50$');

set(ax, 'XLim', [0,24], 'YLim', [-2.0,2.0], 'ZLim', [-2.0,2.0]);
set(ax, 'YDir', 'reverse');
set(ax, 'XTick', 0:12:24, 'YTick', -2:2:2, 'ZTick', -2:2:2);
set(ax, 'Fontsize', 20);

%% Plot orbits
orbit_t = linspace(0, 24, length(orbit)+1);
orbit_t = orbit_t(1:end-1);
plot3(ax, orbit_t, orbit(:,1), orbit(:,2), 'Linewidth', line_width, ...
      'Color', 'black', 'Linestyle', STYs{1});
plotTrajArrow(ax, 800, orbit_t, orbit, 'black', line_width)
plot3(ax, orbit_t, saddle_orbit(:,1), saddle_orbit(:,2), ...
      'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{2});
plot3(ax, orbit_t, unstable_orbit(:,1), unstable_orbit(:,2), ...
      'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{3});

% Plot original manifolds
ind = stable_man(:,1).^2 + stable_man(:,2).^2 <= plot_radius^2;
plot3(ax, 24*ones(sum(ind),1), stable_man(ind,1), stable_man(ind,2), ...
  'Linewidth', line_width, 'Color', dark_blue);
ind = saddle_man(:,1).^2 + saddle_man(:,2).^2 <= plot_radius^2;
plot3(ax, 24*ones(sum(ind),1), saddle_man(ind,1), saddle_man(ind,2), ...
  'Linewidth', line_width, 'Color', light_blue);



%[fig,ax] = postProcessFig(fig, ax, 3, 2, 3, 'manifolds_3D_I_50','B',ax_pos);

function plotTrajArrow(ax, ind, T, orbit, color, line_width)
  t  = T(ind);
  dt = T(ind+1) - t;
  pt = orbit(ind,:);
  dv = orbit(ind+1,:) - pt;
  quiver3( ax, t, pt(1), pt(2), dt, dv(1), dv(2), 150, ...
    'Linewidth', line_width, 'Color', color, 'MaxHeadSize', 1000.0);
end
