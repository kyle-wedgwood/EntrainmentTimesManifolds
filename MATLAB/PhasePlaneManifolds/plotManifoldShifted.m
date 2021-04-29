% Plot manifolds integrated by tShift updated for NewNewB
% 4.11.2020
% Kyle Wedgwood

close all; clear; clc;

p = setDefaultParameters();
p = setDayLengthParameters(p, 12);

orbit_folder = '~/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds/';

stable_orbit_filename = [orbit_folder, 'Orbit_I_50_UZ3.dat'];
stable_orbit = load(stable_orbit_filename);
saddle_orbit_filename = [orbit_folder, 'Orbit_I_50_UZ4.dat'];
saddle_orbit = load(saddle_orbit_filename);
unstable_orbit_filename = [orbit_folder, 'Orbit_I_50_UZ5.dat'];
unstable_orbit = load(unstable_orbit_filename);
man_name = 'Manifold_I_50_UZ*';
man_ind = [3,7];
files = dir(['manifolds/' man_name]);
fig_filename = 'shifted_manifolds_I_50_tau_24.2_N_12';

p.I = 50;
p.tShift = 0.0;

% Prepare figure
[fig,ax] = setupFigure();

xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);

light_ind = 1:(length(stable_orbit)-1)/2;
dark_ind  = light_ind(end)+1:length(stable_orbit);
plot(ax, stable_orbit(light_ind,1), stable_orbit(light_ind,2), ...
     'Linewidth', line_width, ...
     'Color', grey, 'HandleVisibility', 'off');
plot(ax, stable_orbit(dark_ind,1), stable_orbit(dark_ind,2), ...
     'Linewidth', line_width, ...
     'Color', 'black', 'HandleVisibility', 'off');
plot(ax, saddle_orbit(light_ind,1), saddle_orbit(light_ind,2), ...
     'Linewidth', line_width, 'Linestyle', STYs{2}, ...
     'Color', grey, 'HandleVisibility', 'off');
plot(ax, saddle_orbit(dark_ind,1), saddle_orbit(dark_ind,2), ...
     'Linewidth', line_width, 'Linestyle', STYs{2}, ...
     'Color', 'black', 'HandleVisibility', 'off');
plot(ax, unstable_orbit(light_ind,1), unstable_orbit(light_ind,2), ...
     'Linewidth', line_width, 'Linestyle', STYs{3}, ...
     'Color', grey, 'HandleVisibility', 'off');
plot(ax, unstable_orbit(dark_ind,1), unstable_orbit(dark_ind,2), ...
     'Linewidth', line_width, 'Linestyle', STYs{3}, ...
     'Color', 'black', 'HandleVisibility', 'off');

odeopts = odeset('AbsTol', 1e-12, 'RelTol', 1e-10);

% tShift = 0:1:22; % forwards
tShift = -(0:1:22); % backwards
type = 1;
npts = 1000;

cm = lines(10);
cm(2,:) = [];
cm = [cm,ones(9,1)];
cm = repmat(cm, [1,1,2]);
cm(:,4,2) = 0.2*ones(9,1);

for i = man_ind
  temp_man = load(['manifolds/' files(i).name]);
  for j = 1:length(temp_man)-1
    if norm([temp_man(j+1,1)-temp_man(j,1); temp_man(j+1,2)-temp_man(j,2)]) > 0.2
      break;
    end
  end
  man = flipud(temp_man(1:j,:));
  
  temp_man = load(['manifolds/' files(i+1).name]);
  for j = 1:length(temp_man)-1
    if norm([temp_man(j+1,1)-temp_man(j,1); temp_man(j+1,2)-temp_man(j,2)]) > 0.2
      break;
    end
  end
  man = [man; temp_man(1:j,:)];
  
  t = linspace(0, 1, size(man, 1))';
  man = interp1(t, man, linspace(0, 1, npts)');

  % output variable
  A_shift = nan*ones(length(man), length(tShift));
  C_shift = nan*ones(length(man), length(tShift));

  % Integrate manifolds
  for j = 1:length(man)
    [t,y] = ode45(@Kronauer2SingleFS, tShift, [man(j,1),man(j,2)], odeopts, p);

    % Stick value in man
    for k = 1:length(t)
      A_shift(j,k) = y(k,1);
      C_shift(j,k) = y(k,2);
    end
    
  end

  for j = 1:6:length(tShift)
    h = plot(ax, A_shift(:,j), C_shift(:,j), ...
            'Color', cm((j-1)/6+1,:,type), 'Linestyle', '-', ...
            'DisplayName', sprintf('$t_s$ = %d', j-1), ...
            'Linewidth', 0.5*line_width);
    if type == 2
      set(h, 'HandleVisibility', 'off');
    end
  end

  type = type+1;

end

% legend(ax, 'Location', 'NorthEast');
text(ax, 1.1, -0.85, '$t_s$ = 0', 'Interpreter', 'latex');
text(ax, -1.3, -1.5, '$t_s$ = 6', 'Interpreter', 'latex');
text(ax, -1.95, 0.85, '$t_s$ = 12', 'Interpreter', 'latex');
text(ax, 0.5, 1.5, '$t_s$ = 18', 'Interpreter', 'latex');
axis square
[fig,ax] = postProcessFig(fig, ax, 3, 1, 3, fig_filename, 'A', ax_pos);