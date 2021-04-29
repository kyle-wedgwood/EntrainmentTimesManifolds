%% Plot (A,C) solutions

close all; clear; clc;

p = setDefaultParameters();
p = setDayLengthParameters(p, 12);

p.I = 100;
p.tShift = 0.0;
p.taux = 24.2;

orbit_folder = '~/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds/';


OrbitFileNames = {'Orbit_I_UZ_100_UZ2.dat';'Orbit_I_UZ_100_UZ6.dat';'Orbit_I_UZ_100_UZ8.dat'};

man_name = 'Manifold_I_UZ_100_UZ2_'
man_ind = [3,1];

% Prepare figure
[fig,ax] = setupFigure();

xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);

% load and plot po's
for mm= 1:length(OrbitFileNames)
    stable_orbit_filename = [orbit_folder, OrbitFileNames{mm}];
    stable_orbit = load(stable_orbit_filename);
    
    light_ind = 1:(length(stable_orbit)-1)/2;
    dark_ind  = light_ind(end)+1:length(stable_orbit);
    
    plot(ax, stable_orbit(light_ind,1), stable_orbit(light_ind,2), ...
        'Linewidth', line_width, 'Linestyle', STYs{mm}, ...
        'Color', grey, 'HandleVisibility', 'off');
    plot(ax, stable_orbit(dark_ind,1), stable_orbit(dark_ind,2), ...
        'Linewidth', line_width,  'Linestyle', STYs{mm}, ...
        'Color', 'black', 'HandleVisibility', 'off');
end

% load and plot manifolds
temp_man = load([orbit_folder man_name num2str(man_ind(1)) '.dat']);
for j = 1:length(temp_man)-1
    if norm([temp_man(j+1,1)-temp_man(j,1); temp_man(j+1,2)-temp_man(j,2)]) > 0.2
        break;
    end
end
man = flipud(temp_man(1:j,:));

temp_man = load([orbit_folder man_name num2str(man_ind(2)) '.dat']);
for j = 1:length(temp_man)-1
    if norm([temp_man(j+1,1)-temp_man(j,1); temp_man(j+1,2)-temp_man(j,2)]) > 0.2
        break;
    end
end
man = [man; temp_man(1:j,:)];


plot(man(:,1),man(:,2)); hold on

%%
% pick a point
man_point =  180;%length(man);
plot(man(man_point,1),man(man_point,2),'o')

% continue the first point on the manifold forward in time
y0 = man(man_point,:)';
odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6, 'Refine', 10);
no_cycles = 20%1e6;

[t,y] = ode45(@Kronauer2SingleFS, [0:24:no_cycles*24], y0, odeopts, p);

plot(ax, y(:,1), y(:,2),  '.', 'Color', 'b', 'Markersize', 15);
y0 = y(end,:)';
%plot(ax, y0(1), y0(2),  '.', 'Color', 'b', 'Markersize',15);

%%
xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);