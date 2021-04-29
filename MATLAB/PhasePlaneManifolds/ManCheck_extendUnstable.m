%% Plot (A,C) solutions

clear; clc;

close all; 
% Prepare figure
[fig,ax] = setupFigure();

xlabel(ax, '$A$');
ylabel(ax, '$C$', 'Rotation', 0);
set(ax, 'XLim', [-2,2], 'YLim', [-2,2]);
set(ax, 'XTick', -2:2, 'YTick', -2:2);



p = setDefaultParameters();
p = setDayLengthParameters(p, 17);

p.I = 50;
p.tShift = 0.0;
p.taux = 24.2;


orbit_folder = '~/Dropbox/Kyle-Jen-Casey/code/Auto_conts/NewNewB/Manifolds/';

OrbitFileNames = {'Orbit_I_50_UZ3.dat';'Orbit_I_50_UZ4.dat';'Orbit_I_50_UZ5.dat'};

man_name = 'Manifold_I_50_N_17-00_UZ2_'
man_ind = [2,4];
% mans to do: {'2_2','2_4','3_2','3_6'}

%man_name = 'Manifold_I_50_UZ5_'
%man_ind = [1,3,5,7];

% load and plot po's
for mm= 1:length(OrbitFileNames)
    stable_orbit_filename = [orbit_folder, OrbitFileNames{mm}];
    stable_orbit = load(stable_orbit_filename);
    
    light_ind = 1:(length(stable_orbit)-1)/2;
    dark_ind  = light_ind(end)+1:length(stable_orbit);
    
    plot(ax, stable_orbit(light_ind,1), stable_orbit(light_ind,2), ...
        'Linewidth', line_width, 'Linestyle', STYs{mm}, ...
        'Color', grey);
    plot(ax, stable_orbit(dark_ind,1), stable_orbit(dark_ind,2), ...
        'Linewidth', line_width,  'Linestyle', STYs{mm}, ...
        'Color', 'black');
end

%% load and plot manifolds
for kk= 1:length(man_ind)
man = load([orbit_folder man_name num2str(man_ind(kk)) '.dat']);
plot(ax, man(:,1),man(:,2),'Linewidth', line_width+1 ); 

%
% pick a point
man_point = 1;%length(man);
plot(man(man_point,1),man(man_point,2),'o')

% continue the first point on the manifold forward in time
y0 = man(man_point,:)';
odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6, 'Refine', 10);
no_cycles = 100%1e6;

[t,y] = ode45(@Kronauer2SingleFS, [0:24:no_cycles*24], y0, odeopts, p);

plot(ax, y(:,1), y(:,2),  '.', 'Color', 'b', 'Markersize', 15);
y0 = y(end,:)';
%plot(ax, y0(1), y0(2),  '.', 'Color', 'b', 'Markersize',15);


% save extended manifold
% 
new_name = [orbit_folder man_name num2str(man_ind(kk)) '_c.dat'];
fileID = fopen(new_name,'w');
B = [y(:,1), y(:,2)]';
fprintf(fileID,'%12.8f %12.8f\n',B);
            
end          


