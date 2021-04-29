% Wrapper to plot 3D phase portraits for new new B function
% 28.9.2020
% Kyle Wedgwood

% Load parameters
p = setDefaultParameters();
p = setDayLengthParameters(p, 12);

% Load orbits
time = load('orbits/Time.dat');
stable_orbit = load('orbits/Orbit_lux_50_UZ2.dat');
saddle_orbit = load('orbits/Orbit_lux_50_UZ3.dat');
unstable_orbit = load('orbits/Orbit_lux_50_UZ4.dat');

% Reshape orbits
light_ind = find(time <= 12);
dark_ind  = find(time >= 12);
stable_orbit = [stable_orbit(light_ind,1:2), stable_orbit(dark_ind,1:2)];
saddle_orbit = [saddle_orbit(light_ind,1:2), saddle_orbit(dark_ind,1:2)];
unstable_orbit = [unstable_orbit(light_ind,1:2), unstable_orbit(dark_ind,1:2)];

orbits = {@(t) interp1(time(light_ind), stable_orbit, t), ...
          @(t) interp1(time(light_ind), saddle_orbit, t), ...
          @(t) interp1(time(light_ind), unstable_orbit, t)};
        
add_title = false;
[ax,fig,orbits] = plotPoincareMap3D(p, orbits, add_title);
[fig,ax] = postProcessFig(fig, ax, 3, 2, 3, 'battenberg_I_50', 'D');