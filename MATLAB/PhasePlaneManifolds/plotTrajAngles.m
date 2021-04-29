x% Plot angle of trajectories of Poincare map in (A,C) plane
% 28.10.2020
% Kyle Wedgwood

close all; clear; clc;

p = setDefaultParameters();
p = setDayLengthParameters(p, 12);

% Adjust parameters
p.tShift = 0;
p.I = 1000;
p.taux = 25.2;

npts = 512;

a = linspace(-2, 2, npts);
c = linspace(-2, 2, npts);
[A,C] = meshgrid(a, c);

odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
FA = zeros(npts, npts);
FC = zeros(npts, npts);
phi = zeros(npts, npts);
nu = zeros(npts, npts);

for i = 1:npts
  a0 = A(:,i);
  c0 = C(:,i);
  [t,y] = ode45(@Kronauer2SingleFS, [0 24], [a0;c0], odeopts, p);
  a1 = y(end,1:npts)';
  c1 = y(end,npts+1:2*npts)';
  
  FA(:,i) = a1-a0;
  FC(:,i) = c1-c0;
  phi(:,i) = atan2(c1-c0, a1-a0);
  nu(:,i) = norm([a1-a0; c1-c0]);
  
  fprintf('Done %d of %d.\n', i, npts);
    
end

% Adjust phi to [0,2pi]
% ind = phi < 0.0;
% phi(ind) = phi(ind)+2*pi;

%% Plot figure
fig = figure;
ax = axes(fig);
hold(ax, 'on');
imagesc(ax, a, c, phi);
phasemap(18); phasebar;
xlabel(ax, 'A');
ylabel(ax, 'C', 'Rotation', 0);
set(ax, 'Xlim', [A(1),A(end)], 'Ylim', [C(1),C(end)]);
set(ax, 'Ydir', 'normal');
set(ax, 'Fontsize', 20);

% Add orbit
folder = '../EntrainmentTimes/orbits/';
% orbit = load([folder 'stable_orbit_FS_I_1000_N_12_taux_24.2.dat']);
% plot(ax, orbit(:,1), orbit(:,2), 'Linewidth', 4, 'Color', 'red');
