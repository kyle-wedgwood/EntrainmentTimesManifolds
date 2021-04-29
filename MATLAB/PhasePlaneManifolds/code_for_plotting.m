%% some bits of code to format the figures

% In set up
set(groot,'defaultAxesTickLabelInterpreter','latex'); % set the axis font to latex 
set(groot, 'defaultLegendInterpreter','latex'); % and for the legend

STYs = {'-','--',':'};  % Lines styles, stable, saddle and unstable
line_width = 2; % for plots
marker_size = 8;

fcol  = [0.5 0 0.5]; % for FOLD
pdcol  = [0.5 0.5 0]; % for Period doubling
nscol = [0 0.5 0.5]; % for Neimark-Sacker

% during plotting
figure; hold on
xlabel('$x$','interpreter','latex') % and for the axis labels
ylabel('$y$','interpreter','latex')
% plot some stuff...

% At the end to set everything the right size and save
page_width = 19; % or text width, this will be journal dependent.
W = page_width/3; % or 2 ifyou want half page or 1 for full page..
H = page_width/3;

fnam = 'My_awesome_fig';
set(gca,'FontSize',10,'linewidth',1) % font size 10 and axis line width 1
set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[W,H])
set(gcf,'paperposition',[0,0,W,H])
print(['figures/' fnam],'-depsc') % save as eps, or any other format

