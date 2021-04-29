function [fig,ax] = setupFigure()

  set(groot,'defaultAxesTickLabelInterpreter', 'latex'); % set the axis font to latex 
  set(groot, 'defaultLegendInterpreter', 'latex'); % and for the legend

  STYs = {'-','--',':'};  % Lines styles, stable, saddle and unstable
  line_width = 2; % for plots
  marker_size = 8;

  % Line colors
  grey = 0.75*ones(1, 3);
  blue = [0.0 0.45 0.74];
  red  = [0.64 0.08 0.18];

  fcol  = [0.5 0 0.5]; % for FOLD
  pdcol  = [0.5 0.5 0]; % for Period doubling
  nscol = [0 0.5 0.5]; % for Neimark-Sacker

  ax_pos = [0.15 0.15 0.85 0.70];

  vars = whos;
  for i = 1:length(vars)
    assignin('caller', vars(i).name, eval(vars(i).name));
  end

  fig = figure;
  ax = axes(fig);
  hold(ax, 'on');

  % Set interprete as latex
  labs = {'Xlabel', 'Ylabel', 'Zlabel', 'Title'};
  for i = 1:length(labs)
    lab = get(ax, labs{i});
    set(lab, 'Interpreter', 'latex');
  end  

end