function [fig,ax] = postProcessFig(fig, ax, ncol, w_scale, nrow, fnam, label, ax_pos)

  page_width = 18; % or text width, this will be journal dependent. Frontiers = 180mm
  W = w_scale*page_width/ncol; % or 2 ifyou want half page or 1 for full page..
  H = page_width/nrow; % needs to be modified for 2 para bif diags

  box(ax, 'on');
  
  set(ax, 'FontSize', 10, 'linewidth', 1) % font size 10 and axis line width 1
  set(fig, 'paperunits', 'centimeters')
  set(fig, 'papersize', [W,H])
  set(fig, 'paperposition', [0,0,W,H])
  
  set(ax.XLabel, 'Rotation', 0);
  set(ax.YLabel, 'Rotation', 0); % needs to be ignored for 1 para bif diags

  % Sort out colorbar
  % Reposition colorbar
  try 
    c = ax.Colorbar;
    pos = get( c, 'Position');
    pos(1) = pos(1)-3*pos(3);
    pos(3) = pos(3)/2;

    if exist('ax_pos', 'var')
      pos(2) = ax_pos(2);
      pos(4) = ax_pos(4);
    end

    set( c, 'Position', pos);
    set( c, 'Linewidth', 0.5);
    set( c, 'TickLabelInterpreter', 'latex');
    set( get( c, 'label'), 'String', '$\mathcal{T}$', 'Interpreter', 'latex', ...
                           'Rotation', 0);
    pos = get( get( c, 'label'), 'Position');
    pos(1) = 2*pos(1);
    pos(2) = 0.55*c.Limits(2);
    set( get( c, 'label'), 'Position', pos);

    c.FontSize = 10;
  catch
    disp('No colorbar detected.');
  end
  
  % Add figure label
  if nargin > 5
    set(ax,'position', ax_pos)  
    annotation('textbox',[.0 .7 0.3 .3],'String', ['(' label ')'], ...
               'FitBoxToText', 'on', 'FontSize', 11, ...
               'Linestyle', 'none', 'interpreter', 'latex');
  end

  figure(fig); pause(1);
  print([fnam '.eps'], '-depsc', '-painters', '-r300'); % save as eps, or any other format
%   print([fnam '.png'], '-dpng', '-r300'); % save as eps, or any other format

end