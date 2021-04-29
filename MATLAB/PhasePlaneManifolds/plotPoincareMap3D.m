% Plot phase space structures for B_max = 0.1 in 3D
% Kyle Wedgwood
% 31.10.2019

function [ax,fig,orbits] = plotPoincareMap3D(p, orbits, add_title)

  close all; clc;

  if nargin < 1
    p = setDefaultParameters();
  end
  
  colors = { 'red', 'blue', 'black'}; % superceded by linestyles
  colors = { 'black', 'black', 'black'};

  [fig,ax] = setupFigure();

  stability_type = 'stable';
  if nargin > 1
    [stable_sol,start_point,lambda,V] = findPeriodicOrbitDoubleFS( p, stability_type, orbits{1});
  else
    [stable_sol,start_point,lambda,V] = findPeriodicOrbitDoubleFS( p, stability_type);
  end
  plot3( ax, stable_sol.x, stable_sol.y(1,:), stable_sol.y(2,:), ...
    'Linewidth', line_width, 'Color', grey, 'Linestyle', STYs{1});
  pt = 13; % plotting first part of line black to not replace marker
  plot3( ax, stable_sol.x(1:pt), stable_sol.y(1,1:pt), stable_sol.y(2,1:pt), ...
    'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{1});
  plot3( ax, stable_sol.x+12, stable_sol.y(3,:), stable_sol.y(4,:), ...
    'Linewidth', line_width, 'Color', 'black', 'Linestyle', STYs{1});
  pt = 20;
  plot3( ax, zeros( size( stable_sol.x)), stable_sol.y(1,:), stable_sol.y(2,:), ...
    'Linewidth', 1.0, 'Color', grey, 'Linestyle', STYs{1});
  plot3( ax, zeros( 1, pt), stable_sol.y(1,1:pt), stable_sol.y(2,1:pt), ...
    'Linewidth', 1.0, 'Color', 'black', 'Linestyle', STYs{1});
  plot3( ax, zeros( size( stable_sol.x)), stable_sol.y(3,:), stable_sol.y(4,:), ...
    'Linewidth', 1.0, 'Color', 'black', 'Linestyle', STYs{1});
%   plotPointAndEigenvectors( ax, start_point, lambda, V, colors, 2*marker_size);
  
  stability_type = 'saddle';
  
  % Saddle orbit
  if nargin > 1
    [saddle_sol,start_point,lambda,V] = findPeriodicOrbitDoubleFS( p, stability_type, orbits{2});
  else
    [saddle_sol,start_point,lambda,V] = findPeriodicOrbitDoubleFS( p, stability_type);
  end
  plot3( ax, saddle_sol.x, saddle_sol.y(1,:), saddle_sol.y(2,:), ...
    'Linestyle', STYs{2}, 'Linewidth', line_width, 'Color', grey);
  plot3( ax, saddle_sol.x+12, saddle_sol.y(3,:), saddle_sol.y(4,:), ...
    'Linestyle', STYs{2}, 'Linewidth', line_width, 'Color', 'black');
  plot3( ax, zeros( size( saddle_sol.x)), saddle_sol.y(1,:), saddle_sol.y(2,:), ...
    'Linestyle', STYs{2}, 'Linewidth', 1.0, 'Color', grey);
  plot3( ax, zeros( size( saddle_sol.x)), saddle_sol.y(3,:), saddle_sol.y(4,:), ...
    'Linestyle', STYs{2}, 'Linewidth', 1.0, 'Color', 'black');
%   plotPointAndEigenvectors( ax, start_point, lambda, V, colors, 2*marker_size);
  
  % Unstable orbit
  if nargin > 1
    [unstable_sol,start_point,lambda,V] = findPeriodicOrbitDoubleFS( p, stability_type, orbits{3});
  else
    [unstable_sol,start_point,lambda,V] = findPeriodicOrbitDoubleFS( p, stability_type);
  end
  plot3( ax, unstable_sol.x, unstable_sol.y(1,:), unstable_sol.y(2,:), ...
    'Linestyle', STYs{3}, 'Linewidth', line_width, 'Color', grey);
  plot3( ax, unstable_sol.x+12, unstable_sol.y(3,:), unstable_sol.y(4,:), ...
    'Linestyle', STYs{3}, 'Linewidth', line_width, 'Color', 'black');
  plot3( ax, zeros( size( unstable_sol.x)), unstable_sol.y(1,:), unstable_sol.y(2,:), ...
    'Linestyle', STYs{3}, 'Linewidth', 1.0, 'Color', grey);
  plot3( ax, zeros( size( unstable_sol.x)), unstable_sol.y(3,:), unstable_sol.y(4,:), ...
    'Linestyle', STYs{3}, 'Linewidth', 1.0, 'Color', 'black');
%   plotPointAndEigenvectors( ax, start_point, lambda, V, colors, 2*marker_size);
  
  % Add patches
  h = gobjects(3,1);
  limit = 1.5;
  X = 12*ones( 1, 4);
  Y = limit*[-1,-1,1,1];
  Z = limit*[-1,1,1,-1];
  h(1) = fill3( ax, X, Y, Z, 0.7*ones( 1, 3));
  
  X = [12 12 24 24];
  Y = limit*ones( 1, 4);
  Z = limit*[-1,1,1,-1];
  h(2) = fill3( ax, X, Y, Z, 0.7*ones( 1, 3));
  
  X = [12 12 24 24];
  Y = limit*[-1,1,1,-1];
  Z = limit*ones( 1, 4);
  h(3) = fill3( ax, X, Y, Z, 0.7*ones( 1, 3));
  
  set( h, 'FaceAlpha', 0.4, 'EdgeColor', 'black');
  
  % Identify end points with start points
  pt = plot3( ax, 0, stable_sol.y(1,1), stable_sol.y(2,1), ...
    'Marker', 'o', 'Markersize', marker_size, 'Color', 'black', ...
    'MarkerFaceColor', 'black');pause(1);
  for sz = 1:-0.3:0.2
    pt = plot3( ax, 0, saddle_sol.y(1,1), saddle_sol.y(2,1), ...
      'Marker', 'o', 'Markersize', sz*marker_size, 'Color', 'black');
  end
  pt = plot3( ax, 0, unstable_sol.y(1,1), unstable_sol.y(2,1), ...
    'Marker', 'o', 'Markersize', marker_size, 'Color', 'black');

  plot3( ax, 24, stable_sol.y(3,end), stable_sol.y(4,end), ...
    'Marker', 'o', 'Markersize', marker_size, 'Color', 'black', ...
    'MarkerFaceColor', grey);
  for sz = 1:-0.3:0.2
    pt = plot3( ax, 24, saddle_sol.y(3,end), saddle_sol.y(4,end), ...
      'Marker', 'o', 'Markersize', sz*marker_size, 'Color', grey);
  end
  plot3( ax, 24, unstable_sol.y(3,end), unstable_sol.y(4,end), ...
    'Marker', 'o', 'Markersize', marker_size, 'Color', grey);
  
  % Plot arrows on trajectories
%   plotTrajArrow( ax, 100, stable_sol, 'black', line_width);
%   plotTrajArrow( ax, 200, saddle_sol, 'blue');
%   plotTrajArrow( ax, 150, unstable_sol, 'red');
  
  axis image;
  set( ax, 'XLim', [0,24], 'YLim', limit*[-1 1], 'ZLim', limit*[-1 1]);
  set( ax, 'XTick', 0:4:24, 'YTick', -1:1, 'ZTick', -1:1);
  xlabel( ax, 'Time (h)');
  ylabel( ax, '$A$');
  zlab = zlabel( ax, '$C$', 'Rotation', 0); % 'Position',);
  if add_title
    title( ax, sprintf( '$I: %4.f$', p.I));
  end
  
  set( ax, 'YDir', 'reverse');
  grid( ax, 'on');
  
  set( ax, 'Fontsize', 20.0);
  set( ax, 'View',  [-45.7000   26.2307]); pause(1);
  pos = get( zlab, 'Position');
  set( zlab, 'Position', pos + [0,0,0.4]);

  
  orbits = { stable_sol, saddle_sol, unstable_sol};

end

function plotPointAndEigenvectors( ax, pt, lambda, V, colors, marker_size)

  eps = 0.2;

  for i = 1:length( lambda)
    if abs( lambda(i)) < 1
      plotAttractingArrows( ax, pt, V(:,i), eps);
    else
      plotRepellingArrows( ax, pt, V(:,i), eps);
    end
  end

  plot3( ax, 0, pt(1), pt(2), 'Marker', '.', ...
    'Markersize', marker_size, 'Color', colors{ sum( abs( lambda) < 1)+1});
  
end

function plotAttractingArrows( ax, pt, v, eps)
  for val = -eps:2*eps:eps
    alt_pt = pt+1.1*val*v;
    quiver3( ax, 0, alt_pt(1), alt_pt(2), 0, -val*v(1), -val*v(2), 'AutoScale', 'on', ...
      'Linewidth', 2.0, 'Color', 'black', 'MaxHeadSize', 400.0);
  end
end

function plotRepellingArrows( ax, pt, v, eps)
  for val = -eps:2*eps:eps
    quiver3( ax, 0, pt(1), pt(2), 0, val*v(1), val*v(2), 'AutoScale', 'on', ...
      'Linewidth', 2.0, 'Color', 'red', 'MaxHeadSize', 2.0);
  end
end

function plotTrajArrow( ax, ind, sol, color, line_width)
  t  = sol.x(ind);
  dt = 0.01;
  pt = sol.y(:,ind);
  dv = deval(sol, t+dt) - pt;
  quiver3( ax, t, pt(1), pt(2), dt, dv(1), dv(2), 30, ...
    'Linewidth', line_width, 'Color', color, 'MaxHeadSize', 10.0);
end