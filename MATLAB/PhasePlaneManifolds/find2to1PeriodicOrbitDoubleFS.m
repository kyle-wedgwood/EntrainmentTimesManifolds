% Attempting to find 2:1 periodic orbit
% Returns single point along orbit
% Kyle Wedgwood
% 4.12.2020

function x0 = find2to1PeriodicOrbitDoubleFS(p)

  close all;
  
  plot_flag = true;
  
  % Pick initial condition
  x0 = [1,0.5];
  fun = @(x) integrateKronauerFS(x, p, 1) - x;
  
  x0 = fsolve(fun, x0, optimset('Display', 'iter'));
  
  if nargin > 2
    initialGuess = @(t) deval(orbit,t);
  else
    % fill in this with numerically integrated guess
    odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
    u0 = [1;0];
    p = setDayLengthParameters(p, p.day_length);
    tspan = [0 24];
    while true
      sol = ode45(@Kronauer2SingleFS, tspan, u0, odeopts, p);
      if norm(sol.y(:,1)-sol.y(:,end)) < 0.001
        break;
      else
        u0 = sol.y(:,end);
      end
    end
    initialGuess = @(t) [deval(sol,t); deval(sol,t+p.day_length)];
  end
  
  npts = 5;
  time = linspace(0, 12, npts);
  solinit = bvpinit(time, initialGuess);
  bvpopts = bvpset('AbsTol', 1e-12, 'RelTol', 1e-8, 'Vectorized', true);
  
  sol = bvp4c(@(t,u) KronauerDoubleFS(t, u, p), ...
              @(ya,yb) boundaryConditions(ya,yb), solinit, bvpopts);
  fprintf('Periodic orbit found.\n');
  start_point = sol.y(1:2,1);
  
  % Find Floquet multipliers
  if nargout > 2
    [lambda,V,var_sol] = findFloquetMultipliersDoubleFS(sol, @KronauerJacobianFS, p);
    for i = 1:length( lambda)
      fprintf( 'Floquet multiplier %d: %.4f\n',  i, lambda(i));
    end
  end

  if plot_flag
    %% Plot solution + data
    phase_fig = figure;
    phase_ax = axes( phase_fig);
    hold( phase_ax);
        
    plot(phase_ax, solinit.y(1,:), solinit.y(2,:), 'bo', solinit.y(3,:), solinit.y(4,:), 'bo');
    plot(phase_ax, sol.y(1,:), sol.y(2,:), 'Linewidth', 4.0, 'Color', 'red');
    plot(phase_ax, sol.y(3,:), sol.y(4,:), 'Linewidth', 4.0, 'Color', 'black');
    plot(phase_ax, sol.y(1,1), sol.y(2,1), 'Marker', 'o', 'Markersize', 20, 'Color', 'red');
    plot(phase_ax, sol.y(3,1), sol.y(4,1), 'Marker', 'o', 'Markersize', 20, 'Color', 'black');
    
    xlabel(phase_ax,'A');
    ylabel(phase_ax, 'C');
    title(phase_ax, 'Two point BVP solution');
    set(phase_ax, 'Fontsize', 20);
    
    % Now add Floquet multipliers
    if nargout > 2
      for i = 1:length(lambda)
        if abs(lambda(i)) < 1
          plotAttractingArrows(phase_ax, sol.y(1:2,1), V(:,i), 0.2);
        else
          plotRepellingArrows(phase_ax, sol.y(1:2,1), V(:,i), 0.2);
        end
      end
    end
  end
  
end

function x = integrateKronauerFS(x, p, no_cycles)
  % Integrate Kronauer model 24 hours n times
  odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
  for n = 1:no_cycles
    [~,y] = ode45(@Kronauer2SingleFS, [0 24], x, odeopts, p);
    x = y(end,:)';
  end
end

function plotAttractingArrows(ax, pt, v, eps)
  for val = -eps:2*eps:eps
    alt_pt = pt+1.1*val*v;
    quiver(ax, alt_pt(1), alt_pt(2), -val*v(1), -val*v(2), 'AutoScale', 'on', ...
      'Linewidth', 2.0, 'Color', 'black', 'MaxHeadSize', 2.0);
  end
end

function plotRepellingArrows(ax, pt, v, eps)
  for val = -eps:2*eps:eps
    quiver(ax, pt(1), pt(2), val*v(1), val*v(2), 'AutoScale', 'on', ...
      'Linewidth', 2.0, 'Color', 'red', 'MaxHeadSize', 2.0);
  end
end
 
