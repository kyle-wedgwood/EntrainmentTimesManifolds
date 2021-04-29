function [T,A,C] = integrateManifoldForwards(manifold, p, no_pts)

  t = linspace(0, 1, size(manifold, 1))';
  manifold = interp1(t, manifold, linspace(0, 1, no_pts)');

  % Now integrate manifolds
  tspan = [0 24];
  odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6, 'Refine', 10);
  [t,y] = ode45(@Kronauer2SingleFS, tspan, ...
    [manifold(:,1); manifold(:,2)], odeopts, p);

  no_pts = size(y,2)/2;
  T = repmat(t, [1,no_pts]);
  A = y(:,1:no_pts);
  C = y(:,no_pts+1:2*no_pts);

end