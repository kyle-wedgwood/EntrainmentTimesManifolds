function [A,C] = shiftManifold(manifold, p, no_pts, shift)

  t = linspace(0, 1, size(manifold, 1))';
  manifold = interp1(t, manifold, linspace(0, 1, no_pts)');
  
  odeopts = odeset('AbsTol', 1e-10, 'RelTol', 1e-8, ...
                   'Refine', 10, 'Events', @domainEscape);
  A = nan*ones(size(manifold, 1), 1);
  C = nan*ones(size(manifold, 1), 1);
  
  for i = 1:size(manifold, 1)
    [t,y,te,ye,ie] = ode45(@Kronauer2SingleFS, [0 shift],...
      [manifold(i,1); manifold(i,2)], odeopts, p);
    A(i) = y(end,1);
    C(i) = y(end,2);
  end

end

function [val,isterm,dir] = domainEscape(~, y, p)

  N = size(y, 1)/2;
  
  val    = y(1:N).*y(1:N)+y(N+1:2*N).*y(N+1:2*N)-p.domain_radius^2;
  isterm = ones(N, 1);
  dir    = ones(N, 1);

end