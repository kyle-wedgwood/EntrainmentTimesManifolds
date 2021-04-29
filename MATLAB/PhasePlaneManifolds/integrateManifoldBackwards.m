function [T,A,C] = integrateManifoldBackwards(manifold, p, no_pts)

  t = linspace(0, 1, size(manifold, 1))';
  manifold = interp1(t, manifold, linspace(0, 1, no_pts)');
  
  odeopts = odeset('AbsTol', 1e-10, 'RelTol', 1e-8, ...
                   'Refine', 10, 'Events', @domainEscape);
  dt = 0.1;
  tspan = -(0:dt:24)';
  A = nan*ones(length(tspan), size(manifold, 1));
  C = nan*ones(length(tspan), size(manifold, 1));
  
  for i = 1:size(manifold, 1)
    [t,y,te,ye,ie] = ode45(@Kronauer2SingleFS, tspan,...
      [manifold(i,1); manifold(i,2)], odeopts, p);
    no_pts = length(t);
    A(1:no_pts,i) = y(:,1);
    C(1:no_pts,i) = y(:,2);
  end

  no_pts = size(y,2)/2;
  T = 24+repmat(tspan, [1,size(manifold, 1)]);

end

function [val,isterm,dir] = domainEscape(~, y, p)

  N = size(y, 1)/2;
  
  val    = y(1:N).*y(1:N)+y(N+1:2*N).*y(N+1:2*N)-p.domain_radius^2;
  isterm = ones(N, 1);
  dir    = ones(N, 1);

end