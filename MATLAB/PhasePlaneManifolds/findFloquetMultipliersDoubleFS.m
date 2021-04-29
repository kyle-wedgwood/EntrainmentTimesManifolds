function [nu,v,var_sol] = findFloquetMultipliersDoubleFS(orbit, jacHandle, p)

  debug = false;

  I = eye(2);
  options = odeset('RelTol', 1e-6, 'Abstol', 1e-8);

  fun = @(t,y) variationalEquation(t, y, orbit, ...
    @(u,p) jacHandle(u, p), p, [1,2]);
  [~,sol_1] = ode45(fun, orbit.x, I(:), options);
  phi_1 = makeMatrix(sol_1(end,:)');

  fun = @(t,y) variationalEquation(t, y, orbit, ...
    @(u,p) jacHandle(u, p), p, [3,4]);
  [~,sol_2] = ode45(fun, orbit.x, I(:), options);
  phi_2 = makeMatrix(sol_2(end,:)');
  
  var_sol = [sol_1, sol_2];

  phi = phi_2*phi_1;
  [v,nu] = eig(phi);
  
  nu = diag(nu);
  
  if debug
    [v_1,nu_1] = eig(phi_1);
    nu_1 = diag(nu_1);
    
    [v_2,nu_2] = eig( phi_2);
    nu_2 = diag( nu_2);
  end
  
end