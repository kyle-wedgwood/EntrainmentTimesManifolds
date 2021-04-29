function J = KronauerJacobianFS(u, p)
  alpha = p.alpha_0*sqrt(p.I/p.I_0);
  n_inf = alpha/(alpha+p.beta);
  
  B = p.G*alpha*(1-n_inf)*(1-0.4*u(1))*(1-0.4*u(2));
  B_A = -0.4*B/(1-0.4*u(1));
  B_C = -0.4*B/(1-0.4*u(2));
  
  J = (pi/12)*[ p.mu*(1-4.0*u(1)^2) - p.k*u(2)*B_A, ...
                -((24/(0.99669*p.taux))^2+p.k*B) - p.k*B_C*u(2);
                1.0 + B_A , B_C ];
end