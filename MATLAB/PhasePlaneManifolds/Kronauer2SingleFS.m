function dudt = Kronauer2SingleFS(t, u, p)
  N = size(u,1)/2;
  A = u(1:N);
  C = u(N+1:2*N);
  
  ft = (p.length_scaling+sin(2.0*(pi/24.0)*(t-p.tShift-p.length_shift)) > 0);
  
  alpha = p.alpha_0*sqrt(p.I/p.I_0);
  n_inf = alpha/(alpha+p.beta);
  B     = p.G*alpha*ft*(1-n_inf)*(1-0.4*C).*(1-0.4*A);
  
  dudt = pi/12.0*[ ( p.mu*(A-(4/3)*A.^3)-C.*(((24/(0.99669*p.taux))^2+p.k*B)));
                   ( A+B)];
end