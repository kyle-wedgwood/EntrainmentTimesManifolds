function dudt = Kronauer2SingleFS_lightbox(t, u, p)
  N = size(u,1)/2;
  A = u(1:N);
  C = u(N+1:2*N);
  
  ft = (p.length_scaling+sin(2.0*(pi/24.0)*(t-p.tShift-p.length_shift)) > 0);
  
  switch p.type
    case 'delta'
      alpha = p.alpha_0*sqrt((p.I+p.dI)/p.I_0);
      B_light = p.G*alpha*(1-alpha/(alpha+p.beta))*(1-0.4*C).*(1-0.4*A);

      alpha = p.alpha_0*sqrt(p.dI*(p.dI > 0)/p.I_0);
      B_dark  = p.G*alpha*(1-alpha/(alpha+p.beta))*(1-0.4*C).*(1-0.4*A);

      B = B_dark + (B_light-B_dark).*ft;

    case 'fixed_light'
      alpha = p.alpha_0*sqrt((p.I + (p.I_fixed-p.I)*(p.I_fixed > p.I))/p.I_0);
      B_light = p.G*alpha*(1-alpha/(alpha+p.beta))*(1-0.4*C).*(1-0.4*A);

      alpha = p.alpha_0*sqrt(p.I_fixed/p.I_0);
      B_dark  = p.G*alpha*(1-alpha/(alpha+p.beta))*(1-0.4*C).*(1-0.4*A);

      B = B_dark + (B_light-B_dark).*ft;

  end
  
  dudt = pi/12.0*[ ( p.mu*(A-(4/3)*A.^3)-C.*(((24/(0.99669*p.taux))^2+p.k*B)));
                   ( A+B)];
end