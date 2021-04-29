function dudt = KronauerDoubleFS( ~, u, p)
  A1 = u(1);
  C1 = u(2);
  A2 = u(3);
  C2 = u(4);
  
  alpha = p.alpha_0*sqrt(p.I/p.I_0);
  n_inf = alpha/(alpha+p.beta);
  B_dark = 0.0;
  B_light = p.G*alpha*(1-n_inf)*(1-0.4*C1).*(1-0.4*A1);
  
  day_adjust = p.day_length/12.0;
  night_adjust = (24-p.day_length)/12.0;
         
  dudt = [ day_adjust*(pi/12)*( p.mu*(A1-(4/3)*A1.^3)-C1.*(((24/(0.99669*p.taux))^2+p.k*B_light)));
           day_adjust*(pi/12)*( A1+B_light);
           night_adjust*(pi/12)*( p.mu*(A2-(4/3)*A2.^3)-C2.*(((24/(0.99669*p.taux))^2+p.k*B_dark)));
           night_adjust*(pi/12)*( A2+B_dark)];
end