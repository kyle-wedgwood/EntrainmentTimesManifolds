function p = setDefaultParameters()

  p.B_min  = 0.0;
  p.B_max  = 0.1;
  p.mu     = 0.23;
  p.taux   = 24.2;
  p.k      = 0.55;
  p.tShift = 18;
  p.day_length = 12;
  
  % New parameters for fast-slow system
  p.G = 33.75;
  p.alpha_0 = 0.05;
  p.I_0 = 9500.0;
  p.I = 50.0; % lux
  p.beta = 0.0075;
  
end