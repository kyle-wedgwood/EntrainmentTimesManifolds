function dPhiDt = variationalEquation(t, y, orbit, jacHandle, pars, idx)

  phi = makeMatrix(y);
  pt  = deval(orbit, t, idx);
  J   = feval(jacHandle, pt, pars);

  temp = J*phi;

  dPhiDt = temp(:);

end