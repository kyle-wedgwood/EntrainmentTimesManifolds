% Formulate numerical Jacobian

function J = JacNumGen(y, fhandle)

  n = length(y);
  I = eye(n);
  epsi = 1e-4;
  
  J = zeros(n,n);
  f0 = fhandle(y);
  
  for i = 1:n
    J(:,i) = (fhandle(y+epsi*I(:,i)) - f0)/epsi;
  end

end