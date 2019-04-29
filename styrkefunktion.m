function beta = styrkefunktion(delta,z,sigma)
  % STYRKEFUNKTION.m 
  % 
  %
%1 - Phi( (z_crit*sigma - delta)/sigma) + Phi( (-z_crit*sigma - delta)/sigma
    
beta = 1 - normcdf((z.*sigma - delta)./sigma) + normcdf( (-z.*sigma - delta)./sigma);

end