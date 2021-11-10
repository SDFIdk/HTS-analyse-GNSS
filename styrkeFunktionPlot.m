function fig = styrkeFunktionPlot(delta,z,sigma,alpha,titlestring)
  %STYRKEFUNKTIONPLOT.m
  % makes a plot of the statistical power and puts 'titlestring'
  % as the title.
  %inputs:
  %     delta: List of x-values from which a y-value is calculates ([-5,5,500]
  %     z    : 
  %     sigma: 
  %     alpha: significance level
  % Formula:
  % 1 - Phi( (z_crit*sigma - delta)/sigma) + Phi( (-z_crit*sigma - delta)/sigma
  %
  %outputs : returns a figure handle (fig) that can be displayed or stored.
  
  fig = figure(1,"visible","off");
  set(gcf, 'Position', get(0, 'Screensize'));
    
  plot(delta,styrkefunktion(delta,z,sigma));
  title(titlestring);
  ylabel('Probability of rejection');
  xlabel('delta [mm/year]');
  
end