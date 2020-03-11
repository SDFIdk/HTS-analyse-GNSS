function fig = styrkeFunktionPlot(delta,z,sigma,alpha,titlestring)
  % Makes a plot of the statistical power and puts 'titlestring' as the title.
  %inputs:
  %     delta: linspace(-5,5,500);
  %     z    : stats.z_crit (Z-test, calculated in linreg on line 172).
  %     sigma: stats.std_known ( Standard error of slope (known), calculated in linreg on line 166)
  %     alpha: significance level
  %     titlestring: ex. BALN, Statistical Power for Z-test at significance level 0.05, N = 1, N_{binned} = 1 (defined on line 460, TimeSeriesAnalysis).
  % Formula:
  % 1 - Phi( (z_crit*sigma - delta)/sigma) + Phi( (-z_crit*sigma - delta)/sigma
  % where Phi == normcdf.
  %
  %outputs : returns a figure handle (fig) that can be displayed or stored.
  
  fig = figure(1,"visible","on");
  set(gcf, 'Position', get(0, 'Screensize'));
    
  plot(delta,styrkefunktion(delta,z,sigma));
  title(titlestring);
  ylabel('Probability of rejection');
  xlabel('delta [mm/year]');
  
end