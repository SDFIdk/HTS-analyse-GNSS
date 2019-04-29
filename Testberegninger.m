% Testberegninger.m

alpha = 0.05;
sigma0_all_list = [4.0456 3 5];
filename_output = 'outputs\\testberegninger.csv';

figures = 1; % 0 no figures, 1 show, 2 save, 3 show and save
testberegninger_folder = 'testberegninger\\';
figure_format = 'pdf';

% Scenarios List of years to consider [start year: interval : end year]
T = {[2005:3:2032];...
     [2005:6:2029];...
     [2005:9:2032];...
     [2005:3:2011 2017:6:2029];...
     [2005:3:2011 2020 2029];...
     [2005:3:2014 2020:6:2032];...
     [2005:3:2014 2023 2032]};

% Template for output    
template = "Scenario,sigma(%0.2f),sigma(%0.2f),sigma(%0.2f),rel(%0.2f),rel(%0.2f),rel(%0.2f) \n";
output_string = sprintf(template,sigma0_all_list,sigma0_all_list);

% Datapoints for 1st axis in statistical power plots.
delta = linspace(-2,2,500);  

% Critical value for Z-test
z_crit = norminv(1-(alpha)/2);

if figure > 0 
  fig = figure(1,"visible","off");
  set(gcf, 'Position', get(0, 'Screensize'));
end
    
     
% For each scenario do:
for i = 1:length(T)
  for j = 1:length(sigma0_all_list)
    t = T{i};
    sigma0_all = sigma0_all_list(j);
    % Vandermonde matrix
    A = horzcat(ones(length(t),1),t');
  
    fprintf('Scenario %i \n',i);
    SIGMA = (A'*(sigma0_all.^(-2).*eye(length(t)))*A)^(-1);
    
    % Estimated Standar Error of Slope, (p.127) Chapter 5: Regression models
    sigma_b_est = sqrt(SIGMA(2,2))
  
    if i == 1
      sigma_b_est_ref(j) = sigma_b_est;
      sigma_b_rel = 1;
    else
      sigma_b_rel = sigma_b_est/sigma_b_est_ref(j);
    end
    temp_sigmas(j) = sigma_b_est;
    temp_rels(j) = sigma_b_rel;
    
  end
  output_string = [output_string sprintf("%i, %f, %f, %f, %f, %f, %f \n",...
                   i, temp_sigmas, temp_rels)];
                   

    
  %fig = styrkeFunktionPlot(delta,z_crit,sigma_b_est,alpha,...
          %styrke_title_string);
          
  plot(delta,styrkefunktion(delta,z_crit,sigma_b_est),['-;Scenario #' num2str(i) ';']);
  hold on;        
end

title(['Test Calculations for scenarios 1-7. ' 'Statistical Power for Z-test at significance level ' ...
       num2str(alpha)]);

ylabel('Probability of rejection');
xlabel('delta [mm/year]');          
%legend(fig,"location","southoutside");
%legend(fig,"boxoff");

%Show strength curve figure
if (figures == 1) || (figures == 3)
  set(fig, 'visible', 'on');
end

hold off;          
%Save strength curve figure
if (figures == 2) || (figures == 3)
  %set(gcf,'PaperUnits','normalized');
  %set(gcf,'PaperPosition', [0 0 1 1]);
  filename = sprintf('%s.%s','Scenarios',figure_format);
  print(fig, [testberegninger_folder filename]);
end 

fileID = fopen(filename_output,'w');
fprintf(fileID,output_string);
fclose(fileID);

fprintf('.CSV File written.\n')