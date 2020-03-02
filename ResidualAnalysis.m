%RESIDUALANALYSIS
%
%
% 
% oldjo@sdfe.dk, November 2017

% Filename (Set here if not already specified in calling Script)
if ~exist('filename_input')
   filename_input = 'inputs\\5D_tidsserier_20180215.csv';
end
% Note: this must be a comma seperated file with the following columns
% REFNR,GPSNR,XKOOR,YKOOR,ZKOOR,XRMS,YRMS,ZRMS,JNR_BSIDE,SYS,EPOCH,IN_DATE
filename_5d_points = 'inputs\\5d_points_sorted.csv';

% output filenames
filename_output_csv = 'residuals\\residuals.csv';
resid_fig_folder = 'residuals\\';

% Alpha
if ~exist('alpha')
  alpha = 0.05;
else
  fprintf(['alpha already defined as ' num2str(alpha) '\n']);
end

% Binning (pre binning of observations, see 'binning.m' for the inner workings)
do_binning = 1; %1: on, 0: off.
binsize = 14; % Binsize in days

% Generate figures (0: no, 1: display, 2: save and display).
residuals_figures = 0; %Denne del virker ikke uden funktionen "residuals", fra statistics pakken, som ikke er implementeret i Octave endnu. Prøv i Matlab?
%numbins = 25; %number of bins in histograms
residuals_closeFigures = 0;
Nmin = 3;

% Figure format: (png, pdf, jpg... etc.)
figure_format = 'pdf'; 

% Formatstring (change this if the input files get different fields.
T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
[REFNR, GPSNR, XKOOR, YKOOR, ZKOOR, SYS, EPOCH] =...
textread(filename_input, T, "delimiter", ",", "headerlines", 1);
[GPSNR_5D] = textread(filename_5d_points, "%s");

% Find Unique entries
[Y, I, J] = unique (REFNR, "first");

%Initialize strings, arrays, cells etc.
sigma_0_sqrd_top_h = [];
sigma_0_sqrd_top_n = [];
sigma_0_sqrd_top_e = [];
sigma_0_sqrd_bottom_h = [];
sigma_0_sqrd_bottom_n = [];
sigma_0_sqrd_bottom_e = [];
residuals_h = [];
residuals_n = [];
residuals_e = [];
residuals_m_h = cell(1, 1024*1024); %cell to store residuals catagorized by number of observatiosn
residuals_m_n = residuals_m_h;
residuals_m_e = residuals_m_n;
sigma_0_sqrd_top_m_h = cell(1, 1024*1024); %... similar to the above...
sigma_0_sqrd_top_m_n = sigma_0_sqrd_top_m_h;
sigma_0_sqrd_top_m_e = sigma_0_sqrd_top_m_h;
sigma_0_sqrd_bottom_m_h = cell(1, 1024*1024); %... ditto ...
sigma_0_sqrd_bottom_m_n = sigma_0_sqrd_bottom_m_h;
sigma_0_sqrd_bottom_m_e = sigma_0_sqrd_bottom_m_h;

output_string = "";
graphics_toolkit ("gnuplot"); %fixes crash with some Intel drivers
output_string_csv = "m, mu, var\n";

%Transformation
[EASTINGS, NORTHINGS, HEIGHTS] = cartesian_to_UTM32Eetrs89(XKOOR,YKOOR,ZKOOR);
%save('preresid.m')
%load('preresid.m')
% For each unique in GPSNR do:
for i = 1:length(Y);
  refnr = Y(i);
  gpsnr = GPSNR(I(i));
  
  is5Dpoint = 0;
  for j = 1:length(GPSNR_5D)
    if gpsnr{} == GPSNR_5D(j){}
      is5Dpoint = 1;
    end
  end

  %% Only use data from points designated in 5d_points file.
  if is5Dpoint
    
    if i < length(Y) %length(Y) = 642
      xkoords = XKOOR(I(i):(I(i+1)-1));
      ykoords = YKOOR(I(i):(I(i+1)-1));
      zkoords = ZKOOR(I(i):(I(i+1)-1));
      epochs = EPOCH(I(i):(I(i+1)-1));
      eastings = EASTINGS(I(i):(I(i+1)-1));
      northings = NORTHINGS(I(i):(I(i+1)-1));
      heights = HEIGHTS(I(i):(I(i+1)-1));
    else %Robusthed
      xkoords = XKOOR(I(i):end);
      ykoords = YKOOR(I(i):end);
      zkoords = ZKOOR(I(i):end);
      epochs = EPOCH(I(i):end);
      eastings = EASTINGS(I(i):end);
      northings = NORTHINGS(I(i):end);
      heights = HEIGHTS(I(i):end);
    end
    
    
    
    %% Convert epoch yyyymmdd to datenums (since Jan 1, zear ZERO)
    epochs = datenum(epochs,'yyyymmdd');
    num_measurements = length(epochs);
      
    %% Sorting 
    [epochs,Ind] = sort(epochs);
    xkoords = xkoords(Ind);
    ykoords = ykoords(Ind);
    zkoords = zkoords(Ind);
    eastings = eastings(Ind);
    northings = northings(Ind);
    heights = heights(Ind);
 
    %% Binning (see binning.m for details)
    if do_binning
      [eastings, ~] = binning(eastings,epochs,binsize);
      [northings, ~] = binning(northings,epochs,binsize);
      [heights, ~] = binning(heights,epochs,binsize);
      [xkoords, ~] = binning(xkoords,epochs,binsize);
      [ykoords, ~] = binning(ykoords,epochs,binsize);
      [zkoords, epochs] = binning(zkoords,epochs,binsize);
    end
    num_measurements_binned = length(epochs);
    
    %% Transformation (moved outside loop)
    %[eastings, northings, heights] = cartesian_to_UTM32Eetrs89(xkoords,ykoords,zkoords);
    
    %Convert from m/day to mm/year
    epochs = epochs./365.25;
    heights = 1000.*heights;
    northings = 1000.*northings;
    eastings = 1000.*eastings;

    %% Means
    epochs0 = mean(epochs);
    heights0 = mean(heights);
    northings0 = mean(northings);
    eastings0 = mean(eastings);
    
    %% Remove means from data  
    epochs = epochs - epochs0;
    heights = heights - heights0;
    northings = northings - northings0;
    eastings = eastings - eastings0;
    
    m = length(epochs);

    % Only save residuals if at least "Nmin" measurements
    if m >= Nmin;
      % Fit using custom 'linreg' function
      [bh_custom, statsh_custom] = linreg(epochs, heights, alpha);
      [bn_custom, statsn_custom] = linreg(epochs, northings, alpha);
      [be_custom, statse_custom] = linreg(epochs, eastings, alpha);

      %Save to total residual vector [mm] (får residualerne fra linreg, statsh_custom.resid = stats i linreg med undervariable, som defineret ovenfor, lin 161).
      residuals_h = [residuals_h; statsh_custom.resid]; %residuals_h etc. er initialiseret men ikke defineret endnu. (l. 55).
      residuals_n = [residuals_n; statsn_custom.resid];
      residuals_e = [residuals_e; statse_custom.resid];
      
      %Save to 'm' residual vector (cell)
      residuals_m_h{m} = [residuals_m_h{m} ; statsh_custom.resid];
      residuals_m_n{m} = [residuals_m_n{m} ; statsn_custom.resid];
      residuals_m_e{m} = [residuals_m_e{m} ; statse_custom.resid];

      % For calculating sigma0_all (Pooled Standard Deviation)
      %Heights;
      sigma_0_sqrd_top_h = [sigma_0_sqrd_top_h; statsh_custom.sigma_0_sqrd_top]; %Definerer top og bottom for hver residual vektor, 2 for hver, 4 i alt.
      sigma_0_sqrd_bottom_h = [sigma_0_sqrd_bottom_h; statsh_custom.sigma_0_sqrd_bottom]; %stats.sigma_0_sqrd_bottom og top kommer begge fra linreg, og er "The corrected sample standard deviation".
      sigma_0_sqrd_top_m_h{m} = [sigma_0_sqrd_top_m_h{m}; statsh_custom.sigma_0_sqrd_top]; %De er præcis de samme stats værdier som bliver læst ind i begge.
      sigma_0_sqrd_bottom_m_h{m} = [sigma_0_sqrd_bottom_m_h{m}; statsh_custom.sigma_0_sqrd_bottom];
      
      %Northings;
      sigma_0_sqrd_top_n = [sigma_0_sqrd_top_n; statsn_custom.sigma_0_sqrd_top];
      sigma_0_sqrd_bottom_n = [sigma_0_sqrd_bottom_n; statsn_custom.sigma_0_sqrd_bottom];
      sigma_0_sqrd_top_m_n{m} = [sigma_0_sqrd_top_m_n{m}; statsn_custom.sigma_0_sqrd_top];
      sigma_0_sqrd_bottom_m_n{m} = [sigma_0_sqrd_bottom_m_n{m}; statsn_custom.sigma_0_sqrd_bottom];
      
      %Eastings;
      sigma_0_sqrd_top_e = [sigma_0_sqrd_top_e; statse_custom.sigma_0_sqrd_top];
      sigma_0_sqrd_bottom_e = [sigma_0_sqrd_bottom_e; statse_custom.sigma_0_sqrd_bottom];
      sigma_0_sqrd_top_m_e{m} = [sigma_0_sqrd_top_m_e{m}; statse_custom.sigma_0_sqrd_top];
      sigma_0_sqrd_bottom_m_e{m} = [sigma_0_sqrd_bottom_m_e{m}; statse_custom.sigma_0_sqrd_bottom];
    end
  end
end
%Done with For each unique in GPSNR do, 170 for loop.

%Calcualte sigma0_all (Pooled standard deviation)
%sigma0_all = sqrt(sum(sigma_0_sqrd_top)./sum(sigma_0_sqrd_bottom));

%Standard deviation pooled for hvert enkelt station (ikke sikker på man kan gøre det sådan dog)
sigma0_ind_h = sqrt(sigma_0_sqrd_top_h./sigma_0_sqrd_bottom_h); %Hvorfor summen i nævneren? Nævneren er antal frihedsgrader, som kun kan være en scalar.
sigma0_ind_n = sqrt(sigma_0_sqrd_top_n./sigma_0_sqrd_bottom_n);
sigma0_ind_e = sqrt(sigma_0_sqrd_top_e./sigma_0_sqrd_bottom_e);

%Dette er så pooled standard deviation for samtlige punkter, for northings, heights og eastings.
sigma0_all_h = sqrt(sum(sigma_0_sqrd_top_h)./sum(sigma_0_sqrd_bottom_h)); %Hvorfor summen i nævneren? Nævneren er antal frihedsgrader, som kun kan være en scalar.
sigma0_all_n = sqrt(sum(sigma_0_sqrd_top_n)./sum(sigma_0_sqrd_bottom_n));
sigma0_all_e = sqrt(sum(sigma_0_sqrd_top_e)./sum(sigma_0_sqrd_bottom_e));

m_vec = find(~cellfun('isempty',residuals_m_h));
numbins = length(m_vec);

%% Remove empty cells
residuals_m_h = residuals_m_h(~cellfun('isempty',residuals_m_h));
sigma_0_sqrd_top_m_h = sigma_0_sqrd_top_m_h(~cellfun('isempty',sigma_0_sqrd_top_m_h));
sigma_0_sqrd_bottom_m_h = sigma_0_sqrd_bottom_m_h(~cellfun('isempty',sigma_0_sqrd_bottom_m_h));

residuals_m_n = residuals_m_n(~cellfun('isempty',residuals_m_n));
sigma_0_sqrd_top_m_n = sigma_0_sqrd_top_m_n(~cellfun('isempty',sigma_0_sqrd_top_m_n));
sigma_0_sqrd_bottom_m_n = sigma_0_sqrd_bottom_m_n(~cellfun('isempty',sigma_0_sqrd_bottom_m_n));

residuals_m_e = residuals_m_e(~cellfun('isempty',residuals_m_e));
sigma_0_sqrd_top_m_e = sigma_0_sqrd_top_m_e(~cellfun('isempty',sigma_0_sqrd_top_m_e));
sigma_0_sqrd_bottom_m_e = sigma_0_sqrd_bottom_m_e(~cellfun('isempty',sigma_0_sqrd_bottom_m_e));

%Experimental code input by Jonathan, 10-02-2020:
tic

if histcalc == 1
figure(1)
bar(sigma0_ind_h)
title('Height standard deviation histogram')
xlabel 'Station nr.'
ylabel 'Standard deviation'
figure(2)
bar(sigma0_ind_n)
title('Northings standard deviation histogram')
xlabel 'Station nr.'
ylabel 'Standard deviation'
figure(3)
bar(sigma0_ind_e)
title('Eastings standard deviation histogram')
xlabel 'Station nr.'
ylabel 'Standard deviation'
end
resT = toc;
sprintf('Histograms calculated in %i seconds\n',resT)

%% Figures
if residuals_figures %Default 0, da residuals funktionen ikke er implementeret i Octave endnu. Måske virker det i Matlab? Tjek.
  % For each M
  for i = 1:length(m_vec) %i == 0???
    fig = figure;%('PaperPositionMode','auto');
    set(gcf,'PaperType','A4', ...
         'paperOrientation', 'landscape','PaperPositionMode','auto', ...
         'paperunits','CENTIMETERS', ...
         'PaperPosition',[.63, .63, 28.41, 19.72]); 
##    if i == 1
      %[h,binCenters] = hist(residuals,numbins);
      %Experimental
      C=residuals_m_h; 
    maxLengthCell=max(cellfun('size',C,1));  %finding the longest vector in the cell array
    for i=1:length(C)
        for j=cellfun('size',C(i),2)+1:maxLengthCell
             C{i}(j)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
        end
    end
    A=cell2mat(C); %A is your matrix
      [h,binCenters] = hist(A,numbins);
      %
##      mu_residuals = mean(residuals);
##      var_residuals = var(residuals);
      mu_residuals = mean(A);
      var_residuals = var(A);

      %sigma0_m = sigma0_all;
      titlestring = ['m >= ' num2str(Nmin) " \n"];
            
      output_string_csv = [output_string_csv, "all,",...
                           num2str(mu_residuals),", ",...
                           num2str(var_residuals), "\n"];
##    else 
##      [h,binCenters] = hist(residuals_m{i},numbins);
##      titlestring = ['m = ' num2str(m_vec(i)) " \n"];
##      mu_residuals = mean(residuals_m{i});
##      var_residuals = var(residuals_m{i});
##      sigma0_m = sqrt(sum(sigma_0_sqrd_top_m{i})./sum(sigma_0_sqrd_bottom_m{i}));
##      
##      output_string_csv = [output_string_csv, int2str(m_vec(i)),", ",...
##                           num2str(mu_residuals),", ",...
##                           num2str(var_residuals),"\n"];
##    end

    %Normalize histogram
    h = h/(sum(h)); %Det er ikke sådan du normaliserer?

    bar(binCenters,h, 'DisplayName', 'Residuals'); 
    hold on;
        
    x = linspace(min(binCenters),max(binCenters),100);
##    y = normpdf(x,mu_residuals,sqrt(var_residuals))*(binCenters(2)-binCenters(1));
      y = normpdf(x)*(binCenters(2)-binCenters(1));
##    y2 = normpdf(x,mu_residuals,sigma0_m)*(binCenters(2)-binCenters(1));
    plot(x,y,'k','linewidth',2);
    plot(x,y2,'g','linewidth',2);

    titlespec = "\mu = %f5 mm, \sigma = %f5 mm";
    title([titlestring sprintf(titlespec,mu_residuals,sqrt(var_residuals))]);
    
    ylabel('Normalized counts');
    xlabel('Residual [mm]');
    if residual_figures == 2
      if i == 1
        filename = sprintf('all.%s',figure_format);
      else
        filename = sprintf('m%i.%s',m_vec(i),figure_format);
      end
      saveas(fig,[resid_fig_folder filename],figure_format)
      if closeFigures
        close(fig);
      end
    end
  end
end
fprintf('Calculations done.\n');
% Save residuals
%save("outputs\\residuals.mat","residuals","Nmin");

fileID = fopen(filename_output_csv,'w');
fprintf(fileID,output_string_csv);
fclose(fileID);

fprintf('.CSV File written.\n')