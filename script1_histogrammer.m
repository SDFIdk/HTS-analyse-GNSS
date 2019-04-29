%SCRIPT1_HISTOGRAMMER
%
%
%
% oldjo@sdfe.dk, November 2017



% Filename
filename_input = 'inputs\\5D_tidsserier_20170504_trimmed.csv';
% Note: this must be a comma seperated file with the following columns
% REFNR,GPSNR,XKOOR,YKOOR,ZKOOR,XRMS,YRMS,ZRMS,JNR_BSIDE,SYS,EPOCH,IN_DATE
filename_5d_points = 'inputs\\5d_points_sorted.csv';
uplifts_filename = 'inputs\\altgps_uplift';

% output filenames
filename_output_csv = 'residuals\\residuals.csv';
fig_folder = 'residuals\\';

% (1) The epoch will be in the middle of the interval of measurements
middle_epoch = 1;
% (0) Set epoch to evaluate manually. Format yyyymmdd
epoch_eval = '20070101'; 


% Alpha
alpha = 0.05;

%Bining
binning = 0; %1: on, 0: off.
binsize = 14; % Binsize in days

% Generate figures (0: no, 1: display, 2: save and display).
figures = 2;
numbins = 25; %number of bins in histograms
closeFigures = 0;
Nmin = 3;

% Figure format: (png, pdf, jpg... etc.)
figure_format = 'pdf'; 

% Load 'optim' package for 'polyconf'
%pkg load optim

% Load 'statistics' package for 'regress'
%pkg load statistics

% Formatstring (change this if the input files get different fields.
T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
[REFNR, GPSNR, XKOOR, YKOOR, ZKOOR, SYS, EPOCH] =...
textread(filename_input, T, "delimiter", ",", "headerlines", 1);
[GPSNR_5D] = textread(filename_5d_points, "%s");

% Find Unique entries
[Y, I, J] = unique (REFNR, "first");

%Initialize strings etc.
residuals = [];
residuals_m = cell(1, 1024*1024);
output_string = "";
graphics_toolkit ("gnuplot"); %fixes crash with some Intel drivers

output_string_csv = "m, mu, var\n";

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
  
  %% Specify gps name below to only look at that specific one.
  %if (gpsnr{} == "SEED")
  
  %% Only use data from points designated in 5d_points file.
  if is5Dpoint
    
  if i < length(Y)
    xkoords = XKOOR(I(i):(I(i+1)-1));
    ykoords = YKOOR(I(i):(I(i+1)-1));
    zkoords = ZKOOR(I(i):(I(i+1)-1));
    epochs = EPOCH(I(i):(I(i+1)-1));
  else
    xkoords = XKOOR(I(i):end);
    ykoords = YKOOR(I(i):end);
    zkoords = ZKOOR(I(i):end);
    epochs = EPOCH(I(i):end);
  end
  
  sys = SYS(I(i));
  
  %% Convert epoch yyyymmdd to datenums (since Jan 1, zear ZERO)
  %epochs = datenum(datevec(epochs,'yyyymmdd'));
  epochs = datenum(epochs,'yyyymmdd');
  num_measurements = length(epochs);
     
  %% Sorting 
  [epochs,Ind] = sort(epochs);
  xkoords = xkoords(Ind);
  ykoords = ykoords(Ind);
  zkoords = zkoords(Ind);
 
  %% Binning
  if binning
    [xkoords, ~] = binning(xkoords,epochs,binsize);
    [ykoords, ~] = binning(ykoords,epochs,binsize);
    [zkoords, epochs] = binning(zkoords,epochs,binsize);
  end
  
  %% Transformation
  [eastings, northings, heights] = cartesian_to_UTM32Eetrs89(xkoords,ykoords,zkoords);
  
  %% Means
  epochs0 = mean(epochs);
  xkoords0 = mean(xkoords);
  ykoords0 = mean(ykoords);
  zkoords0 = mean(zkoords);
  eastings0 = mean(eastings);
  northings0 = mean(northings);
  heights0 = mean(heights);
  
  %% Remove means from data  
  epochs = epochs - epochs0;
  x = xkoords - xkoords0;
  y = ykoords - ykoords0;
  z = zkoords - zkoords0;
  eastings = eastings - eastings0;
  northings = northings - northings0;
  heights = heights - heights0;
  
  %Set epoch eval to middle epoch
  if middle_epoch
    epoch_eval = datestr(epochs0,'yyyymmdd');
    epoch_eval_num = datenum(datevec(epoch_eval,'yyyymmdd'));
  end
  
  %% Fit to polynomial by least squares in each dimenstion
  %[px,Sx] = polyfit(epochs, x, 1);
  %[py,Sy] = polyfit(epochs, y, 1);
  %[pz,Sz] = polyfit(epochs, z, 1);
 
  
  % Fit a third time using custom 'linreg' function
  %[bh_custom, statsh_custom] = linreg(epochs, heights, alpha, uplifts(i));
  [bh_custom, statsh_custom] = linreg(epochs, heights, alpha);
  
  m = length(epochs);
  % Only save residuals if at least "Nmin" measurements
  if m >= Nmin;
    %Save to total residual vector [mm]
    residuals = [residuals; 1000.*statsh_custom.resid];
    %Save to 'm' residual vector
    
    residuals_m{m} = [residuals_m{m} ; 1000.*statsh_custom.resid];

  end
  
  %Vertical velocity mm/anno (from m/day);
  %v_h = bh_custom(2)*1000*365.25;
  
end
end

m_vec = find(~cellfun('isempty',residuals_m));

%% Remove empty cells
residuals_m = residuals_m(~cellfun('isempty',residuals_m));

%% Figures
if figures
  % For each M
  for i = 0:length(m_vec) 
    fig = figure;%('PaperPositionMode','auto');
    set(gcf,'PaperType','A4', ...
         'paperOrientation', 'landscape','PaperPositionMode','auto', ...
         'paperunits','CENTIMETERS', ...
         'PaperPosition',[.63, .63, 28.41, 19.72]); 
    if i == 0
      [h,binCenters] = hist(residuals,numbins);
      mu_residuals = mean(residuals);
      var_residuals = var(residuals);
      titlestring = ['m >= ' num2str(Nmin) " \n"];
      
      output_string_csv = [output_string_csv, "all,",...
                           num2str(mu_residuals),", ",...
                           num2str(var_residuals), "\n"];
    else 
      [h,binCenters] = hist(residuals_m{i},numbins);
      titlestring = ['m = ' num2str(m_vec(i)) " \n"];
      mu_residuals = mean(residuals_m{i});
      var_residuals = var(residuals_m{i});
      
      output_string_csv = [output_string_csv, int2str(m_vec(i)),", ",...
                           num2str(mu_residuals),", ",...
                           num2str(var_residuals),"\n"];
    end

    %Normalize histogram
    h = h/(sum(h)); 

    bar(binCenters,h, 'DisplayName', 'Residuals'); 
    hold on;
    
    
    x = linspace(min(binCenters),max(binCenters),100);
    y = normpdf(x,mu_residuals,sqrt(var_residuals))*(binCenters(2)-binCenters(1));
    plot(x,y,'k','linewidth',2);

    titlespec = "\mu = %f5 mm, \sigma = %f5 mm";
    title([titlestring sprintf(titlespec,mu_residuals,sqrt(var_residuals))]);
    
    ylabel('Normalized counts');
    xlabel('Residual [mm]');
    if figures == 2
      if i == 0
        filename = sprintf('all.%s',figure_format);
      else
        filename = sprintf('m%i.%s',m_vec(i),figure_format);
      end
      saveas(fig,[fig_folder filename],figure_format)
      if closeFigures
        close(fig);
      end
    end
  end
end
% Save residuals

save("outputs\\residuals.mat","residuals","Nmin");
fprintf('Calculations done.\n');

fileID = fopen(filename_output_csv,'w');
fprintf(fileID,output_string_csv);
fclose(fileID);

fprintf('.CSV File written.\n')