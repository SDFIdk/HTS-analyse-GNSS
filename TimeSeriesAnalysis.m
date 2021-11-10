% TIMESERIESANALYSIS
%
% Inputs:
% filename_input: csv file with gpsnr, refnr, X,Y,Z coordinates and
%                 YYYYmmdd epoch (see the format string 'T' below)
%						Note: this must be a comma seperated file with the following columns
% 						    REFNR,GPSNR,XKOOR,YKOOR,ZKOOR,XRMS,YRMS,ZRMS,JNR_BSIDE,SYS,EPOCH,IN_DATE
%  						    File must be sorted by REFNR in ascending order.
%   
%	filename_5d_points: column of four letter gnss stations that are 5d points
%   
%	filename_uplifts: see readUplifts.m for file specification
%
% Most settings are on/off by using 1 or 0;
%
% oldjo@sdfe.dk, November 2017
%
%%############################################################################ 
%% INPUT FILENAMES:
clear;clc;
%% Data file 
filename_input = 'inputs\\5D_tidsserier_20190820.csv';

%% File identifying the 5D points:
filename_5d_points = 'inputs\\5d_points_sorted.csv'; 

%% File with uplift values from DTU in mm/year
filename_uplifts = 'inputs\\altgps_uplift';

%%############################################################################ 
%% SIGNIFICANCE LEVEL:
alpha = 0.05;

%% SIGMA_0_ALL e.g. Pooled Standard Deviation, leave undefined to recalculate.
%sigma0_all = 4.01401;

%%############################################################################ 
%% CALCULATE AND OUTPUT FOR EASTINGS AND NORTHINGS:
%% 1 = yes, 0 = only height is calculated
outputENH = 0; 

%%############################################################################
%% EXCLUSION OPTION:
%% Select which data to use
%%
%% 1: Use all points with 2 or more measurements (or as defined in min_points) 
%%
%% 2: Only use data from points designated in 5d_points file
%% if (is5Dpoint) && (num_measurements_binned > min_points;))
%%
%% 3: Only use data from points designated in 5d_points file and Uplift file
%% if (is5Dpoint && ~isnan(uplifts(i)) && (num_measurements_binned > min_points))
%%
%% 4: Only calculate for given GNSS station
%%
exclusionOption = 1;

%% (When exclusionOption=4, only calculate for this GNSS/5D station)
exclusiveGPSNR = 'VAEG';

%%############################################################################
%% MINIMUM NUMBER OF OBSERVATIONS:
%% Set the minimum number of observation for station to get included in analysis.
%% min_points=2 requires three observations per station.
min_points = 1; 

%%############################################################################
%% BINNING:
%% 1 = on, 0 = off
do_binning = 1; 
%% Binsize in days
binsize = 200; 

%%############################################################################ 
%% FIGURES: 
%% Choose either figures or strength_figures (statistical strength), can't do both (bug due to octave update).
%% 0 = no figures, 1 = show, 2 = save, 3 = show and save
figures = 2; 
strength_figures = 0;

%%############################################################################
%% MISCELLANEOUS SETTINGS: 
%%
%% Combining output plots:
%% IMPORTANT: The two "combine.bat" files in folders "/figures/" and subfolder
%%            "/figures/strength/" use GhostView
%%            Please edit the path in the .bat files or set below to 0
combine_plots = 0;

%% Mark the middle epoch value on plots?
plot_middle_epoch = 1;

%% [mm], define the length of the y-axis, set to 0 for auto
y_limit = 100; 

%% prevents too many open windows, closes them immidiately (=1)
closeFigures = 1; 

%% print out the gps week for the observation (=0)
show_gpsweek = 1;

%% Output folder, format and filenames:
%% Figure format: 'pdf', 'png' etc.
fig_folder = 'figures\\';
strength_fig_folder = 'figures\\strength\\';
figure_format = 'pdf'; 
filename_output_csv = 'outputs\\out.csv';
filename_middle_epoch = 'outputs\\middle_epoch.csv';
filename_output_ENH = 'outputs\\out_ENH.csv';
filename_output_H = 'outputs\\out_H.csv';
filename_output_geoidplot = 'outputs\\out_geoidplot.csv';

%%##############################################################################
%%##############################################################################
%%
%%             PROGRAM BEGINS
%%
%%##############################################################################
%%##############################################################################

pkg load statistics

% Formatstring (change this if the input files get different fields. %* means that the field is ignored in textread.
T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
[REFNR, GPSNR, XKOOR, YKOOR, ZKOOR, SYS, EPOCH] =...
textread(filename_input, T, "delimiter", ",", "headerlines", 1);
[GPSNR_5D] = textread(filename_5d_points, "%s");

graphics_toolkit ("gnuplot"); %fixes crash with some Intel drivers
more off

%Run residual Analysis script if no sigma0_all is specified:
if ~exist('sigma0_all')
  fprintf('Running ResidualAnalysis script....\n')
  residualanalysistimer = tic;
  ResidualAnalysis;
  fprintf(['ResidualAnalysis script completed in '...
            num2str(toc(residualanalysistimer)) ' seconds\n']);
else
  sigma0_all_h = sigma0_all;
  sigma0_all_E = sigma0_all;
  sigma0_all_N = sigma0_all;
  fprintf(['Pooled Standard Deviation already defined,'...
           ' skipping residual analysis...\n']);
end

%Transformation
[EASTINGS, NORTHINGS, HEIGHTS] = cartesian_to_UTM32Eetrs89(XKOOR,YKOOR,ZKOOR);

% Find Unique entries
[Y, I, J] = unique(REFNR, "first");

%Initialize strings etc.
output_string_geoidplot_csv = ['GPSNR, REFNR, N, N_binned,'...  %string, int, int, int  
                         'MiddleEpochDate,'...          %string
                         'MiddleEpochE,'...             %float
                         'MiddleEpochN,'...             %float
                         'MiddleEpochHeight,'...        %float
                         'dh/dt [mm/Year],'...          %float
                         'Uplift DTU [mm/Year],'...     %float
                         'delta dh/dt [mm/Year] (fit-DTU),'...
                         'R^2 (h),'...                  %float ....
                         'StdEst observation (dh/dt),'...
                         'StdPooled observation (h),'...
                         'StdEst (dh/dt),'...
                         'StdPooled (dh/dt),'...
                         't-test (h), z-test (h)\n'];   %int ,int          

output_string_H_csv = ['GPSNR, REFNR, N, N_binned,'...  %string, int, int, int  
                         'MiddleEpochDate,'...          %string
                         'MiddleEpochGPSweek,'...       %int
                         'MiddleEpochE,'...             %float
                         'MiddleEpochN,'...             %float
                         'MiddleEpochHeight,'...        %float
                         'StdMiddleEpochHeightEst,'...  %float
                         'StdMiddleEpochHeightPooled,'...
                         'dh/dt [mm/Year],'...          %float
                         'Uplift DTU [mm/Year],'...     %float
                         'delta dh/dt [mm/Year] (fit-DTU),'...
                         'R^2 (h),'...                  %float ....
                         'StdEst observation (dh/dt),'...
                         'StdPooled observation (h),'...
                         'StdEst (dh/dt),'...
                         'ConfIntervalEstLow (dh/dt),'...
                         'ConfIntervalEstHigh (dh/dt),'...
                         'StdPooled (dh/dt),'...
                         'ConfIntervalLow (dh/dt),'...
                         'ConfIntervalHigh (dh/dt),'...
                         't-test (h), z-test (h)\n'];   %int ,int           
                     
output_string_ENH_csv = ['GPSNR, REFNR, N, N_binned,'...
                         'MiddleEpochDate,'...
                         'MiddleEpochGPSweek,'...
                         'MiddleEpochE,'...
                         'StdMiddleEpochEest,'...
                         'StdMiddleEpochEpooled,'...
                         'MiddleEpochN,'...
                         'StdMiddleEpochNest,'...
                         'StdMiddleEpochNpooled,'... 
                         'MiddleEpochHeight,'...
                         'StdMiddleEpochHeightEst,'...
                         'StdMiddleEpochHeightPooled,'...
                         'dh/dt [mm/Year],'...
                         'Uplift DTU [mm/Year],'...
                         'delta dh/dt [mm/Year] (fit-DTU),'...
                         'R^2 (h),'...
                         'StdEst observation (dh/dt),'...
                         'StdPooled observation (h),'...
                         'StdEst (dh/dt),'...
                         'ConfIntervalEstLow (dh/dt),'...
                         'ConfIntervalEstHigh (dh/dt),'...
                         'StdPooled (dh/dt),'...
                         'ConfIntervalLow (dh/dt),'...
                         'ConfIntervalHigh (dh/dt),'...
                         't-test (h), z-test (h),'...
                         'dE/dt [mm/year],'...
                         'R^2 (E),'...
                         'EstVar (E),'...
                         'StdEst (dE/dt),'...
                         'ConfIntervalEstLow (E),'...
                         'ConfIntervalEstHigh (E),'...
                         'Std_pooled (dE/dt),'...
                         'ConfIntervalLow (E),'...
                         'ConfIntervalHigh (E),'...
                         't-test (E), z-test (E),'...
                         'dN/dt [mm/year],'...
                         'R^2 (N),'...
                         'EstVar (N),'...
                         'StdEst (N),'...
                         'ConfIntervalEstLow (N),'...
                         'ConfIntervalEstHigh (N),'...
                         'Std_pooled (dN/dt),'...
                         'ConfIntervalLow (N),'...
                         'ConfIntervalHigh (N),'...
                         't-test (N), z-test (N)\n'];

formatspec_geoidplot = ['%s, %i, %i, %i, %s, %.5f, %.5f,'...
                  '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'...
                  '%.5f, %.5f, %.5f, %i, %i\n'];                  
                             
formatspec_H = ['%s, %i, %i, %i, %s, %i, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'...
                  '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'...
                  '%.5f, %.5f, %.5f, %i, %i\n'];                  
                  
formatspec_ENH = ['%s, %i, %i, %i, %s, %i, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'...
                  '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'...
                  '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %i,'...
                  '%i, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'...
                  '%i, %i, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'...
                  '%.5f, %i, %i\n'];
output_middle_epoch = 'GPSNR, REFNR, Middle Epoch, X, Y, Z, N, dh/dt [mm/Year]\n';
formatspec_middle_epoch = '%s, %i, %s, %.5f, %.5f, %.5f, %i, %.2f\n';

%Read uplifts from file (see readUplifts.m for input file specification)
uplifts = readUplifts(filename_uplifts,GPSNR(I));
disp(binsize)
% For each unique in GPSNR do:
for i = 1:length(Y);
  refnr = Y(i);
  gpsnr = GPSNR(I(i));
    
  %Check if point is 5D point  
  is5Dpoint = 0;
  for j = 1:length(GPSNR_5D)
    if gpsnr{} == GPSNR_5D(j){}
      is5Dpoint = 1;
    end
  end

  uplift = uplifts(i);    
  if i < length(Y)
    xkoords = XKOOR(I(i):(I(i+1)-1));
    ykoords = YKOOR(I(i):(I(i+1)-1));
    zkoords = ZKOOR(I(i):(I(i+1)-1));
    epochs = EPOCH(I(i):(I(i+1)-1));
    eastings = EASTINGS(I(i):(I(i+1)-1));
    northings = NORTHINGS(I(i):(I(i+1)-1));
    heights = HEIGHTS(I(i):(I(i+1)-1));
  else
    xkoords = XKOOR(I(i):end);
    ykoords = YKOOR(I(i):end);
    zkoords = ZKOOR(I(i):end);
    epochs = EPOCH(I(i):end);
    eastings = EASTINGS(I(i):end);
    northings = NORTHINGS(I(i):end);
    heights = HEIGHTS(I(i):end);
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
  eastings = eastings(Ind);
  northings = northings(Ind);
  heights = heights(Ind);

  %% Binning (see binning.m for details)
  if do_binning == 1
    [eastings, ~] = binning(eastings,epochs,binsize);
    [northings, ~] = binning(northings,epochs,binsize);
    [heights, ~] = binning(heights,epochs,binsize);
    [xkoords, ~] = binning(xkoords,epochs,binsize);
    [ykoords, ~] = binning(ykoords,epochs,binsize);
    [zkoords, epochs] = binning(zkoords,epochs,binsize);
  end
  num_measurements_binned = length(epochs);

  %This switch controls which GPSNR are evaluated.
  switch exclusionOption
    case 1
      condition = (num_measurements_binned > min_points);
    case 2
      condition = ((is5Dpoint) && (num_measurements_binned > min_points));
    case 3
      condition = ((is5Dpoint && ~isnan(uplifts(i)) && (num_measurements_binned > min_points)));
    case 4
      condition = strcmp(gpsnr{},exclusiveGPSNR);
  end
  if condition
  
    %% Middle Epoch
    X0 = mean(xkoords); 
    Y0 = mean(ykoords);
    Z0 = mean(zkoords);
    eastings0 = mean(eastings);
    northings0 = mean(northings);
    epochs0 = mean(epochs);
    heights0 = mean(heights);
  
    %Convert from m/day to mm/year
    epochs = epochs./365.25;
    heights = 1000.*heights;
    eastings = 1000.*eastings;
    northings = 1000.*northings;
  
    fprintf('...');
      
    %Call linreg function
    [b, stats] = linreg(epochs, heights, alpha, uplift, sigma0_all_h);
    if outputENH
      [b_easting, stats_easting] = linreg(epochs, eastings, alpha, 0, sigma0_all_e);
      [b_northing, stats_northing] = linreg(epochs, northings, alpha, 0, sigma0_all_n);
    
      % Save values to output string
      output_string_ENH_csv = [output_string_ENH_csv sprintf(formatspec_ENH,...
             gpsnr{}, refnr, num_measurements, num_measurements_binned,...
		  datestr(epochs0,'dd-mm-YYYY'), gpsweek(epochs0), eastings0,...
      stats_easting.StdMiddleEpochEst,...
      stats_easting.StdMiddleEpochPooled,...
      northings0,...
      stats_northing.StdMiddleEpochEst,...
      stats_northing.StdMiddleEpochPooled,...
      heights0,...
      stats.StdMiddleEpochEst,...
      stats.StdMiddleEpochPooled,...
			b(2), uplift, b(2)-uplift, ...
			stats.Rsqr, stats.sigma_0, sigma0_all_h,...
			stats.std_unknown,...
      stats.confinterval_estimated(1),...
      stats.confinterval_estimated(2),...
			stats.std_known,...
      stats.confinterval_known(1),...
      stats.confinterval_known(2),...
			stats.t_test, stats.z_test_known,...
      b_easting(2),...
      stats_easting.Rsqr, stats_easting.sigma_0, sigma0_all_e,...
			stats_easting.std_unknown,...
      stats_easting.confinterval_estimated(1),...
      stats_easting.confinterval_estimated(2),...
			stats_easting.std_known,...
      stats_easting.confinterval_known(1),...
      stats_easting.confinterval_known(2),...
			stats_easting.t_test, stats_easting.z_test_known,...
      b_northing(2),...
      stats_northing.Rsqr, stats_northing.sigma_0, sigma0_all_n,...
			stats_northing.std_unknown,...
      stats_northing.confinterval_estimated(1),...
      stats_northing.confinterval_estimated(2),...
			stats_northing.std_known,...
      stats_northing.confinterval_known(1),...
      stats_northing.confinterval_known(2),...
			stats_northing.t_test, stats_northing.z_test_known)];
    else  
      output_string_H_csv = [output_string_H_csv sprintf(formatspec_H,...
             gpsnr{}, refnr, num_measurements, num_measurements_binned,...
		  datestr(epochs0,'dd-mm-YYYY'), gpsweek(epochs0),...
      eastings0,...
      northings0,...
      heights0,...
      stats.StdMiddleEpochEst,...
      stats.StdMiddleEpochPooled,...
			b(2), uplift, b(2)-uplift, ...
			stats.Rsqr, stats.sigma_0, sigma0_all_h,...
			stats.std_unknown,...
      stats.confinterval_estimated(1),...
      stats.confinterval_estimated(2),...
			stats.std_known,...
      stats.confinterval_known(1),...
      stats.confinterval_known(2),...
			stats.t_test, stats.z_test_known)];
    end

	  %Make list to print to the output file "middle_epoch.csv"
	  output_middle_epoch = [output_middle_epoch sprintf(formatspec_middle_epoch,...
                         gpsnr{}, refnr, datestr(epochs0,'YYYYmmdd'), X0, Y0, Z0, num_measurements, b(2))];
 
    output_string_geoidplot_csv = [output_string_geoidplot_csv sprintf(formatspec_geoidplot,...
    gpsnr{}, refnr, num_measurements, num_measurements_binned,...
		datestr(epochs0,'dd-mm-YYYY'),...
    eastings0,...
    northings0,...
    heights0,...
    b(2), uplift, b(2)-uplift, ...
		stats.Rsqr, stats.sigma_0, sigma0_all_h,...
		stats.std_unknown,...
    stats.std_known,...
    stats.t_test, stats.z_test_known)];
 
    fprintf(['\n gpsnr: ' gpsnr{} ' done! \n'])
       
  end
  
  %Only plot figures if there are 2 or more points
  if ((figures > 0) && condition)
  
    %Set title for plots if any.
    titlestring = ['gpsnr: ' gpsnr{} ', n_{measurements}: ',...
      num2str(num_measurements) ', n_{binned}: ',...
      num2str(num_measurements_binned) ];
    fig = linregplot(epochs,heights,alpha,uplift,sigma0_all_h,show_gpsweek,...
                     plot_middle_epoch,y_limit,titlestring);
    set(fig, 'Position', get(0, 'Screensize'));
       
    %Show figures
    if (figures == 1) || (figures == 3)
      set(fig, 'visible', 'on');
    end
  
    %Save figures    
    if (figures == 2) || (figures == 3)
      filename = sprintf('%s.%s',gpsnr{},figure_format);
      set(gcf,'PaperUnits','normalized');
      set(gcf,'PaperPosition', [0 0 1 1]);
      print(fig, [fig_folder filename]);
    end
    
    if closeFigures
      close(fig);
    end
  end
  if ((strength_figures > 0) && condition)
    %z = abs(b(2)-uplift)/stats.std_known;
    delta = linspace(-5,5,500);
    
    styrke_title_string = [gpsnr{} ', Statistical Power for Z-test at significance ' ...
    'level ' num2str(alpha) ', N = ' num2str(num_measurements) ', N_{binned} = ' ...
    num2str(num_measurements_binned)];
    
    %fprintf(['delta is ' num2str(delta)]);
    fprintf(['stats.z_crit is ' num2str(stats.z_crit) '\n']);
    fprintf(['stats.std_known is ' num2str(stats.std_known) '\n']);
    fprintf(['alpha is ' num2str(alpha) '\n']);
    
    fig = styrkeFunktionPlot(delta,stats.z_crit,stats.std_known,alpha,...
          styrke_title_string);     
    %Show strength curve figures
    if (strength_figures == 1) || (strength_figures == 3)
      set(fig, 'visible', 'on');
    end
    
    %Save strength curve figures    
    if (strength_figures == 2) || (strength_figures == 3)
      filename = sprintf('%s_%s.%s',gpsnr{},'strength',figure_format);
      set(gcf,'PaperUnits','normalized');
      set(gcf,'PaperPosition', [0 0 1 1]);
      print(fig, [strength_fig_folder filename ]);
    end
          
    if closeFigures
      close(fig);
    end
  end
end


%Combine pdf files into one
if combine_plots > 0
  if (figures == 2) || (figures == 3)
    cd 'figures'
    [~,output] = system("combine.bat");
    cd '..'
  end
  if (strength_figures == 2) || (strength_figures == 3)
    cd 'figures'
    cd 'strength'
      [~,output] = system("combine.bat");
    cd '..'
    cd '..'
  end
end
 

if outputENH
  fileID = fopen(filename_output_ENH,'w');
  fprintf(fileID,output_string_ENH_csv);
  fclose(fileID);
  fprintf('.CSV File written.\n')  
else
  fileID = fopen(filename_output_H,'w');
  fprintf(fileID,output_string_H_csv);
  fclose(fileID);
  fprintf('.CSV File written.\n')  
end

fileID = fopen(filename_middle_epoch,'w');
fprintf(fileID,output_middle_epoch);
fclose(fileID);
fprintf('.CSV File written.\n')

fileID = fopen(filename_output_geoidplot,'w');
fprintf(fileID,output_string_geoidplot_csv);
fclose(fileID);
fprintf('.CSV File written.\n')  