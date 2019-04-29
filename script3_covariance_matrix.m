%SCRIPT2_SIGMA_0.m
%
% Script to get covariance matrices from the regression of "5D-net" data
%
% oldjo@sdfe.dk, November 2017
%

% If sigma_0 is not defined
try (sigma_0_all) catch sigma_0_all = 0.00401401223635878; %Found in script2
end

graphics_toolkit ("gnuplot"); %fixes crash with some Intel drivers

% Filename
filename_input = 'inputs\\5D_tidsserier_20170504_trimmed.csv';
% Note: this must be a comma seperated file with the following columns
% REFNR,GPSNR,XKOOR,YKOOR,ZKOOR,XRMS,YRMS,ZRMS,JNR_BSIDE,SYS,EPOCH,IN_DATE
filename_5d_points = 'inputs\\5d_points_sorted.csv';
uplifts_filename = 'inputs\\altgps_uplift';

% output filenames
filename_output_csv = 'outputs\\sigma_0.csv';
filename_output_koords = 'outputs\\coordinates.txt';
fig_folder = 'figures\\';

% Alpha
alpha = 0.05;

% Binsize in days
binning = 0;
binsize = 14;

% Generate figures (0: no, 1: display, 2: save and display).
figures = 2;
closeFigures = 0; %1;
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

%Init strings etc.
output_string_csv = "GPSNR,sigma_0 [mm],Number of Measurements\n";
residuals = [];
sigma_0_sqrd_top = [];
sigma_0_sqrd_bottom = [];

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
    %eastings0 = mean(eastings);
    %northings0 = mean(northings);
    heights0 = mean(heights);
    
    epochs = epochs - epochs0;
    x = xkoords - xkoords0;
    y = ykoords - ykoords0;
    z = zkoords - zkoords0;
    %eastings = eastings - eastings0;
    %northings = northings - northings0;
    heights = heights - heights0;
   
    % Fit a third time using custom 'linreg' function
    [bx_custom, statsx_custom] = linreg(epochs, z, alpha);%, uplifts(i));
    [by_custom, statsy_custom] = linreg(epochs, y, alpha);%, uplifts(i));
    [bz_custom, statsz_custom] = linreg(epochs, z, alpha);%, uplifts(i));
    [bh_custom, statsh_custom] = linreg(epochs, heights, alpha, 0 ,sigma_0_all);%, uplifts(i));
    statsh_custom.SIGMA3
    statsh_custom.SIGMA4
    % Only save residuals if at least "Nmin" measurements
    if length(epochs) >= Nmin;
      residuals = [residuals; statsh_custom.resid];
      sigma_0_sqrd_top = [sigma_0_sqrd_top; statsh_custom.sigma_0_sqrd_top];
      sigma_0_sqrd_bottom = [sigma_0_sqrd_bottom; statsh_custom.sigma_0_sqrd_bottom];
    end
  
  
    %Vertical velocity mm/anno (from m/day);
    v_h = bh_custom(2)*1000*365.25;
 
    % Don't record from stations with only one measurement
  if num_measurements >= Nmin
 
    formatSpec_csv = "%s,%.5f,%i\n";
    output_string_csv = [output_string_csv sprintf(formatSpec_csv,...
                         char(gpsnr), 1000.*statsh_custom.sigma_0,...
                         length(epochs))];
    end
  end
end

% Save residuals
save("outputs\\residuals.mat","residuals","Nmin");

%sigma_0 = sqrt((residuals'*eye(length(residuals))*residuals)./...
%               (length(residuals)-2));

sigma_0 = sqrt(sum(sigma_0_sqrd_top)/sum(sigma_0_sqrd_bottom));
               
               
               
output_string_csv = [output_string_csv sprintf(formatSpec_csv,...
                     "ALL_5D", sigma_0.*1000,...
                     length(residuals))];   


fileID = fopen(filename_output_csv,'w');
fprintf(fileID,output_string_csv);
fclose(fileID);

fprintf('.CSV File written.\n')
