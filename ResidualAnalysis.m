%RESIDUALANALYSIS
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

% alpha
if ~exist('alpha')
  alpha = 0.05;
else
  fprintf(['alpha already defined as ' num2str(alpha) '\n']);
end

% Binning (pre binning of observations, see 'binning.m' for the inner workings)
do_binning = 1; %1: on, 0: off.
binsize = 14; % Binsize in days

% Generate figures (0: no, 1: display, 2: save and display).
residuals_figures = 0; 
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

      %Save to total residual vector [mm] (gets residuals from linreg, statsh_custom.resid = stats in linreg with subvariables, as defined above).
      residuals_h = [residuals_h; statsh_custom.resid]; %residuals_h etc. is initialized but not defined yet (1. 55).
      residuals_n = [residuals_n; statsn_custom.resid];
      residuals_e = [residuals_e; statse_custom.resid];
      
      %Save to 'm' residual vector (cell)
      residuals_m_h{m} = [residuals_m_h{m} ; statsh_custom.resid];
      residuals_m_n{m} = [residuals_m_n{m} ; statsn_custom.resid];
      residuals_m_e{m} = [residuals_m_e{m} ; statse_custom.resid];

      % For calculating sigma0_all (Pooled Standard Deviation)
      %Heights;
      sigma_0_sqrd_top_h = [sigma_0_sqrd_top_h; statsh_custom.sigma_0_sqrd_top]; %Defines top and bottom for each residual vector, 2 each, 4 in total
      sigma_0_sqrd_bottom_h = [sigma_0_sqrd_bottom_h; statsh_custom.sigma_0_sqrd_bottom]; %stats.sigma_0_sqrd_bottom and top both comes from linreg, and are "The corrected sample standard deviation".
      sigma_0_sqrd_top_m_h{m} = [sigma_0_sqrd_top_m_h{m}; statsh_custom.sigma_0_sqrd_top]; %It is exactly the same stats values that are read into both.
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
%Done with "For each unique in GPSNR do", of the 170 for loop.

%Calculate sigma0_all (Pooled standard deviation)
%sigma0_all = sqrt(sum(sigma_0_sqrd_top)./sum(sigma_0_sqrd_bottom));

%Standard deviation pooled for each individual station
sigma0_ind_h = sqrt(sigma_0_sqrd_top_h./sigma_0_sqrd_bottom_h);
sigma0_ind_n = sqrt(sigma_0_sqrd_top_n./sigma_0_sqrd_bottom_n);
sigma0_ind_e = sqrt(sigma_0_sqrd_top_e./sigma_0_sqrd_bottom_e);

%This is then the pooled standard deviation for all points, for northings, heights and eastings.
sigma0_all_h = sqrt(sum(sigma_0_sqrd_top_h)./sum(sigma_0_sqrd_bottom_h)); 
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

%By Jonathan, 10-02-2020:
tic
% Generates histograms for Heights, Northings and Eastings standard deviations.
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

fprintf('Calculations done.\n');

fileID = fopen(filename_output_csv,'w');
fprintf(fileID,output_string_csv);
fclose(fileID);

fprintf('.CSV File written.\n')