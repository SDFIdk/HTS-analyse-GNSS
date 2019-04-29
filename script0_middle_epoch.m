%SCRIPT0_MIDDLE_EPOCH
%
%
%
% oldjo@sdfe.dk, November 2017



% Filename
%filename_input = 'inputs\\5D_tidsserier_20170504_trimmed.csv';
filename_input = 'inputs\\HHLM.csv';

% Note: this must be a comma seperated file with the following columns
% REFNR,GPSNR,XKOOR,YKOOR,ZKOOR,XRMS,YRMS,ZRMS,JNR_BSIDE,SYS,EPOCH,IN_DATE
filename_5d_points = 'inputs\\5d_points_sorted.csv';
%uplifts_filename = 'inputs\\altgps_uplift';

% output filenames
filename_output_csv = 'middle_epoch\\out.csv';
fig_folder = 'middle_epoch\\';

% (1) The epoch will be in the middle of the interval of measurements
middle_epoch = 1;
% (0) Set epoch to evaluate manually. Format yyyymmdd
%epoch_eval = '20070101'; 

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

% Formatstring (change this if the input files get different fields.
T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
[REFNR, GPSNR, XKOOR, YKOOR, ZKOOR, SYS, EPOCH] =...
textread(filename_input, T, "delimiter", ",", "headerlines", 1);
[GPSNR_5D] = textread(filename_5d_points, "%s");

% Find Unique entries
[Y, I, J] = unique (REFNR, "first");

%Initialize strings etc.
graphics_toolkit ("gnuplot"); %fixes crash with some Intel drivers

output_string_csv = "REFNR,GPSNR,MIDDLE_EPOCH,XKOOR,YKOOR,ZKOOR,H,Number of measurements, V_H (mm/year)\n";

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
  %if is5Dpoint
  if true
    
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

  %Cartesian
  [px,Sx] = linreg(epochs, x);
  [py,Sy] = linreg(epochs, y);
  [pz,Sz] = linreg(epochs, z);
  
  %UTM
  [pe,Se] = linreg(epochs, eastings);
  [pn,Sn] = linreg(epochs, northings);
  [ph,Sh] = linreg(epochs, heights);
   
  m = length(epochs);
  
  %Vertical velocity mm/anno (from m/day);
  v_h = ph(2)*1000*365.25;
  fprintf('%f.5 mm/year\n',v_h);
  height = polyval([ph(2),ph(1)],epoch_eval_num-epochs0)+heights0;
  fprintf('%f.5 height\n',height);
  
  if figures > 0
    figure;  
    plot(epochs+epochs0,heights+heights0);
    hold on;
    plot(epochs+epochs0, polyval([ph(2),ph(1)],epochs)+heights0);
    plot(epoch_eval_num,height,'*');
    titlestring = strjoin([ gpsnr ' (' num2str(refnr) ') ']);
    title(titlestring);
  end
  
end
end