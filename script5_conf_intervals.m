%SCRIPT5_conf_intervals.m
% oldjo@sdfe.dk, November 2017

graphics_toolkit ("gnuplot"); %fixes crash with some Intel drivers

alpha = 0.05;
sigma0 = 4.01401;

% Filename
filename_input = 'inputs\\5D_tidsserier_20170504_trimmed.csv';
% Note: this must be a comma seperated file with the following columns
% REFNR,GPSNR,XKOOR,YKOOR,ZKOOR,XRMS,YRMS,ZRMS,JNR_BSIDE,SYS,EPOCH,IN_DATE
filename_5d_points = 'inputs\\5d_points_sorted.csv';
uplifts_filename = 'inputs\\altgps_uplift';

%%Figures
figures = 2; %1 show, 2 show and save
fig_folder = 'conf_intervals\\figures\\'
figure_format = 'pdf';
closeFigures = 1;
% output filenames
filename_output_csv = 'conf_intervals\\out.csv';

%Bining
binning = 0; %1: on, 0: off.
binsize = 14; % Binsize in days

% Formatstring (change this if the input files get different fields.
T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
[REFNR, GPSNR, XKOOR, YKOOR, ZKOOR, SYS, EPOCH] =...
textread(filename_input, T, "delimiter", ",", "headerlines", 1);
[GPSNR_5D] = textread(filename_5d_points, "%s");

% Find Unique entries
[Y, I, J] = unique (REFNR, "first");

%Initialize strings etc.
output_string_csv = "GPSNR, REFNR, N, Beta1, ConfEstLow, ConfEstHigh, ConfLow, ConfHigh\n";
formatspec_csv = "%s, %i, %i,%.5f,%.5f,%.5f,%.5f,%.5f\n";

uplifts = readUplifts(uplifts_filename,GPSNR(I));

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
  %if (gpsnr{} == "KVND")
  
  %% Only use data from points designated in 5d_points file.
  if (is5Dpoint && ~isnan(uplifts(i)))
    uplift = uplifts(i);
    
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
  
  %Convert from m/day to mm/year
  epochs = epochs./365.25;
  heights = 1000.*heights;
  
  %Call linreg function
  [b, stats] = linreg(epochs, heights, alpha, 0, sigma0);
  
  %Set title for plots if any.
  titlestring = ['gpsnr: ' gpsnr{} ', n_{measurements}: ' num2str(num_measurements)];
  
  known = stats.confinterval_known;
  esti = stats.confinterval_estimated;
  output_string_csv = [output_string_csv sprintf(formatspec_csv,...
             gpsnr{}, refnr, num_measurements, b(2), esti(1), esti(2), known(1), known(2))];
  if (figures > 0 && length(epochs) > 2)  
    epochs = epochs + epochs0/365.25;
    heights = heights +heights0*1000;
    fig = linregplot(epochs,heights,alpha,uplift,sigma0,titlestring);
    title(titlestring);
    set(gcf, 'Position', get(0, 'Screensize'));
    
    
    if figures == 2
      filename = sprintf('%s.%s',gpsnr{},figure_format);
      saveas(fig,[fig_folder filename],figure_format)
      if closeFigures
        close(fig);
      end
    end
  end
end
end


 
fileID = fopen(filename_output_csv,'w');
fprintf(fileID,output_string_csv);
fclose(fileID);

fprintf('.CSV File written.\n')