%SCRIPT_5D_NET_cartesian_middle_epoch.m
%
% Script to make linear regression of "5D-net" data and output file compatible
% with GNSS trans in cartesian coordinates.
% Heights are in utm32Eetrs89
%
% Regression is carried out individually for each of the cartesian coordinates
%
% oldjo@sdfe.dk, August-October 2017
%
% Requires Octave Forge packages "io" and "statistics"
% >> pkg install -forge io
% >> pkg install -forge statistics
% >> pkg install -forge struct
% >> pkg install -forge optim


graphics_toolkit ("gnuplot"); %fixes crash with some Intel drivers

% Filename
filename_input = 'inputs\\5D_tidsserier_20170504_trimmed.csv';
% Note: this must be a comma seperated file with the following columns
% REFNR,GPSNR,XKOOR,YKOOR,ZKOOR,XRMS,YRMS,ZRMS,JNR_BSIDE,SYS,EPOCH,IN_DATE
filename_5d_points = 'inputs\\5d_points_sorted.csv';
uplifts_filename = 'inputs\\altgps_uplift';

% output filenames
filename_output = 'outputs\\out.txt';
filename_output_csv = 'outputs\\out.csv';
filename_output_koords = 'outputs\\coordinates.txt';
fig_folder = 'figures\\';


% Epoch to evaluate. Format yyyymmdd
% Not used. The epoch will be in the middle of the interval of measurements
%epoch_eval = '20070101'; 

% Alpha
alpha = 0.05;

% Binsize in days
binsize = 14;

% Generate figures (0: no, 1: display, 2: save and display).
figures = 2;
closeFigures = 0; %1;
Nmin = 3;

% Figure format: (png, pdf, jpg... etc.)
figure_format = 'pdf'; 

% Load 'optim' package for 'polyconf'
pkg load optim

% Load 'statistics' package for 'regress'
pkg load statistics

% Formatstring (change this if the input files get different fields.
T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
[REFNR, GPSNR, XKOOR, YKOOR, ZKOOR, SYS, EPOCH] =...
textread(filename_input, T, "delimiter", ",", "headerlines", 1);
[GPSNR_5D] = textread(filename_5d_points, "%s");

% Find Unique entries
[Y, I, J] = unique (REFNR, "first");

% Read the uplifts
uplifts = readUplifts(uplifts_filename,GPSNR(I));

% Convert into meters per day
uplifts = uplifts./(1000*365.25);

%Init strings etc.
residuals = [];
koords_string = "";
output_string = "";
formatSpec = '# crt_igs08 %s \n';

output_string_csv = "REFNR,GPSNR,MIDDLE_EPOCH,XKOOR,YKOOR,ZKOOR,H,...
Number of measurements, V_H (mm/year), lower, upper, V_H (model), t_score \n";

%output_string = sprintf(formatSpec,epoch_eval);
%epoch_eval_num = datenum(datevec(epoch_eval,'yyyymmdd'))

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
  [xkoords, ~] = binning(xkoords,epochs,binsize);
  [ykoords, ~] = binning(ykoords,epochs,binsize);
  [zkoords, epochs] = binning(zkoords,epochs,binsize);
  
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
    
  epochs = epochs - epochs0;
  x = xkoords - xkoords0;
  y = ykoords - ykoords0;
  z = zkoords - zkoords0;
  eastings = eastings - eastings0;
  northings = northings - northings0;
  heights = heights - heights0;
  
  
  %Set epoch eval to middle epoch
  epoch_eval = datestr(epochs0,'yyyymmdd');
  epoch_eval_num = datenum(datevec(epoch_eval,'yyyymmdd'));
  
  %% Fit to polynomial by least squares in each dimenstion
  [px,Sx] = polyfit(epochs, x, 1);
  [py,Sy] = polyfit(epochs, y, 1);
  [pz,Sz] = polyfit(epochs, z, 1);
  %[pE,SE] = polyfit(epochs, eastings, 1);
  %[pN,SN] = polyfit(epochs, northings, 1);
  %[ph,Sh] = polyfit(epochs, heights, 1);
  
  
  % Fit a third time using custom 'linreg' function
  [bh_custom, statsh_custom] = linreg(epochs, heights, alpha, uplifts(i));
  
  % Only save residuals if at least "Nmin" measurements
  if length(epochs) >= Nmin;
    residuals = [residuals; statsh_custom.resid];
  end
  
  
  %Vertical velocity mm/anno (from m/day);
  v_h = bh_custom(2)*1000*365.25;
  
  %% Figures
  if ((figures > 0) && (length(epochs) >= Nmin))
    fig = figure;%('PaperPositionMode','auto');
    set(gcf,'PaperType','A4', ...
         'paperOrientation', 'landscape','PaperPositionMode','auto', ...
         'paperunits','CENTIMETERS', ...
         'PaperPosition',[.63, .63, 28.41, 19.72]);
    t = linspace(epochs(1),epochs(end),25);
    %subplot(1,3,1)
    %[xf,dxf] = polyconf(px,t,Sx,'pi');
    %plot(t,xf,'g-;Fit;',t,xf+dxf,'g-.;95% Confidence Interval;',...
    %  t,xf-dxf,'g.;;',epochs,x,'xr;Data;');
    plot(epochs,heights0+heights,'xr;Data;');
    hold on;
    plot(epoch_eval_num-epochs0,heights0+polyval(bh_custom([2,1])', epoch_eval_num-epochs0),'*c')
    plot(t,heights0+polyval(bh_custom([2,1])', t),'-k')
    legend('off');
    titlestring = strjoin([ gpsnr ' (' num2str(refnr) ') ']);
    
    titlespec = "%s, Slope: %0.5f mm/year";
    title(sprintf(titlespec, titlestring, v_h));
    ylabel('h-anomally [m]')
    xlabel('Relative Epoch [Years]')
    datetick('x','yyyy')
    
    str1 = ["Slope of regression line: ",num2str(v_h)," mm/year"];
    str2 = ["Std.dev.: ",num2str(statsh_custom.sigma_b)," mm/year"];
    str3 = ["95% CI: [",num2str(statsh_custom.bmin)," mm/year; ",...
            num2str(statsh_custom.bmax)," mm/year]"];
    str4 = ["Std.dev.: ",num2str(1000.*statsh_custom.s_n)," mm"];
    
    text(epochs(1),min(heights)+heights0+0.001,str1);
    text(epochs(1),min(heights)+heights0+0.00075,str2);
    text(epochs(1),min(heights)+heights0+0.0005,str3);
    

    if figures == 2
      filename = sprintf('%i_%s.%s',refnr,epoch_eval,figure_format);
      saveas(fig,[fig_folder filename],figure_format)
      if closeFigures
        close(fig);
      end
    end
  end
  
    % Don't record from stations with only one measurement
  if num_measurements > 1
  
  x_eval = polyval(px, epochs0) + xkoords0;
  y_eval = polyval(py, epochs0) + ykoords0;
  z_eval = polyval(pz, epochs0) + zkoords0;

  output_string = [output_string sprintf(formatSpec,epoch_eval)];
  
  %# crt_igs08  19890101
  %1  3513649.63200 m   778954.53800 m   5248201.78400 m
  %-1z
  formatSpec2 = '%i %0.5f m %0.5f m %0.5f m \n';
  
  output_string = [output_string sprintf(formatSpec2, refnr ,x_eval,y_eval,...
  z_eval)];
  
  %Cartesian to geo
  [easting_eval, northing_eval, h_eval] = cartesian_to_UTM32Eetrs89(x_eval,...
    y_eval,z_eval);
  formatSpec_csv = '%i, %s, %s, %.5f, %.5f, %.5f, %.5f, %i, %.5f, %.5f, %.5f, %.5f, %.5f\n';
  output_string_csv = [output_string_csv sprintf(formatSpec_csv, refnr,...
    char(gpsnr),epoch_eval,x_eval,y_eval,z_eval,h_eval,num_measurements,...
    v_h, statsh_custom.bmin.*(1000*365.25),...
    statsh_custom.bmax.*(1000*365.25),...
    uplifts(i).*(1000*365.25), statsh_custom.t_score)];
    
  formatSpec_koords = '%s %.5f %.5f %.5f\n';  
  koords_string = [koords_string sprintf(formatSpec_koords, char(gpsnr),...
    x_eval,y_eval,z_eval)];
    
  end
  end
end

% Save residuals

save("outputs\\residuals.mat","residuals","Nmin");

fprintf('Calculations done.\n');

output_string = [output_string '-1z'];
fileID = fopen(filename_output,'w');
fprintf(fileID,output_string);
fclose(fileID);

fprintf('Text File written.\n')

fileID = fopen(filename_output_csv,'w');
fprintf(fileID,output_string_csv);
fclose(fileID);

fprintf('.CSV File written.\n')

fileID = fopen(filename_output_koords,'w');
fprintf(fileID,koords_string);
fclose(fileID);