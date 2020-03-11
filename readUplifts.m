function uplifts = readUplifts(filename,gpsnames)
  % DESCRIPTION: %
  % reads in all 146 rows of stations from the DTU Space file, and overwrites these values.
  % For each entry in gpsnames (i.e. those stations we want to compare with DTU Space), 
  % it runs a comparison for all entries from DTU Space and finds the station that fits.
  % If gpsnames = gpsnr (sdfe station == DTU Space station), it moves on. Otherwise, if
  % no match is found (i.e. isempty(ind)), the uplift value of this station is then set to NaN, 
  % as there is nothing to compare with. Then uplifts(i) is set as UPLIFT(ind).
  % The function ends up with a vector (uplifts) who contains NaN values in all the places where there
  % was no corresponding station found from DTU Space, and with a value (from UPLIFTS vector, capital, as calculated in altgps_uplift), 
  % which contains the value from DTU Space.
  % So UPLIFTS is a vector containing all uplift values from DTU Space, but only in the fields where a corresponding
  % value was found from the SDFE stations. 
  % These are then the values that are later compared with in TimeSeriesAnalysis (uplift = uplifts(i) ). This uplift
  % value is then read when calling the linreg function [b, stats] = linreg(epochs, heights, alpha, uplift, sigma0_all_h);.
  % %
  
  %inputs:
  %filename = filename_uplifts; 
  %gpsnames = GPSNR(I);
  
  % Formatstring (change this if the input files get different fields.
  T = "%s %f %f %f";
  [GPSNR, XKOOR, YKOOR, UPLIFT] =...
    textread(filename, T); 
    
  uplifts = NaN(length(gpsnames),1);  
  
  %for each of the supplied gpsnames
  for i = 1:length(gpsnames)
    %clear last result
    ind = []; 
    
    %try to locate uplift value for desired gpsname(nr)
    %ind = find(GPSNR == gpsnames(i));
    ind = [];
    for j = 1:length(GPSNR)
      if GPSNR{j} == gpsnames{i}
        ind = j;
      end
    end
      
    
    %Check if the name came up
    if isempty(ind)
      uplifts(i) = NaN; %assign NaN when no value found
    else  
      uplifts(i) = UPLIFT(ind);  
    end
  end
end
