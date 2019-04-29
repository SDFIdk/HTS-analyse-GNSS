function uplifts = readUplifts(filename,gpsnames)
  
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
      
    
    %Check if name came up
    if isempty(ind)
      uplifts(i) = NaN; %assign NaN when no value found
    else  
      uplifts(i) = UPLIFT(ind);  
    end
  end
end