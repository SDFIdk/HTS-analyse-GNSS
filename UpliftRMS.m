%UpliftRMS
%Supporting function that can be used to test data with. Not used in TimeSeriesAnalysis.

% Calculates RMS between fit and uplift model as function of number of 
% datapoints and displays the results.
% Uses data from out_H.csv

%EXAMPLE:
%
%>>UpliftRMS.m
%RMS(m>=3) = 0.992786 mm/year, n = 88
%RMS(m>=4) = 0.794984 mm/year, n = 71
%RMS(m>=5) = 0.773245 mm/year, n = 43
%RMS(m>=6) = 0.394553 mm/year, n = 21
%RMS(m>=7) = 0.331300 mm/year, n = 9
%RMS(m>=8) = 0.370921 mm/year, n = 6
%RMS(m>=9) = 0.208400 mm/year, n = 1

filename_input_csv = 'outputs\\out_H.csv';
% Formatstring (change this if the input files get different fields.

T = "%*s %f %f %*f %*s %*i %*f %*f %*f %*f %*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f"; 
[GPSNR N alfa ref] =...
textread(filename_input_csv, T, "delimiter", ",", "headerlines", 1);

% Prune for NaNs
idx = (~isnan(alfa) & ~isnan(ref));
GPSNR = GPSNR(idx); %GPSNR == REFNR
N = N(idx); %Number of points
alfa = alfa(idx); %dh/dt [mm/year]
ref = ref(idx); %Uplift values from DTU Space [mm/year]

m = unique(N);
RMS = NaN(length(m),1);
n = NaN(length(m),1);

%Makes an RMS value for each series that contains N amounts of points.
%So one RMS calculation for all stations that has 5 points, one for all that has 6 points, etc.
for i = 1:length(m)
  n(i) = length(alfa(N >= m(i)));
  RMS(i) = sqrt((1./n(i))*sum((alfa(N >= m(i)) - ref(N >= m(i))).^2)); %alfa = predicted (regressionsanalyse), ref = known (DTU Space).
  
  if m(i) > 2  
    fprintf('RMS(m>=%i) = %f mm/year, n = %i\n',m(i),RMS(i),n(i));
  end
end

