%UpliftRMS

% Calculates RMS between fit and uplift model as function of number of 
% datapoints and displays the results in the editor.
% Uses data from out_H.csv (produce of TimeSeriesAnalysis.m)

%EXAMPLE:
%
%>>UpliftRMS.m
%
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
GPSNR = GPSNR(idx);
N = N(idx);
alfa = alfa(idx);
ref = ref(idx);

m = unique(N);
RMS = NaN(length(m),1);
n = NaN(length(m),1);

for i = 1:length(m)
  n(i) = length(alfa(N >= m(i)));
  RMS(i) = sqrt((1./n(i))*sum((alfa(N >= m(i)) - ref(N >= m(i))).^2));
  if m(i) > 2  
    fprintf('RMS(m>=%i) = %f mm/year, n = %i\n',m(i),RMS(i),n(i));
  end
end

