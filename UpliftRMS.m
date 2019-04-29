%UpliftRMS


filename_input_csv = 'outputs\\out.csv';
% Formatstring (change this if the input files get different fields.
%T = "%d %s %f %f %f %*f %*f %*f %*d %s %s %*s";
T = "%s %*d %d %*s %*d %*f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*d %*d";
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

