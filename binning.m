function [xout, tout] = binning(x,t,binsize)
%BINNING  Bins temporarily close values together with normalized weight
%
% Runs sequentially through the data and sums and normalized points together
% that are within 'binsize' apart in t.
% This potentially results in actual binsizes of 2*binsize for dense data.
% 
% Inputs
%   x = data (column vector)
%   t = time (column vector)
%   binsize = 'size' of bin in units of t
%
% Outputs
%   x = data binned
%   x = time binned
%
% 21/09/2017 Óli D. Jóhannsson



num = ones(size(x)); %vector to keep track of number of measurements in bin
j = 1;
k = length(t);
while j < k
  %if the next t value is within 'binsize' days do:
  if t(j+1) <= (t(j) + binsize)
    if (j+2) <= length(t)
      %If there are at least two more elements in 't'
      x = [x(1:j-1); (num(j)*x(j)+x(j+1))/(num(j)+1); x(j+2:length(x))];
      t = [t(1:j-1); (num(j)*t(j)+t(j+1))/(num(j)+1); t(j+2:length(t))];
      %update weight vector
      num = [num(1:j-1); (num(j)+num(j+1)) ;num(j+2:length(num))]; 
    else
      %If j+1 is the last element in 't'
      x = [x(1:j-1); (num(j)*x(j)+x(j+1))/(num(j)+1) ];
      t = [t(1:j-1); (num(j)*t(j)+t(j+1))/(num(j)+1) ];
      %update weight vector
      num = [num(1:j-1); (num(j)+num(j+1)) ];
    end
    %calculate new length
    k = length(t);
    j = j - 1; %decrease counter
  end
  %increment counter
  j = j+1;
end

xout = x;
tout = t;