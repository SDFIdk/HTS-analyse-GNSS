function [out] = gpsweek(date)
 % Calculates GPSweek from date
  refdate = '20170101';
  refgpsweek = 1930;
  refdatenum = datenum(datevec(refdate,'yyyymmdd'));
  weekdiff = floor((date - refdatenum)./7);
  out = refgpsweek + weekdiff;
 
end