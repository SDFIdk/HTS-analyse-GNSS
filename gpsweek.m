function [out] = gpsweek(date)
  %DESCRIPTION:
  %When called through TimeSeriesAnalysis, it outputs a number in the datenum format, of the week difference 
  %between the input date and the reference date, which is the 1st of january, 2017.
  %This week difference is then added to the reference week, which is 1930.
  %If a station has a negative datenum number, it is because this station was measured before 1st. January 2017.
  %The reason the reference week is set to 1930 is because this is the officially defined reference week internationally.
  
 % Calculates GPSweek from date
  refdate = '20170101';
  refgpsweek = 1930; %Amount of weeks since 1st of January, 1980 until 1st. of January 2017.
  refdatenum = datenum(datevec(refdate,'yyyymmdd'));
  weekdiff = floor((date - refdatenum)./7);
  out = refgpsweek + weekdiff;
 
end