function [out] = gpsweek_to_date(gpsweek)
  %GPSWEEK_TO_DATE gpsweek -> datenum
  %Turns a gpsweek into 1st day of the week in datenum format
  %Identical to the gpsweek program, except it calculates in days instead of weeks, and turns this number into a datenum format (epoch).
 
  refdate = '20170101';
  refgpsweek = 1930;
  refdatenum = datenum(datevec(refdate,'yyyymmdd'));
  
  daydiff = 7*(gpsweek-refgpsweek)
  
  out = refdatenum + daydiff
  epochs = datenum(refdate,'yyyymmdd')
end