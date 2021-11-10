function [out] = gpsweek_to_date(gpsweek)
  %GPSWEEK_TO_DATE gpsweek -> datenum
  %  turns a gpsweek into 1st day of the week in datenum format
 
  refdate = '20170101';
  refgpsweek = 1930;
  refdatenum = datenum(datevec(refdate,'yyyymmdd'));
  
  daydiff = 7*(gpsweek-refgpsweek)
  
  
  out = refdatenum + daydiff
  epochs = datenum(refdate,'yyyymmdd')
end