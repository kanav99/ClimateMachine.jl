% Parse netCDF time unit.
%
% [t0,f] = nctimeunits(u)
%
% Parse the netCDF time unit u and returns the time offset (days since 31 
% December 1 BC, as datenum) and scaling factor f (in days).
% See the netCDF CF convention for the structure of the time units.
% http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#time-coordinate
% Also: http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM/CalendarDateTime.html

function [t0,f] = nctimeunits(u)

% years in days for udunits
% http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#time-coordinate
year_in_days =  365.242198781;


l = strfind(u,'since');

if length(l) ~= 1
  error(['time units sould expect one "since": "' u '"']);
end

period = strtrim(lower(u(1:l-1)));
reference_date = strtrim(u(l+6:end));

if strcmp(period,'millisec') || strcmp(period,'msec')
  f = 1/(24*60*60*1000);
elseif strcmp(period,'second') || strcmp(period,'seconds') ...
      || strcmp(period,'s') || strcmp(period,'sec')
  f = 1/(24*60*60);
elseif strcmp(period,'minute') || strcmp(period,'minutes') ...
      || strcmp(period,'min')
  f = 1/(24*60);
elseif strcmp(period,'hour') || strcmp(period,'hours') ...
      || strcmp(period,'hr')
  f = 1/24;
elseif strcmp(period,'day') || strcmp(period,'days')
  f = 1;
elseif strcmp(period,'week') || strcmp(period,'weeks')
  f = 1/(24*60*60*7);
elseif strcmp(period,'year') || strcmp(period,'years') ...
      strcmp(period,'yr')
  f = year_in_days;
elseif strcmp(period,'month') || strcmp(period,'months') ...
      strcmp(period,'mon')
  f = year_in_days/12;
else
  error(['unknown units "' period '"']);
end


if strcmp(reference_date,'1900-01-01 00:00:0.0')
  t0 = datenum(1900,1,1);
else
  try
    t0 = datenum(reference_date,'yyyy-mm-dd HH:MM:SS');
  catch
    try
      t0 = datenum(reference_date,'yyyy-mm-ddTHH:MM:SS');
    catch    
      try
        t0 = datenum(reference_date,'yyyy-mm-dd');
      catch
        error(['date format is not recogized ' reference_date])
      end
    end
  end
end
