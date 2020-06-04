% create an example NetCDF file with the name filename and given data

function ncarray_example_file(filename,data)

  dtype = 'double';
  sz = size(data);

  % Variables
  nccreate(filename,'lon','Format','classic','Datatype',dtype,...
           'Dimensions',{'x',sz(1), 'y',sz(2)});
  ncwriteatt(filename,'lon','long_name','Longitude')
  ncwriteatt(filename,'lon','units','degrees_east')
  
  nccreate(filename,'lat','Datatype',dtype,'Dimensions',{'x',sz(1), 'y',sz(2)});
  ncwriteatt(filename,'lat','long_name','Latitude')
  ncwriteatt(filename,'lat','units','degrees_north')
  
  nccreate(filename,'time','Datatype',dtype,'Dimensions',{'time',1});
  ncwriteatt(filename,'time','long_name','Time')
  ncwriteatt(filename,'time','units','days since 1858-11-17 00:00:00 GMT')
  
  nccreate(filename,'SST','Datatype',dtype,'Dimensions',...
           {'x',sz(1), 'y',sz(2), 'time',1});
  ncwriteatt(filename,'SST','missing_value',single(9999))
  ncwriteatt(filename,'SST','_FillValue',single(9999))
  ncwriteatt(filename,'SST','units','degC')
  ncwriteatt(filename,'SST','long_name','Sea Surface Temperature')
  ncwriteatt(filename,'SST','coordinates','lat lon')

  ncwrite(filename,'SST',data);  

  % Copyright (C) 2012,2013,2015 Alexander Barth <barth.alexander@gmail.com>
  %
  % This program is free software; you can redistribute it and/or modify
  % it under the terms of the GNU General Public License as published by
  % the Free Software Foundation; either version 2 of the License, or
  % (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful,
  % but WITHOUT ANY WARRANTY; without even the implied warranty of
  % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  % GNU General Public License for more details.
  %
  % You should have received a copy of the GNU General Public License
  % along with this program; If not, see <http://www.gnu.org/licenses/>.

