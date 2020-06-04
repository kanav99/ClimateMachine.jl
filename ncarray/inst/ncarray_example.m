
% Tutorial for using ncArray

% It is advised to run this script in an empty directory.
% It will delete and overwrite files named file1.nc, file2.nc and file3.nc.

% size of the example data (2x3)

n = 3;
m = 2;

% create 3 files (file1.nc, file2.nc,...) with a 2x3 variable called SST
data = zeros(n,m);

disp('create example files: file1.nc, file2.nc, file3.nc')

for i = 1:3  
    data(:) = i;
    files{i} = sprintf('file%d.nc',i);
    delete(files{i});
    ncarray_example_file(files{i},data);
end


% Using ncArray

SST = ncArray('file1.nc','SST');

disp('load the entire file')
data = SST(:,:,:);

disp('get the attribute units')
units = SST.units;


disp('load a particular value');
data = SST(3,2,1);

% Using ncCatArray

disp('concatenate the files over the 3rd dimension (here time)')

SST = ncCatArray(3,{'file1.nc','file2.nc','file3.nc'},'SST');

% or just 
% SST = ncCatArray(3,'file*.nc','SST');

disp('load all 3 files');
data = SST(:,:,:);

disp('load a particular value (1,2,1) of the 3rd file');
data = SST(1,2,3);
