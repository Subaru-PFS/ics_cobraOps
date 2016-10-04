function dircells = dir2cell(option_string)
% DIR2CELL reads in a list of files and puts them in a cell array
%
% USAGE: dir2cell('*.cdf')
%
% APPLICATION EXAMPLE:
%  files = dir2cell(file_descriptor);
%  nfiles = length(files);
%  for jj=1:nfiles
%    ...
%  end
%
% 06/30/2006: written by Peter Mao
% 12/20/2007: modified to be platform independent by Tak Kunihiro

[pathstr, name, ext] = fileparts(option_string);
files = dir(option_string);

for ii_file = 1:length(files)
    filenames{ii_file} = fullfile(pathstr,files(ii_file).name);
    times(ii_file)     = datenum(files(ii_file).date);
end

[times_t index] = sort(times);

for ii_file = 1:length(files)
    filenames_t{ii_file} = filenames{index(ii_file)};
end

dircells = filenames_t;
return;

