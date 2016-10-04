function lscells = ls2cell(option_string)
% LS2CELL reads in a list of files and puts them in a cell array
%
% USAGE: ls2cell('*.cdf')
%
% APPLICATION EXAMPLE:
%  files = ls2cell(file_descriptor);
%  nfiles = length(files);
%  for jj=1:nfiles
%    ...
%  end

  
try
    files = ls('-1L',option_string);
catch
    files = '';
end
delimits = [0 findstr(char(10), files)];
filename = {};
for jj = 1:length(delimits)-1
  filename{jj} = files( delimits(jj)+1 : delimits(jj+1)-1 );
end
lscells = filename;
return;
