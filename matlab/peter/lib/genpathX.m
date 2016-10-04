function pathstr = genpathX(basedir,varargin)
% path = genpathX(dir, varargin)
%
% generate path under dir while excluding directories that match
% expressions given as varargin arguments

path_raw = genpath(basedir);

for jj=1:size(varargin,2)
  path_raw = regexprep(path_raw,[basedir '[^:]*?' varargin{jj} '.*?:'],'');
end

pathstr = path_raw;