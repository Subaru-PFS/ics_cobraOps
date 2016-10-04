function cdf=readallcdf(varargin)
% READALLCDF: read all cdf files in current directory
%
% optional arg: 'ls2cell', [argument to ls], ie '*.cdf'
%               'var', [data structure], ie cdf
%               'plot' (declare this to run plotcdf)

ls_arg = '*.cdf';
ndata = 0;
makeplots = 0;
for jj = 1:2:length(varargin)
  if ~isempty(regexp('var',['^' varargin{jj}]))
    cdf = varargin{jj+1};
  elseif ~isempty(regexp('ls2cell',['^' varargin{jj}]))
    ls_arg = varargin{jj+1};
  elseif ~isempty(regexp('plot',['^' varargin{jj}]))
    makeplots = 1;
    if (jj < length(varargin))
      jj = jj - 1;
    end
  else
    disp('unrecognized argument to readallcdf');
    help readallcdf;
    return
  end
end

files = ls2cell(ls_arg);
load_files = logical(ones(size(files)));

if exist('cdf','var')
  % determine which files have already been read in
  dataset = fields(cdf);
  ndata = length(dataset);
  for kk=1:ndata
    for jj=1:length(load_files)
      if ~load_files(jj)
        continue
      end
      if ~isempty(strfind(files{jj},cdf.(dataset{kk}).header.Filename))
        load_files(jj) = 0;
      end
    end
  end
end

OKtoLoad = find(load_files);
if (length(OKtoLoad)==0)
  disp('READALLCDF: nothing to do!');
  return
end

total_datasets = length(OKtoLoad) + ndata;

for kk = 1:length(OKtoLoad) 
  jj = OKtoLoad(kk);
  ndata = ndata + 1;
  if total_datasets<10
    datasetID = sprintf('a%d',ndata);
  elseif total_datasets < 100
    datasetID = sprintf('a%02d',ndata);
  elseif total_datasets < 1000
    datasetID = sprintf('a%03d',ndata);
  end
  cmd = sprintf('cdf.%s = readcdf(''%s'');',datasetID,files{jj});
  disp(cmd);
  eval(cmd);
  if makeplots
    plotcdf(cdf.(datasetID));
    input('next...');
  end
end
