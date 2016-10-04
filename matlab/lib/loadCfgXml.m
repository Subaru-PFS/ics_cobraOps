function CobraConfig = loadCfgXml(cfgFilePath,cfgFile)
% CobraConfig = loadCfgXml(cfgFilePath,cfgFile)
% dialog box opens if path or file are not specified.

% $$$ if nargin==0
% $$$     % Load config xml through dialog box
% $$$     [cfgFile cfgFilePath] = uigetfile('*.xml','Select cobra config xml to load');
% $$$ end

if ~exist('cfgFilePath','var'), cfgFilePath = '.'; end
if ~exist('cfgFile','var')
  warning off ;
  try
    cfgfiles = dir2cell('*.xml');
  catch
    cfgfiles = {};
  end
  warning on ;
  %% if there is one file in the default directory, and the directory was specified, then prompt the
  %% user on the command line to use it.  only an 'n' will trigger the GUI.
  %
  % with the directory unspecified, unique file in the directory is taken.
  if length(cfgfiles) == 1
    if nargin > 0
      cli_answer = input(sprintf('Use %s? [Y|n] ',cfgfiles{1}), 's');
    else
      cli_answer = 'Y';
      fprintf(1,'Using %s...\n',cfgfiles{1});
    end
    if strcmp(cli_answer,'n')
      [cfgFile cfgFilePath] = uigetfile('*.xml','Select cobra config xml to load');
    else
      %% this is the default for 1 file in the specified directory.
      cfgFile = cfgfiles{1};
    end
    %% if there are no files or more than one file, then trigger
    %% the GUI.
  else
    [cfgFile cfgFilePath] = uigetfile('*.xml','Select cobra config xml to load');
  end
end
   
if ~exist('cfgFile','var') | cfgFile == 0
  disp('No file loaded');
  CobraConfig = [];
  return;
end

CobraConfig = xml2struct(fullfile(cfgFilePath,cfgFile));
CobraConfig.cfgFile = cfgFile;
CobraConfig.cfgPath = cfgFilePath;
