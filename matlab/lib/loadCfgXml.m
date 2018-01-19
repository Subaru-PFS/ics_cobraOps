function CobraConfig = loadCfgXml(cfgFilePath,cfgFile)
% CobraConfig = loadCfgXml(cfgFilePath,cfgFile)
% dialog box opens if path or file are not specified.
%
% If there is only one xml file in the current directory, then use
% that one.  If there are multiple xml files and the cfgFile is not
% specified, then bring up the dialog box UNLESS "usedXMLFile.xml"
% is in the directory, in which case we use it.

if ~exist('cfgFilePath','var')
    cfgFilePath = '.'; 
end
if ~exist('cfgFile','var') 
    if exist(fullfile(cfgFilePath,'usedXMLFile.xml'), 'file')
        cfgFile = 'usedXMLFile.xml';
    else
        warning off ;
   
        try
            cfgfiles = dir2cell([cfgFilePath '/*.xml']);
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
            end
            if strcmp(cli_answer,'n')
                [cfgFile cfgFilePath] = uigetfile('*.xml','Select cobra config xml to load');
            else
                %% this is the default for 1 file in the specified directory.
                fprintf(1,'Using %s...\n',cfgfiles{1});
                [fdir, fname, fext] = fileparts(cfgfiles{1});
                cfgFile = [fname fext];
            end
        else
            %% if there are no files or more than one file, then
            %% trigger the GUI.
            [cfgFile cfgFilePath] = uigetfile('*.xml','Select cobra config xml to load');
        end
    end
end

if ~exist('cfgFile','var') | cfgFile == 0
    disp('No file loaded');
    CobraConfig = [];
    return;
end

if ~exist(fullfile(cfgFilePath,cfgFile))
    [cfgFile cfgFilePath] = uigetfile('*.xml','File not found -- select cobra config xml to load');
end

CobraConfig = xml2struct(fullfile(cfgFilePath,cfgFile));
CobraConfig.cfgFile = cfgFile;
CobraConfig.cfgPath = cfgFilePath;
