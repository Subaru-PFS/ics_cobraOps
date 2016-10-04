clear al
close all
clc

fprintf('\n\n')
fprintf('\t Choose an older XML that has the appropriate refpos values\n\n')
oldXML = loadCfgXml
fprintf('\t Choose the current XML that needs refpos values\n\n')
currentXML = loadCfgXml

currentXML.ARM_DATA.refpos.Text = oldXML.ARM_DATA.refpos.Text;

origFilename = currentXML.cfgFile;
[filepath,justname,file_ext] = fileparts(origFilename);
newFilename = [justname '_refpos' file_ext];

cobraCfg2xml(currentXML,newFilename);