function CobraConfig = modifyDampingJoeEdit(CobraConfigOrig, dataDir, maxSteps, revUpdate, fwdUpdate)
% Calculates damping coefficent and updates motor map
% CobraConfigOrig : Input testbed configuration (including motor maps)
% dataDir         : load mat files from this directory
% maxSteps        : ignore convergences that take more than this number of steps.
% revUpdate       : 
% fwdUpdate       : 

% [PHM] major update 7/10/2014, code cleanup.  output verified on
% TargetConvergence/07_09 data set.
%keyboard;
%% EXECUTION CODE
% Load all mat files in the dataDir
convData = loadmats('*.mat',dataDir);
eval(unpackstruct(convData));
clear convData;

S = whos;
Names = {S.name};
SInd = find(~cellfun(@isempty,regexp(Names,'^mId_\d*_pId_\d*_str','match')));
jj = 1;
% $$$ fh1 = figure('Name','Error vs Request');
% $$$ fh2 = figure('Name','Error Residuals');
CW = lines(length(SInd));
CobraConfig = CobraConfigOrig;
alpha_beta_fid = fopen('alphabeta.txt','w');
slope_fid = fopen('slopes.txt','w');
fprintf(slope_fid,'pId,Theta_fwd,Theta_rev,Phi_fwd,Phi_rev\n');

for ii=SInd
    data0 = eval(Names{ii});
    pId = data0(1).pid;
% $$$     pId = regexp(Names(ii),'(?<=pId_)\d*','match');
% $$$     pId = pId{1}{1};
% $$$     legendArr{jj} = num2str(pId); % only for err vs req plots
    
    data1 = analyzeTC_Joe(data0, CobraConfig, maxSteps);   
end
fclose(alpha_beta_fid);
fclose(slope_fid);
fprintf(1,'\n');