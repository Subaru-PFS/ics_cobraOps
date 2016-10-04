%% Convergence Testing Comparison Script
% This script takes the data directories you provide it which have mat
% files from previously run targConvAnlys script and compares the
% cummulative convergences of the different tests which those directories
% correspond to.

% NOTE: THIS MUST BE RUN FROM ITS DIRECTORY

clear all
close all

%% INPUTS
% Which folders in Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence are going
% to be compared. These must all contain the results of analyzeTargetRun
% script

% Specify which group is being analyzed
Group = 'B';

% Group Definitions
switch Group
    case 'A'
        dataDirs = {
            '06_11_14_17_27_59_TargetRun'
            '06_12_14_17_42_30_TargetRun'
            '06_16_14_10_36_21_TargetRun'
            '06_16_14_23_18_51_TargetRun'
            '06_17_14_11_21_32_TargetRun'
            '06_17_14_22_41_35_TargetRun'
            '06_18_14_21_43_16_TargetRun'
            '06_19_14_13_44_30_TargetRun'
            '06_23_14_12_19_58_TargetRun'
            '06_23_14_17_00_09_TargetRun'
            '06_24_14_08_51_13_TargetRun'
            '06_24_14_17_56_21_TargetRun'};
        cobraMap = [121	123 152 143 149 104 142 103 112]; %SET A
        figbase = 1000;
        
    case 'B'
        dataDirs = {
            '06_26_14_17_18_33_TargetRun'
            '06_27_14_09_41_28_TargetRun'
            '07_01_14_10_12_08_TargetRun'
            '07_01_14_17_09_05_TargetRun'};
        cobraMap = [117	109	145	150	101	102	118	116	106]; %SET B
        figbase = 2000;
        
    case 'C'
        dataDirs = {
            '07_09_14_16_37_46_TargetRun'
            '07_10_14_16_51_12_TargetRun'
            '07_11_14_19_29_30_TargetRun'
            '07_14_14_08_25_24_TargetRun'
            '07_14_14_21_47_40_TargetRun'
            '07_15_14_09_06_05_TargetRun'
            '07_17_14_04_35_18_TargetRun'
            '07_17_14_10_47_26_TargetRun'
            '07_17_14_22_06_23_TargetRun'};
        cobraMap = [105 141 151 146 108 144 122 147 110]; %SET C
        figbase = 3000;
        
    case 'C2'
        dataDirs = {
            '08_25_14_08_24_46_TargetRun'
            '08_25_14_15_04_42_TargetRun'
            '08_26_14_08_25_26_TargetRun'
            '08_26_14_13_32_16_TargetRun'
            '08_27_14_08_24_33_TargetRun'
            '08_27_14_14_17_49_TargetRun'
            '08_27_14_22_49_46_TargetRun'
            '08_28_14_11_25_26_TargetRun'
            '08_28_14_16_40_02_TargetRun'
            '08_28_14_22_25_26_TargetRun'
            '08_29_14_17_03_31_TargetRun'
            '09_02_14_13_23_46_TargetRun'
            '09_03_14_08_38_30_TargetRun'
            '09_03_14_15_18_47_TargetRun'};
        cobraMap = [105 141 151 146 108 144 122 147 110]; %SET C2
        figbase = 5000;
        
    case 'D'
        dataDirs = {
            '07_19_14_11_58_58_TargetRun'
            '07_20_14_10_16_05_TargetRun'
            '07_20_14_18_09_04_TargetRun'
            '07_21_14_07_48_04_TargetRun'
            '07_21_14_23_31_12_TargetRun'
            '07_22_14_14_41_38_TargetRun'
            '07_23_14_15_28_24_TargetRun'
            '07_24_14_14_52_33_TargetRun'
            '07_24_14_22_34_42_TargetRun'
            '07_25_14_20_13_03_TargetRun'
            '07_26_14_09_46_51_TargetRun'
            '07_26_14_17_29_32_TargetRun'
            '07_27_14_16_03_29_TargetRun'
            '07_27_14_22_04_53_TargetRun'
            '07_28_14_08_16_24_TargetRun'
            '07_28_14_13_44_03_TargetRun'
            '07_28_14_21_57_52_TargetRun'
            '07_29_14_08_33_41_TargetRun'
            '07_29_14_13_32_51_TargetRun'
            '07_29_14_21_30_55_TargetRun'
            '07_30_14_08_44_40_TargetRun'
            '07_30_14_14_34_51_TargetRun'
            '07_30_14_22_33_06_TargetRun'};
        cobraMap = [121	123	145	143	149	108	142	116	112]; %SET D
        figbase = 4000;
        
end
 
 
% Specify prefix of files to screen for.
loadFilter = '5um';

%% EXECUTION (DANGER IF YOU MESS WITH THIS)

% Initialize vars and other parms
nTests = length(dataDirs);
 

stickySpots = [];
locOfSpots = [];

for ii=1:nTests
     testid = sprintf('test_%d',ii) ;
     stickySpots.(testid) = [];

    dataDir = dataDirs{ii};
    QQ = loadmats([loadFilter '_mId*.mat'],dataDir);
    eval(unpackstruct(QQ)); 
    fields = fieldnames(QQ);
    %Go over the positioners
    for i = 1:numel(fields)
        if mod(i,2) == 0
          targetArray =  QQ.(fields{i});
          posid = sprintf('pos_%d_J1',i/2) ;
          stickySpots.(testid).(posid) = [];
          % Go over the targets
          for t = 1:length(targetArray)
              % What is a Sticky Spot? : 
              % The movement per step is extremly small (< 0.02 deg/step)
              % Go over the moves
              thisTargetsMap = [];
              for move = 2:length(targetArray(1).iter) 
                  %Filter out all small moves (smaller 10 steps)
                  if abs(targetArray(t).J1_S(move)) > 5
                      %Motor Map Value J1
                      movedAngle = (targetArray(t).J1(move)- targetArray(t).J1(move-1))*180/pi;
                      commandedSteps = targetArray(t).J1_S(move);
                      mapValue = movedAngle/commandedSteps;
                      thisTargetsMap(1,move) = mapValue*sign(commandedSteps);
                      thisTargetsMap(2,move) = movedAngle; 
                  end
                  %Motor Map Value J2
              end
              %is there a consecutive row of low map values?
              isStuck = 0; % 0 not, 1 fwd, 2 rvs
              stuckCount = 0;
              for kk = 1: size(thisTargetsMap,2)
                  if abs(thisTargetsMap(2,kk)) < 5 %only smaller 5 deg moved value.
                  if  0.01 > thisTargetsMap(1,kk) >  0;
                      if isStuck == 1 % Does not exclude stuckness at two distinct points.
                        stuckCount = stuckCount +1; 
                      end
                       isStuck = 1;
                  elseif 0 > thisTargetsMap(1,kk) >  -0.01;
                        if isStuck == 2 % Does not exclude stuckness at two distinct points.
                           stuckCount = stuckCount -1 ; 
                        end
                      isStuck = 2;
                  else
                      isStuck = 0;
                  end
                  end
              end
              if not(ismember(1, targetArray(t).status)) % Only really sticky if target hasnt converged.
              if stuckCount < -3 || stuckCount > 3
                  stickySpots.(testid).(posid)(1,end+1) = stuckCount;
                  stickySpots.(testid).(posid)(2,end) = targetArray(t).J1(9)*180/pi;
                  stickySpots.(testid).(posid)(3,end) = t;
                %  stickySpots.(posid)(3,t) = thisTargetsMap;
              %   locOfSpots(i/2
              end
              end
              
          end
        end
    end
    % Go over the positioners
for i = 1:numel(fields)
        if mod(i,2) == 0
          targetArray =  QQ.(fields{i});
          posid = sprintf('pos_%d_J2',i/2) ;
          stickySpots.(testid).(posid) = [];
          % Go over the targets
          for t = 1:length(targetArray)
              % What is a Sticky Spot? : 
              % The movement per step is extremly small (< 0.02 deg/step)
              % Go over the moves
              thisTargetsMap = [];
              for move = 2:length(targetArray(1).iter) 
                  %Filter out all small moves (smaller 10 steps)
                  if abs(targetArray(t).J2_S(move)) > 3
                      %Motor Map Value J2
                      movedAngle = (targetArray(t).J2(move)- targetArray(t).J2(move-1))*180/pi;
                      commandedSteps = targetArray(t).J2_S(move);
                      mapValue = movedAngle/commandedSteps;
                      thisTargetsMap(1,move) = mapValue*sign(commandedSteps);
                      thisTargetsMap(2,move) = movedAngle;
                      
                    
                  end
                  %Motor Map Value J2
              end
              %is there a consecutive row of low map values?
              isStuck = 0; % 0 not, 1 fwd, 2 rvs
              stuckCount = 0;
              for kk = 1: size(thisTargetsMap,2)
                  if abs(thisTargetsMap(2,kk)) < 5 %only smaller 5 deg moved value.
                  if  0.01 > thisTargetsMap(1,kk) >  0;
                      if isStuck == 1 % Does not exclude stuckness at two distinct points.
                        stuckCount = stuckCount +1; 
                      end
                       isStuck = 1;
                  elseif 0 > thisTargetsMap(1,kk) >  -0.01;
                        if isStuck == 2 % Does not exclude stuckness at two distinct points.
                           stuckCount = stuckCount -1 ; 
                        end
                      isStuck = 2;
                  else
                      isStuck = 0;
                  end
                  end
              end
              if not(ismember(1, targetArray(t).status)) % Only really sticky if target hasnt converged.
              if stuckCount < -3 || stuckCount > 3
                  stickySpots.(testid).(posid)(1,end+1) = stuckCount;
                  stickySpots.(testid).(posid)(2,end) = targetArray(t).J2(9)*180/pi;
                  stickySpots.(testid).(posid)(3,end) = t;
                %  stickySpots.(posid)(3,t) = thisTargetsMap;
              end
              end
              
          end
        end
end
%     for jj = 1:9
%         
%       if isempty(regexp(newname,'_str'))
%             dnames{dc} = newname;
%             dc = dc+1;
%       end
%     end
end
 
