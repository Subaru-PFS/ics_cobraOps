fidDirs = {'D:\PfsTests\06_11_14_17_27_59_TargetRun\Log'
'D:\PfsTests\06_16_14_10_36_21_TargetRun\Log'
'D:\PfsTests\06_18_14_21_43_16_TargetRun\Log'
'D:\PfsTests\06_23_14_12_19_58_TargetRun\Log'
'D:\PfsTests\06_23_14_17_00_09_TargetRun\Log'
'D:\PfsTests\06_26_14_17_18_33_TargetRun\Log'
'D:\PfsTests\06_30_14_15_57_58_dotObscurePhi\Log'
'D:\PfsTests\07_02_14_09_07_18_TargetRun\Log'
'D:\PfsTests\07_10_14_16_51_12_TargetRun\Log'
'D:\PfsTests\07_14_14_08_25_24_TargetRun\Log'
'D:\PfsTests\07_17_14_04_35_18_TargetRun\Log'
'D:\PfsTests\07_20_14_10_16_05_TargetRun\Log'
'D:\PfsTests\07_23_14_15_28_24_TargetRun\Log'
'D:\PfsTests\07_26_14_17_29_32_TargetRun\Log'
'D:\PfsTests\07_29_14_13_32_51_TargetRun\Log'
'D:\PfsTests\07_30_14_22_33_06_TargetRun\Log'
'D:\PfsTests\07_31_14_09_57_25_homeRpt\Log'
'D:\PfsTests\08_05_14_13_40_31_thetaTilt\Log'};

moveToDir = 'C:\Users\sage\Desktop\Dropbox\PFS_EM\TEST_RESULTS\FiducialMeasurements';

for ii=1:length(fidDirs)
    folderName = regexp(fidDirs{ii},'(\d+_)+.+(?=\\)','match');
    newname = strcat('FID_',folderName{1},'.txt');
    command = sprintf('copy /Y %s\\Fiducials.txt %s\\%s',fidDirs{ii},moveToDir,newname);
    [status,cmdout] = system(command);
    disp(cmdout)
    
end


