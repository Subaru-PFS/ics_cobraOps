clear all
close all

dirs = {
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_02_14_11_47_52_Targets_NoDamping\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_02_14_17_49_31_CustomDamping\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_03_14_16_15_49_Targets\output_log.txt';
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_04_14_20_23_02_Targets\output_log.txt';
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_05_14_18_36_53_Targets\output_log.txt';
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_07_14_15_41_24_Targets\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_09_14_18_06_55_Targets\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_10_14_14_30_50_Targets\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_13_14_07_42_20_Targets\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_27_14_21_01_19_Targets\output_log.txt';    
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_16_14_15_20_40_Targets\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_20_14_10_44_35_Targets\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_28_14_11_13_33_Targets\output_log.txt';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_29_14_21_58_02_Targets\output_log.txt'
};




c=1;
for ii=1:length(dirs)
    try
        T = getLogTemps(dirs{ii});

        Tdata(c).logfile = dirs{ii};
        Tdata(c).targets = T.targets';
        Tdata(c).sensor = T.sensor';
        Tdata(c).ambient = T.ambient';
        
        c = c+1;
        
    catch err
        disp(err)
    end
end
    







th1 = figure('name','Test Sensor Temps');
CW = lines(length(dirs));
for ii=1:length(Tdata)
    Tdata(ii).desc = regexp(Tdata(ii).logfile,'\d*_\d*_\d*_\d*_\d*_\d*','match');
    ph(ii) = plot(Tdata(ii).targets,Tdata(ii).sensor,'o','color',CW(ii,:));
    hold on
end

lh1 = legend(ph,[Tdata.desc]);
set(lh1,'interpreter','none')
xlabel('Target #')
ylabel('degC')
title('Sensor Temperature')












th2 = figure('name','Test Ambient Temps');
CW = lines(length(dirs));
for ii=1:length(Tdata)
    Tdata(ii).desc = regexp(Tdata(ii).logfile,'\d*_\d*_\d*_\d*_\d*_\d*','match');
    ph(ii) = plot(Tdata(ii).targets,Tdata(ii).ambient,'o','color',CW(ii,:));
    hold on
end

lh1 = legend(ph,[Tdata.desc]);
set(lh1,'interpreter','none')
xlabel('Target #')
ylabel('degC')
title('Ambient Temperature')
























