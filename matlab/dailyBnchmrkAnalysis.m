% Daily Benchmark Analysis
clear all
close all



% Provide directory containing images
ImgDir = 'D:\PfsTests\04_10_14_08_32_05\Images/'; %Must end with /
dateStamp = char(regexp(ImgDir,'\d\d_\d\d_\d\d(?=_\d\d_\d\d_\d\d)','match'));
configFile = 'EMconfig_03-31-14.mat';


%% Theta 1 center metrology
% Get Centroids from images
savedir = 'C:\Users\sage\Desktop\Dropbox\PFS_EM\TEST_RESULTS\Metrology/';
dataName = strcat('thetaCenters_',dateStamp);
ImgFilter = 'centerms1_*.fits'; % make sure to specify fits or ppm 
assgn2one = false;
getCentroids;

% thisfile = fullfile(savedir,strcat(dataName,'.mat'));
% thetaData = S1CentersAnlys(thisfile,dataName,configFile);
% save(thisfile,'thetaData','-append')

%% Phi metrology
% Get Centroids from images
savedir = 'C:\Users\sage\Desktop\Dropbox\PFS_EM\TEST_RESULTS\Metrology/';
dataName = strcat('phiMtrlgy_',dateStamp);
ImgFilter = 'centerms2*.fits'; % make sure to specify fits or ppm 
assgn2one = false;
getCentroids;

% thisfile = fullfile(savedir,strcat(dataName,'.mat'));
% phiData = S2CentersAnlys(thisfile,dataName,configFile);
% save(thisfile,'phiData','-append')