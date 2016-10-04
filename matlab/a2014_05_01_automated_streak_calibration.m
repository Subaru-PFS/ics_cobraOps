%% Center Metrology Script

%addpath('C:/Users/cmorantz/Desktop/PFS/MATLAB/')
addpath('C:/Users/johannes/Dropbox/PFSMatlab/download/')
%addpath('C:\Users\sage\Desktop\Dropbox\PFSMatlab/download/')

close all
clear all
clear data;

stage = 'stage1';
%% Define test directory with folders containing the center metrology images
dotDir = 'D:\PfsTests\04_30_14_07_49_20_Benchmark\Images/';
streakDir = 'D:\PfsTests\02_18_14_11_07_02_PreciseAngleTest10\Images\1fs/';
numStepsLc = 160;

ImgDir = dotDir;
savedir = dotDir;

getCentroids;



numSteps = numStepsLc;
%  imagefile = dir2cell([dotDir '*.ppm']);
centerDirs = dir2cell(dotDir);

prevws = dir2cell([streakDir horzcat('*',stage,'.mat')]);
prevdt = dir2cell([dotDir horzcat('*',stage,'.mat')]);

if numel(prevws) ==0
    images = dir2cell([streakDir '*.ppm']);
      nDir=1; 
    %% Get the StreakData. 
    intensities = [];
    overloads = [];
    imageMatrix = [];
    
    % Angular size of the pie:
    res = 100;
    resolution_angle=2*pi/res;
    resangle = 360/res;
    angles=[resolution_angle:resolution_angle:2*pi];
    
    %% Read Streak Images
    for ll = 1: length(images)
        imagefile = dir2cell([streakDir '*.ppm']);
        imageMatrix=imread(strcat(streakDir, imagefile{ll}));
        imageMatrix=uint32(imageMatrix);
        disp(imagefile);
        [intensity, overload ] = getIntensity2( cobracenter, min_range, max_range, imageMatrix, resolution_angle);
        intensities = vertcat(intensities, intensity);
        overloads = vertcat(overloads, overload);
        imageMatrix = [];
    end
    clear numStepsLc;
    clear prevdt;
    clear dotDir;
   save(horzcat(streakDir,'streakspace',stage,'.mat'));
else
    load(strcat(streakDir, prevws{1}));
    numSteps = numStepsLc;
end
%% GET THE DOTS


if numel(prevdt) ==0
    images = dir2cell([streakDir '*.ppm']);
    data.imageDir = dotDir;
    
    %'C:\Users\johannes\Dropbox\PFSMatlab\Results\2014_01_08_StreakMappingAngles\comparisonBetweenTradDotsAndStreak\09msFWD\dotsImg/';
    data.pId = [3,6,7,58,62];
    % Centers in CCD frame
    data.centers=[631, 623;      % GEN1 @ 1.3
        187,1535;           % FIDUCIAL 1.6
        459,1776;           % FIDUCIAL 1.7
        1414,602;           % GEN2 @ 2.1
        1437,1554           % GEN3 @ 2.5
        ] * [1;i] ;
    data.activeCobras = [58,62];
    data.Fiducials = [1,2,3]; % Same purpose as Fiducials, but more intuitive to use pId
    data.movingCobras = [4,5]; % Same purpose as movingCobras but more intuitive to use pId
    armLengths = [158.3834, 144.0562; % GEN2
        159.4869, 154.8084]; %GEN3
    
    % Get centroids and assign to correct Cobra or fiducial
    data = Get_multiple_centroids(data);
    
    % Apply Horn's method Common Mode subtraction
    %data = cobraCMS(data,true,false,true);
    
    %% ORDER THE DOTS
    data = getS1ArmAngles(data,CobraConfig.(pid).s1Center );
    clear numStepsLc;
    clear prevdt; 
    save(horzcat(dotDir,'wdotspace',stage,'.mat'));
else
    load(strcat(dotDir, prevdt{1}));
    numSteps = numStepsLc;
end

%% Align Coordinate Systems
meansity = [mean(double(imag(intensities)), 1) * 1i + real(intensity)];
extensity = horzcat(double(meansity), double(meansity));
spdds = [];
% bs = zeros(data.Nimg,1000); 
dotSpeeds = cell(data.Nimg);
for jj = 1:data.Nimg
    cfitID = sprintf('cfit0_%d',jj) ;
    centerID=sprintf('img%d',jj) ;
    angles = [];
    angles = data.(pid).(cfitID).thetas + pi;
    angles(end +1) = angles(1) + 2* pi;

    %% Measure Angular Distance and Integrate Streak Portion
    csteps = [];
    imagefile = dir2cell([dotDir '*.ppm']);
    I=imread(strcat(dotDir, imagefile{1}));
    spdd =[];
    temp =[];
    prevAngle = 0;
    
    speedsDot1 = angles(1:end-1)*180/pi + 1i * numSteps*pi./(diff(angles)*180);
    speedsDots = [angles(1:end-1)*180/pi, angles(2:end)*180/pi] + 1i * numSteps*pi./([diff(angles),diff(angles)]*180);
    b = zeros(1,length(angles)*2+5);
    b(1:2:length(angles)*2) = angles*180/pi+ 1i * numSteps*pi./([diff(angles),angles(1)+2*pi-angles(end)]*180);
    b(2:2:length(angles)*2) = [angles(2:end),angles(1)+2*pi]*180/pi + 1i* numSteps*pi./([diff(angles),angles(1)+2*pi-angles(end)]*180);
    del = find(imag(b) >= 20);
    b(del) = [];
    del = find(imag(b) == 0);
    b(del) = [];
    dotSpeeds{jj} = b;
%     b = horzcat(b, zeros(1,1000-length(b)));
%     bs(jj,1:length(b)) = b;
    %figure(5);
    %imshow(I);
    %imagesc(I);
end
% max = floor(7000/numSteps);
% bs = bs(:,1:max);
% bs = horzcat(bs, bs+360);

speeds = [imag(intensities)/(flowRate*3.6) *1i + real(intensities), imag(intensities)/(flowRate*3.6) *1i + real(intensities+360)];
mspeeds = [imag(meansity)/(flowRate*3.6) *1i + real(meansity),imag(meansity)/(flowRate*3.6) *1i + real(meansity+360)];
I = [];
save(horzcat(dotDir,'workspace',stage,'final.mat'));

