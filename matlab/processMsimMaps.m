% processMsimMaps
clear all
close all 

%% INPUTS

% Map to Update? Options are: FAST, SLOW, BOTH
mapToUpdate = 'FAST';

% Provide the config file with which the data was taken with
baseCfg = loadCfgXml;

%How many steps was theta moving between centroids
thetaNsteps = 15;
% thetaNsteps = 120;
%How many steps was phi moving between centroids
phiNsteps = 10;
% phiNsteps = 60;

% Theta FWD Move threshold and persistence to trigger cut
thetaFwMvMin = .2; % Deg
thetaFwPrst  = 3; % Iterations

% Theta REV Move threshold and persistence to trigger cut
thetaRvMvMin = .2;
thetaRvPrst = 3;

% Phi FWD Move threshold and persistence to trigger cut
phiFwMvMin = .2;
phiFwMvPrst = 3;

% Phi REV Move threshold and persistence to trigger cut
phiRvMvMin = .2; 
phiRvMvPrst = 3;

%% EXECUTION (Don't edit below here)
% Get theta maps
thetaFw = makeMsimMaps(1,1,'thetaFwMap',thetaNsteps,thetaFwMvMin,thetaFwPrst);
thetaRv = makeMsimMaps(1,2,'thetaRvMap',thetaNsteps,thetaRvMvMin,thetaRvPrst);

% Get phi maps
phiFw = makeMsimMaps(2,1,'PhiFwMap',phiNsteps,phiFwMvMin,phiFwMvPrst);
phiRv = makeMsimMaps(2,2,'PhiRvMap',phiNsteps,phiRvMvMin,phiRvMvPrst);
% LEGACY, might be necessary to uncomment this to look at old data
% phiFw = makeMsimMaps(2,1,'thetaPhiFwMap',phiNsteps,moveMin,movePersist);
% phiRv = makeMsimMaps(2,2,'thetaPhiRvMap',phiNsteps,moveMin,movePersist);

% Start new config
newConfig = baseCfg;

c=1;
for pid = [thetaFw.pId]
    
    % Store the first value of phiFw results as the minimum phi angle
    
    armid = sprintf('ARM_DATA_%d',pid);
    
    % Save the max/min phi angles to the cobra config structure
    newConfig.ARM_DATA.(armid).phiMax = phiFw(find([phiFw.pId]==pid)).rgns(end-1);
    newConfig.ARM_DATA.(armid).phiMin = phiFw(find([phiFw.pId]==pid)).rgns(1);
    
    fht = figure('name',armid);
    fph1 = plotMotMapsFromStruct(baseCfg.ARM_DATA.(armid).FAST_CALIBRATION_TABLE,'b+--');
    
    subplot(2,2,1)
    tph1 = plot([0;thetaFw(c).rgns],[thetaFw(c).dps(1);thetaFw(c).dps],'rx-');
    legend([fph1.ph(1),tph1],'Old Maps','NewMaps')
    
    subplot(2,2,2)
    plot([0;thetaRv(c).rgns],[thetaRv(c).dps(1);thetaRv(c).dps],'rx-');
    
    subplot(2,2,3)
    plot([0;phiFw(c).rgns],[phiFw(c).dps(1);phiFw(c).dps],'rx-');
    
    subplot(2,2,4)
    plot([0;phiRv(c).rgns],[phiRv(c).dps(1);phiRv(c).dps],'rx-');
    
    savefig(fht,['mmap_' armid '.fig']);
    
    oldorient = getARMval(newConfig,pid,'Global_base_ori_z');
    neworient = oldorient + thetaFw(c).hsoffset;
    if neworient > 360
        neworient = neworient - 360;
    end
%     if thetaFw(c).hsoffset < 180
%         adjustment = thetaFw(c).hsoffset;
%     else
%         adjustment = -(360-thetaFw(c).hsoffset);
%     end
%     neworient = oldorient + adjustment;
    

    newConfig = setARMval(newConfig,pid,'Global_base_ori_z',neworient);
    
    % Set CCW Joint limits to zero
    newConfig = setARMval(newConfig,pid,'Joint1_CCW_limit_angle',0);
    newConfig = setARMval(newConfig,pid,'Joint2_CCW_limit_angle',0);
    
    % Set CW Joint limits to 200 and 400deg for phi/theta respectively
    newConfig = setARMval(newConfig,pid,'Joint1_CW_limit_angle',400);
    newConfig = setARMval(newConfig,pid,'Joint2_CW_limit_angle',200);
    
    switch mapToUpdate
        case 'FAST'
            %FAST MOTOR MAPS
            %Joint 1 FW
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text   = NumArr2mMapStr([length(thetaFw(c).rgns);200;thetaFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = NumArr2mMapStr([length(thetaFw(c).dps);200;thetaFw(c).dps],200);
            % Joint 2 FW
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text   = NumArr2mMapStr([length(phiFw(c).rgns);200;phiFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = NumArr2mMapStr([length(phiFw(c).dps);200;phiFw(c).dps],200);

            % Joint 1 RV
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text   = NumArr2mMapStr([length(thetaRv(c).rgns);200;thetaRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = NumArr2mMapStr([length(thetaRv(c).dps);200;thetaRv(c).dps],200);
            % Joint 2 RV
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text   = NumArr2mMapStr([length(phiRv(c).rgns);200;phiRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = NumArr2mMapStr([length(phiRv(c).dps);200;phiRv(c).dps],200);
        case 'SLOW'
            %SLOW MOTOR MAPS
            %Joint 1 FW
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_fwd_regions.Text   = NumArr2mMapStr([length(thetaFw(c).rgns);200;thetaFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = NumArr2mMapStr([length(thetaFw(c).dps);200;thetaFw(c).dps],200);
            % Joint 2 FW
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_fwd_regions.Text   = NumArr2mMapStr([length(phiFw(c).rgns);200;phiFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = NumArr2mMapStr([length(phiFw(c).dps);200;phiFw(c).dps],200);

            % Joint 1 RV
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_rev_regions.Text   = NumArr2mMapStr([length(thetaRv(c).rgns);200;thetaRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = NumArr2mMapStr([length(thetaRv(c).dps);200;thetaRv(c).dps],200);
            % Joint 2 RV
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_rev_regions.Text   = NumArr2mMapStr([length(phiRv(c).rgns);200;phiRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = NumArr2mMapStr([length(phiRv(c).dps);200;phiRv(c).dps],200);
        case 'BOTH'
            %FAST MOTOR MAPS
            %Joint 1 FW
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text   = NumArr2mMapStr([length(thetaFw(c).rgns);200;thetaFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = NumArr2mMapStr([length(thetaFw(c).dps);200;thetaFw(c).dps],200);
            % Joint 2 FW
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text   = NumArr2mMapStr([length(phiFw(c).rgns);200;phiFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = NumArr2mMapStr([length(phiFw(c).dps);200;phiFw(c).dps],200);

            % Joint 1 RV
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text   = NumArr2mMapStr([length(thetaRv(c).rgns);200;thetaRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = NumArr2mMapStr([length(thetaRv(c).dps);200;thetaRv(c).dps],200);
            % Joint 2 RV
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text   = NumArr2mMapStr([length(phiRv(c).rgns);200;phiRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = NumArr2mMapStr([length(phiRv(c).dps);200;phiRv(c).dps],200);
            %SLOW MOTOR MAPS
            %Joint 1 FW
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_fwd_regions.Text   = NumArr2mMapStr([length(thetaFw(c).rgns);200;thetaFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = NumArr2mMapStr([length(thetaFw(c).dps);200;thetaFw(c).dps],200);
            % Joint 2 FW
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_fwd_regions.Text   = NumArr2mMapStr([length(phiFw(c).rgns);200;phiFw(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = NumArr2mMapStr([length(phiFw(c).dps);200;phiFw(c).dps],200);

            % Joint 1 RV
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_rev_regions.Text   = NumArr2mMapStr([length(thetaRv(c).rgns);200;thetaRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = NumArr2mMapStr([length(thetaRv(c).dps);200;thetaRv(c).dps],200);
            % Joint 2 RV
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_rev_regions.Text   = NumArr2mMapStr([length(phiRv(c).rgns);200;phiRv(c).rgns],200);
            newConfig.ARM_DATA.(armid).SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = NumArr2mMapStr([length(phiRv(c).dps);200;phiRv(c).dps],200);
    end

    plotMotMapsFromStruct(newConfig.ARM_DATA.(armid).FAST_CALIBRATION_TABLE,'go');
    
    c=c+1;
     
end

% Take the last arm data structure and store it as filler for dummy arms
ARM_DATA_FILLER = newConfig.ARM_DATA.(armid);

% Create dummy arm data structures for missing cobras
for pid = 1:max([thetaFw.pId])
    if isempty(find([thetaFw.pId]==pid))
        armid = sprintf('ARM_DATA_%d',pid);
        newConfig.ARM_DATA.(armid) = ARM_DATA_FILLER;
    end
end

% Save new XML with adjusted maps
[xmlfile, xmlfilepath] = uiputfile('*.xml','Save new CobraConfig XML file with new motor maps');
cobraCfg2xml(newConfig,fullfile(xmlfilepath,xmlfile));
        