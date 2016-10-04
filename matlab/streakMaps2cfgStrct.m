function CobraConfig = streakMaps2cfgStrct(CobraConfigOrig, thetaMMapWorkspace, phiMMapWorkspace, nTheta, nPhi, pIdMapping)
% 
% To use this function for only one stage, pass the other stage workspace
% as 0. Eg: CobraConfig = streakMaps2cfgStrct(OriginalCfgStrct, 0, phiMMapWorkspace)
% 
% The pIdMapping array is to be used if the pid# fields of motormap
% structures do not correspond to correct positioners. Remap them with this
% array. Eg if there is no pId4 in reality: pIdMapping = [1,2,3,5,6]
% instead of default [1,2,3,4,5].

if thetaMMapWorkspace ~= 0
    load(thetaMMapWorkspace)
    thetaMotorMapFw = motorMapFw;
    thetaMotorMapRv = motorMapRv;
end;

if phiMMapWorkspace ~= 0
    load(phiMMapWorkspace)
    phiMotorMapFw = motorMapFw;
    phiMotorMapRv = motorMapRv;
end;

if ~exist('pIdMapping','var')
    pIdMapping = 1:length(fields(motorMapFw));
end

% Initialize new structure
CobraConfigFromStreakMaps = CobraConfigOrig;

% Initialize color wheel and counter
CW = lines(length(pIdMapping));
c = 1;

% For each positioner do the following:
for ii=pIdMapping
    
    % Create field strings to reference positioner structure
    cnfgID = ['ARM_DATA_' num2str(ii)];
    pidStr = ['pId' num2str(c)];
    
    % Call up positioner figure
    existingFigure(cnfgID);
    
    % Write theta maps to structure if exist
    if exist('thetaMotorMapFw','var')

        % Check that array for theta contains correct number of elements
        if length(thetaMotorMapFw.(pidStr)(2,:)) ~= nTheta
            disp('Wrong number of elements found for theta. please debug');
            keyboard
        end
        
        % Write stepsize arrays to structure
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = ...
            NumArr2mMapStr([nTheta,100,thetaMotorMapFw.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = ...
            NumArr2mMapStr([nTheta,100,thetaMotorMapRv.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text;
        
        % Write region arrays to structure
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text = ...
            NumArr2mMapStr([nTheta,100,thetaMotorMapFw.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text = ...
            NumArr2mMapStr([nTheta,100,thetaMotorMapRv.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text;
        
        % Plot the maps
        subplot(2,2,1)
        plot(thetaMotorMapFw.(pidStr)(2,3:end),'b--')
        hold on;
        subplot(2,2,2)
        plot(thetaMotorMapRv.(pidStr)(2,3:end),'b--')
        hold on;
        
    end
        
    if exist('phiMotorMapFw','var')
        
        % Check that array for phi contains correct number of elements
        if length(phiMotorMapFw.(pidStr)(2,:)) ~= nPhi
            disp('Wrong number of elements found for phi. please debug');
            keyboard
        end
        
        if length(phiMotorMapFw.(pidStr)(2,:)) ~= nPhi
            keyboard
        end
        
        
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = ...
            NumArr2mMapStr([nPhi,100,phiMotorMapFw.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = ...
            NumArr2mMapStr([nPhi,100,phiMotorMapRv.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text;
        
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text = ...
            NumArr2mMapStr([nPhi,100,phiMotorMapFw.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text = ...
            NumArr2mMapStr([nPhi,100,phiMotorMapRv.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text;
        
        % Plot the maps
        subplot(2,2,3)
        plot(phiMotorMapFw.(pidStr)(2,3:end),'b--')
        subplot(2,2,4)
        plot(phiMotorMapRv.(pidStr)(2,3:end),'b--')
    end

    % Incriment counter
    c = c+1;
end

CobraConfig = CobraConfigFromStreakMaps;
end