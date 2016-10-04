function CobraConfig = modifyDamping(CobraConfigOrig, dataDir, maxSteps, revUpdate, fwdUpdate, slowUpdate, fastUpdate)
% Calculates damping coefficent and updates motor map
% CobraConfigOrig : Input testbed configuration (including motor maps)
% dataDir         : load mat files from this directory
% maxSteps        : ignore convergences that take more than this number of steps.
% revUpdate       : 
% fwdUpdate       : 
% slowUpdate      : update slow map (true or false)
% fastUpdate      : update fast map (true or false)

% [PHM] major update 7/10/2014, code cleanup.  output verified on
% TargetConvergence/07_09 data set.

%% EXECUTION CODE
% Load all mat files in the dataDir
convData = loadmats('*.mat',dataDir);
eval(unpackstruct(convData));
clear convData;

% List all variables in workspace
S = whos;
Names = {S.name};

% Find names who match convergence data format
SInd = find(~cellfun(@isempty,regexp(Names,'^mId_\d*_pId_\d*_str','match')));

% Initialize data set counter
jj = 1;

% Create color wheel for data sets
CW = lines(length(SInd));

% Initialize new cobraConfig to provided one
CobraConfig = CobraConfigOrig;

% Open slope and motor model files for writting to
alpha_beta_fid = fopen('alphabeta.txt','w');
slope_fid = fopen('slopes.txt','w');
fprintf(slope_fid,'pId,Theta_fwd,Theta_rev,Phi_fwd,Phi_rev\n');

% This loop runs for each data set
for ii=SInd
    % Store data structure in data0 var
    data0 = eval(Names{ii});
    
    % Get PID
    pId = data0(1).pid;
    
    % Create string for accessing struture arm data
    armfld = sprintf('ARM_DATA_%d',pId);
    
    % Analyze data0 for error vs request slope and motor model variables
    data1 = analyzeTC_CIT(data0, CobraConfig, maxSteps);

    % 
    J1posSlope = data1.J1.posslope;
    J2posSlope = data1.J2.posslope;
    J1negSlope = data1.J1.negslope;
    J2negSlope = data1.J2.negslope;
    
    updateMaps = (revUpdate || fwdUpdate) && (slowUpdate || fastUpdate);
    
    if updateMaps
        % Make a copy of original config
        CobraConfigOld = CobraConfigOrig;
        
        % Get motor maps from config structure and convert to numeric
        % arrays
        Fast_J1_rev_vmap = getMMap(CobraConfig,pId,'FAST',1,'rev');
        Fast_J1_fwd_vmap = getMMap(CobraConfig,pId,'FAST',1,'fwd');
        Fast_J2_rev_vmap = getMMap(CobraConfig,pId,'FAST',2,'rev');
        Fast_J2_fwd_vmap = getMMap(CobraConfig,pId,'FAST',2,'fwd');
        Slow_J1_rev_vmap = getMMap(CobraConfig,pId,'SLOW',1,'rev');
        Slow_J1_fwd_vmap = getMMap(CobraConfig,pId,'SLOW',1,'fwd');
        Slow_J2_rev_vmap = getMMap(CobraConfig,pId,'SLOW',2,'rev');
        Slow_J2_fwd_vmap = getMMap(CobraConfig,pId,'SLOW',2,'fwd');
        
        % Initialize new motor maps as copies of originals
        newFast_J1_rev_vmap = Fast_J1_rev_vmap;
        newFast_J1_fwd_vmap = Fast_J1_fwd_vmap;
        newFast_J2_rev_vmap = Fast_J2_rev_vmap;
        newFast_J2_fwd_vmap = Fast_J2_fwd_vmap;
        newSlow_J1_rev_vmap = Slow_J1_rev_vmap;
        newSlow_J1_fwd_vmap = Slow_J1_fwd_vmap;
        newSlow_J2_rev_vmap = Slow_J2_rev_vmap;
        newSlow_J2_fwd_vmap = Slow_J2_fwd_vmap;
        
        if fwdUpdate
            if fastUpdate
                % Apply damping correction
                % First 2 entries are size and number of entries in motor map.
                newFast_J1_fwd_vmap(3:end) = Fast_J1_fwd_vmap(3:end) * (1-J1posSlope.fast.a);
                newFast_J2_fwd_vmap(3:end) = Fast_J2_fwd_vmap(3:end) * (1-J2posSlope.fast.a);
                % Save new map to config structure
                CobraConfig = setMMap(CobraConfig,pId,'FAST',1,'fwd',newFast_J1_fwd_vmap);
                CobraConfig = setMMap(CobraConfig,pId,'FAST',2,'fwd',newFast_J2_fwd_vmap);
            end
            if slowUpdate
                % Apply damping correction
                % First 2 entries are size and number of entries in motor map.
                newSlow_J1_fwd_vmap(3:end) = Slow_J1_fwd_vmap(3:end) * (1-J1posSlope.slow.a);
                newSlow_J2_fwd_vmap(3:end) = Slow_J2_fwd_vmap(3:end) * (1-J2posSlope.slow.a);
                % Save new map to config structure
                CobraConfig = setMMap(CobraConfig,pId,'SLOW',1,'fwd',newSlow_J1_fwd_vmap);
                CobraConfig = setMMap(CobraConfig,pId,'SLOW',2,'fwd',newSlow_J2_fwd_vmap);
            end
        end
        
        if revUpdate
            if fastUpdate
                % Apply damping correction
                % First 2 entries are size and number of entries in motor map.
                newFast_J1_rev_vmap(3:end) = Fast_J1_rev_vmap(3:end) * (1-J1negSlope.fast.a);
                newFast_J2_rev_vmap(3:end) = Fast_J2_rev_vmap(3:end) * (1-J2negSlope.fast.a);
                % Save new map to config structure
                CobraConfig = setMMap(CobraConfig,pId,'FAST',1,'rev',newFast_J1_rev_vmap);
                CobraConfig = setMMap(CobraConfig,pId,'FAST',2,'rev',newFast_J2_rev_vmap);
            end
            if slowUpdate
                % Apply damping correction
                % First 2 entries are size and number of entries in motor map.
                newSlow_J1_rev_vmap(3:end) = Slow_J1_rev_vmap(3:end) * (1-J1negSlope.slow.a);
                newSlow_J2_rev_vmap(3:end) = Slow_J2_rev_vmap(3:end) * (1-J2negSlope.slow.a);
                % Save new map to config structure
                CobraConfig = setMMap(CobraConfig,pId,'SLOW',1,'rev',newSlow_J1_rev_vmap);
                CobraConfig = setMMap(CobraConfig,pId,'SLOW',2,'rev',newSlow_J2_rev_vmap);
            end
        end
    end
    
    % Create figure for current positioner
    figure('name',sprintf('PID%d',pId)); 
    
    % Grab positive and negative ranges of joint request data
    PrangeJ1 = [0 max(data1.J1.req)];
    NrangeJ1 = [min(data1.J1.req) 0];
    PrangeJ2 = [0 max(data1.J2.req)];
    NrangeJ2 = [min(data1.J2.req) 0];

    % Grab all data for 1st moves
    move1 = find(data1.J1.iter == 1);
    %%%% THETA %%%%
    % Plot error vs request with first moves highlighted in green circles
    subplot(211)
    plot(data1.J1.req, data1.J1.err0,'.'); hold on;
    plot(data1.J1.req(move1), data1.J1.err0(move1),'go');% first move position
    
    % Plot forward and reverse slopes for fast and slow maps on theta
    plot(PrangeJ1, J1posSlope.fast.a * PrangeJ1, 'r');
    plot(NrangeJ1, J1negSlope.fast.a * NrangeJ1, 'r');
    plot(PrangeJ1, J1posSlope.slow.a * PrangeJ1, 'r--');
    plot(NrangeJ1, J1negSlope.slow.a * NrangeJ1, 'r--');
    xlabel('request angle [rad]'); ylabel('error angle [rad]');
    title(sprintf('err vs. req J1; R_f:%s F_f:%s',...
                  format_data(J1negSlope.fast.a,J1negSlope.fast.sigma),...
                  format_data(J1posSlope.fast.a,J1posSlope.fast.sigma)));
    refline(0,0,0,'k');refline(0,0,Inf,'k');
    hold off;
    
    %%%% PHI %%%%
    % Plot error vs request with first moves highlighted in green circles
    subplot(212)
    plot(data1.J2.req, data1.J2.err0,'.'); hold on;
    plot(data1.J2.req(move1), data1.J2.err0(move1),'go'); % first move position
    
    % Plot forward and reverse slopes for fast and slow maps on theta
    plot(PrangeJ2, J2posSlope.fast.a * PrangeJ2, 'r');
    plot(NrangeJ2, J2negSlope.fast.a * NrangeJ2, 'r');
    plot(PrangeJ2, J2posSlope.slow.a * PrangeJ2, 'r--');
    plot(NrangeJ2, J2negSlope.slow.a * NrangeJ2, 'r--');
    xlabel('request angle [rad]'); ylabel('error angle [rad]');
    title(sprintf('err vs. req J2; R_f:%s F_f:%s',...
                  format_data(J2negSlope.fast.a,J2negSlope.fast.sigma),...
                  format_data(J2posSlope.fast.a,J2posSlope.fast.sigma)));
    refline(0,0,0,'k');refline(0,0,Inf,'k');
    hold off;

    % Check matlab version before saving figure
    if verLessThan('matlab','8')
      disp('modifyDamping (155): image not saved');
    else
      savefig(sprintf('PID%d_err_vs_req.fig',pId))
    end

    % Plot new motor maps
    % FAST MAPS
    fhs(jj) = figure('name', sprintf('ARM_DATA_%d', pId));
    lhm(1) = plotMotMapsFromStruct(CobraConfig.ARM_DATA.(armfld).FAST_CALIBRATION_TABLE,'r');
    lhm(2) = plotMotMapsFromStruct(CobraConfigOrig.ARM_DATA.(armfld).FAST_CALIBRATION_TABLE,'r--');
    lhm(3) = plotMotMapsFromStruct(CobraConfig.ARM_DATA.(armfld).SLOW_CALIBRATION_TABLE,'b');
    lhm(4) = plotMotMapsFromStruct(CobraConfigOrig.ARM_DATA.(armfld).SLOW_CALIBRATION_TABLE,'b--');
    subplot(2,2,1)
    lhmh = [lhm.ph];
    legend(lhmh(1:4:end),'New Fast Map','Old Fast Map','New Slow Map','Old Slow Map')
    
    % Assing new data1 structure to q#
    assignin('base',sprintf('q%d',jj),data1);
    
    % Increment data set counter
    jj=jj+1;

    % Write motor model variables to log file
    try
      txtout = sprintf('%d: %.2f %.2f %.2f %.2f | ',...
                       pId,...
                       data1.J1r.alpha, data1.J1r.beta,...
                       data1.J2r.alpha, data1.J2r.beta);
    catch
      txtout = sprintf('%d: NaN NaN NaN NaN | ', pId);
    end
    fprintf(1,txtout);
    fprintf(alpha_beta_fid,txtout);

    % Write err/req slope to log file
    try
      fprintf(slope_fid,'%d,%.2f,%.2f,%.2f,%.2f\n',...
              pId,J1posSlope.fast.a,J1negSlope.fast.a,J2posSlope.fast.a,J2negSlope.fast.a);
    catch
      fprintf(slope_fid,'%d,NaN,NaN,NaN,NaN\n',pId);
    end
end
fclose(alpha_beta_fid);
fclose(slope_fid);
fprintf(1,'\n');

