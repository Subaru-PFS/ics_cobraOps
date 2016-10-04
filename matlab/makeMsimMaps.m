function logData = makeMsimMaps(stage,direction,logPrefix,Nsteps,moveMin,movePersist)

% Find the positioner centroid logs in current directory and store pids
posFiles = ls([logPrefix '_mId_*pId*.txt']);
mId = cellfun(@str2num,(regexp(reshape(posFiles',1,[]),'(?<=mId_)\d*','match')'));
mpId = cellfun(@str2num,(regexp(reshape(posFiles',1,[]),'(?<=pId_)\d*','match')'));
pId = (mId-1)*57 + mpId;

% For each log file
for ii=1:length(posFiles(:,1))

    logData(ii).mId = mId(ii);
    logData(ii).pId = pId(ii);
    
    logData(ii).rawData = dlmread(posFiles(ii,:),',');
    
    
    
    
    switch stage
        case 1
            stageCol = 3;
            logData(ii).rgns = logData(ii).rawData(:,stageCol);
            
            % THETA FWD
            if direction == 1
                
                % Unwrap the angles
                logData(ii).rgns = unwrap(logData(ii).rgns*pi/180)*180/pi;
                
                % Shift the dataset over so that first value is 0deg
                logData(ii).hsoffset = logData(ii).rgns(1);
                logData(ii).rgns = logData(ii).rgns - logData(ii).rgns(1);
                
                
            
            % THETA REV
            elseif direction == 2

                % Reverse the array
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
                
                % Unwrap the angles
                logData(ii).rgns = unwrap(logData(ii).rgns*pi/180)*180/pi;
                
                % Shift the dataset over so that first value is 0deg
                logData(ii).rgns = logData(ii).rgns - logData(ii).rgns(1);
                
                % Reverse the array back
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
            end
            
            
            %==STALL PERSIST METHOD==
            % Find indices of good movement data
            goodInd = abs(diff(logData(ii).rgns))>moveMin;
            movecutoff = length(logData(ii).rgns)+1;
            for jj = movePersist:length(goodInd)
                if sum(goodInd(jj-movePersist+1:jj))==0
                    movecutoff1 = jj-movePersist+1;
                    break;
                end
            end
            
            %==REACH FINAL VALUE METHOD==
            lastTheta = logData(ii).rgns(end);
            finalInd = find(abs(logData(ii).rgns-lastTheta) < moveMin);
            movecutoff2 = finalInd(1)+1;
            
            
            if movecutoff2>movecutoff1
                disp('End of map triggered by last value')
                movecutoff = movecutoff2;
            else
                disp('End of map triggered by move persistance')
                movecutoff = movecutoff1;
            end
            
            % Cut bad movement data from array
            logData(ii).rgns = logData(ii).rgns(1:movecutoff-1);

            % Find deltas between each angle and divide by steps to obtain motor
            % constant. 
            logData(ii).dps = abs(diff(logData(ii).rgns)/Nsteps);
            
            
            % Extend last motor constant out one more so that array is same size as
            % angle array
            keyboard;
            logData(ii).dps(end+1) = logData(ii).dps(end);
            
            if direction == 2
                % Reverse the arrays
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
                logData(ii).dps = logData(ii).dps(end:-1:1);
            end
                
            
        case 2
            stageCol = 4;
            logData(ii).rgns = logData(ii).rawData(:,stageCol);
            
            % PHI FWD
            if direction == 1
                
                % Unwrap the angles
                logData(ii).rgns = unwrap(logData(ii).rgns*pi/180)*180/pi;
                
                % Shift the dataset over so that first value is 0deg
%                 logData(ii).rgns = logData(ii).rgns - logData(ii).rgns(1);
            
            % PHI REV
            elseif direction == 2
                % Reverse the array
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
                
                % Unwrap the angles
                logData(ii).rgns = unwrap(logData(ii).rgns*pi/180)*180/pi;
                
                % Shift the dataset over so that first value is 0deg
%                 logData(ii).rgns = logData(ii).rgns - logData(ii).rgns(1);
                
                % Reverse the array back
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
            end

            % Find indices of good movement data
            goodInd = abs(diff(logData(ii).rgns))>moveMin;
            movecutoff = length(logData(ii).rgns)+1;
            for jj = movePersist:length(goodInd)
                if sum(goodInd(jj-movePersist+1:jj))==0
                    movecutoff = jj-movePersist+1;
                    break;
                end
            end
            
            % Cut bad movement data from array
            logData(ii).rgns = logData(ii).rgns(1:movecutoff-1);

            % Flip the array if direction is reverse
            if direction == 2
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
            end
            
%             keyboard;
            
            % Check for reversals indicating flip into left-arm config
            reversals = find(diff(logData(ii).rgns)<0);
            if length(reversals)>0
                logData(ii).rgns = logData(ii).rgns(1:reversals(1));
            end
            
            % Flip the array back if direction is reverse
            if direction == 2
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
            end
    
            
            % Find deltas between each angle and divide by steps to obtain motor
            % constant. 
            logData(ii).dps = abs(diff(logData(ii).rgns)/Nsteps);
            if length(logData(ii).dps) == 0
                fprintf('\n length(logData(ii).dps) == 0');
                keyboard;
            end
            
            % Extend last motor constant out one more so that array is same size as
            % angle array
%             try
            % keyboard
            logData(ii).dps(end+1) = logData(ii).dps(end);
%             catch
%                 keyboard
%             end
            
            if direction == 2
                % Flip the arrays
                logData(ii).rgns = logData(ii).rgns(end:-1:1);
                logData(ii).dps = logData(ii).dps(end:-1:1);
            end
            
            
    end
    
    % Check that rgns array is larger than one
    if length(logData(ii).rgns) <= 1
        % 2) Construct an MException object to represent the error.
        err = MException('ResultChk:BadAngleData', ...
            'Cant continue because angle array is of length one or less');

        disp('Current region array is:')
        logData(ii).rgns

        disp('Raw data for this positioner angles is:')
        logData(ii).rawData(:,stageCol)

        throw(err)
    end
    
    switch stage
        case 1
            logData(ii).rgns(end) = 400;
        case 2
            logData(ii).rgns(end) = 200;
    end
            
     
    

    
end
    


return

