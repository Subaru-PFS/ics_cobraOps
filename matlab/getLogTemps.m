function output = getLogTemps(msimLogFile)

% logFilePath = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_16_14_15_20_40_Targets\output_log.txt';

logF = fopen(msimLogFile);

% Get initial line
tline = fgetl(logF);

% Line filters
filedate = regexp(msimLogFile,'\d*_\d*_\d*_\d*_\d*_\d*','match');
fltSvImg = strcat('seq_saveCurrent_RemoteImage.*',filedate);
fltTemp = 'Get_Temp';

% Initiate vars
rexTrg = [];
trgts = [];
sensor = [];
ambient = [];
lookTemp = false;



while ischar(tline)
    
    % Look for temperature line if flag raised
    if lookTemp
        rexTemp = regexp(tline,fltTemp);
        if ~isempty(rexTemp)
            ttemp = str2num(char(regexp(tline,'(?<=(Sensor=)|(ambient=))\d*\.\d*','match')));
            sensor = [sensor ttemp(1)];
            ambient = [ambient ttemp(2)];
            lookTemp = false; % Lower look flag for temp
%             disp(tline)
        end
    end
        
    
    % Look for save image line with target ID
    rexTrg = regexp(tline,fltSvImg);
    if ~isempty(rexTrg{1})
        trgts = [trgts str2num(char(regexp(tline,'(?<=Target_)\d*','match')))];
        lookTemp = true; % Raise flag to look now for temp
%         disp(tline)
    end
    
    % Get next line
    tline=fgetl(logF);
end
% 
fclose(logF);

output.targets = trgts;
output.sensor = sensor;
output.ambient = ambient;

end