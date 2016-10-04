

dataDirs = {
% '08_19_14_09_56_44_ringout'
% '08_19_14_09_43_39_ringout'
% '08_19_14_09_37_27_ringout'
% '08_19_14_09_31_13_ringout'
% '08_19_14_09_27_27_ringout'
'08_06_14_16_48_39_msimMaps'
'07_21_14_22_15_30_msimMaps'
'07_21_14_15_19_50_msimMaps'
'07_18_14_16_57_14_msimMaps'
'07_17_14_17_33_35_msimMaps'
'07_16_14_11_52_47_msimMaps'
'07_16_14_09_14_52_msimMaps'
    };
nTests = length(dataDirs);

for kk=1:nTests
    ts = getLogTemps(fullfile('D:\PfsTests',dataDirs{kk},'log','output_log.txt'));
    temp(kk).sensor = ts.sensor;
    temp(kk).ambient = ts.ambient;
    
    disp(max(ts.sensor));
    
    filedate = regexp(dataDirs{kk},'\d*_\d*_\d*_\d*_\d*_\d*','match');
    
    thisdate = datenum(filedate,'mm_dd_yy_HH_MM_SS');
    
    for ll=1:length(ts.sensor)
        ts.time(ll) = thisdate + 1/24/60/6*(ll-1);
    end
    
    temp(kk).time = ts.time;
    
end


% % Initiate vars
% rexTrg = [];
% trgts = [];
% sensor = [];
% ambient = [];
% lookTemp = true;
% nTests = length(dataDirs);

% for kk=1:nTests
%     
%     msimLogFile = fullfile(dataDirs{kk},'output_log.txt');
% 
%     logF = fopen(msimLogFile);
% 
%     % Get initial line
%     tline = fgetl(logF);
% 
%     % Line filters
%     % filedate = regexp(msimLogFile,'\d*_\d*_\d*_\d*_\d*_\d*','match');
%     % fltSvImg = strcat('seq_saveCurrent_RemoteImage.*',filedate);
%     fltTemp = 'Get_Temp';
% 
% 
% 
% 
% 
%     while ischar(tline)
% 
%         % Look for temperature line if flag raised
%         if lookTemp
%             rexTemp = regexp(tline,fltTemp);
%             if ~isempty(rexTemp)
%                 ttemp = str2num(char(regexp(tline,'(?<=(Sensor=)|(ambient=))\d*\.\d*','match')));
%                 sensor = [sensor ttemp(1)];
%                 ambient = [ambient ttemp(2)];
%     %             lookTemp = false; % Lower look flag for temp
%     %             disp(tline)
%             end
%         end
% 
% 
%     %     % Look for save image line with target ID
%     %     rexTrg = regexp(tline,fltSvImg);
%     %     if ~isempty(rexTrg{1})
%     %         trgts = [trgts str2num(char(regexp(tline,'(?<=Target_)\d*','match')))];
%     %         lookTemp = true; % Raise flag to look now for temp
%     % %         disp(tline)
%     %     end
% 
%         % Get next line
%         tline=fgetl(logF);
%     end
%     % 
%     fclose(logF);
% end
% 
% output.targets = trgts;
% output.sensor = sensor;
% output.ambient = ambient;
