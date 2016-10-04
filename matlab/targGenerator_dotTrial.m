%TARGET GENERATOR
clear all
close all

%% INPUTS
% Get config xml
CobraConfig = loadCfgXml;
% CobraConfig = loadCfgXml('C:\Users\sage\Desktop\Dropbox\PFS_EM\TEST_RESULTS\Metrology','LinkLengths_Centers_062614_2.xml');

% Specify which cobras to create targets for
moduleID = 1;
activeCobras = [3:9];

% Specify number of random targets to generate
numtrg = 1;

reachAround = 0.05; % Specify reach factor. Eg: 0.05 means that 5% of 
% reachable radial distance is excluded for target generation so that there
% is some padding on ideal keepout zones.


%% EXECUTION
fh1 = figure;
for ii = activeCobras
    cnfgID = sprintf('ARM_DATA_%d',ii);
    
    tfile = fopen(sprintf('TargetList_mId_%d_pId_%d.txt',moduleID,ii),'w');
    
    center = str2num(CobraConfig.ARM_DATA.(['ARM_DATA_' num2str(ii)]).KINEMATICS.Global_base_pos_x.Text) + ...
                    str2num(CobraConfig.ARM_DATA.(['ARM_DATA_' num2str(ii)]).KINEMATICS.Global_base_pos_y.Text)*i;
    L1 = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_Link_Length.Text);
    L2 = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_Link_Length.Text);
    hsAngle = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_ori_z.Text);
    phiMax = str2num(char(CobraConfig.ARM_DATA.(cnfgID).phiMax.Text));
    phiMin = str2num(char(CobraConfig.ARM_DATA.(cnfgID).phiMin.Text));
    
    % Inner keepout radius
    Riko = sqrt(L1^2+L2^2-2*L1*L2*cos(phiMin*pi/180));
    % Max reach
    reach = sqrt(L1^2+L2^2-2*L1*L2*cos(phiMax*pi/180));
    
    
    plot(center,'b+');
    hold on;
    cmplx(@circle,center,Riko,'m')
    cmplx(@circle,center,reach,'r')
    
    for ii=1:numtrg
        % Pick a random angle between 0 and 2pi
        gamma = 60 + hsAngle - 17.5;
        % gamma = rand*2*pi;
        % Pick a random radii reachable by arm
        Rtrg = sqrt(L1^2+L2^2-2*L1*L2*cos(60*pi/180));
        % pad = reachAround*reach;
        % Rtrg = Riko + pad + rand*(reach-Riko-2*pad);
        
        
        
        % Calculate x and y components of target from cobra center
        cxT = Rtrg*cos(gamma);
        cyT = Rtrg*sin(gamma);
        
        % Add cobra center to get target x,y in fiducial frame
        target = center + cxT + cyT*i;
        
        figure(1)
        plot(target,'ro','MarkerFaceColor','red','MarkerSize',1.5)
        cmplx(@text,target,num2str(ii));
        
        if ii==numtrg
            fprintf(tfile,'%4.2f,%4.2f',real(target),imag(target));
        else
            fprintf(tfile,'%4.2f,%4.2f\n',real(target),imag(target));
        end
    end
    
%     target = center + L1 + L2*i;
%     fprintf(tfile,'%4.2f,%4.2f',real(target),imag(target));
    
    fclose(tfile);

end

axis equal
