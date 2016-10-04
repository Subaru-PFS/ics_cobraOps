%FLIP COBRA CONFIG

activeCobras = [1:4 6:9];


for ii = activeCobras
    
    cnfgID = sprintf('ARM_DATA_%d',ii);
%     CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_ori_z = 360 - CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_ori_z;
%     CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_y = 2048 - CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_y;
%     CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Pixel_scale = 90;
%     CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_BaseRadious = 2;
%     
%     CobraConfig.ARM_DATA.(cnfgID).DATA_HEADER.Module_Id = 1;
%     CobraConfig.ARM_DATA.(cnfgID).DATA_HEADER.Positioner_Id = ii;
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint1_CW_limit_angle = CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint1_CW_limit_angle + 5;
    

end

% %REMOVE COBRAS
% CobraConfigOLD = CobraConfig;
% 
% CobraConfig.ARM_DATA = struct;
% for ii=activeCobras
%     cnfgID = sprintf('ARM_DATA_%d',ii);
%     CobraConfig.ARM_DATA.(cnfgID) = CobraConfigOLD.ARM_DATA.(cnfgID);
% end
% 
% configSavePath = 'CobraConfig_061014_2_inv.mat';
% save(configSavePath,'CobraConfig')