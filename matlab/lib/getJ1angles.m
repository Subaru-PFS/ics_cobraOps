function output = getJ1angles(data,CobraConfig)

% There should only be one fiber image per positioner per image for this
% function to work.

for ii = data.activeCobras
    fldID = strcat('pId',num2str(ii));
    cnfgID = sprintf('ARM_DATA_%d',ii);
    % Initialize Joint1 array
    data.(fldID).J1 = [];
    data.(fldID).J2 = [];
    data.(fldID).leftArm = false;
    tcenter = CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_x + ...
                    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_y*i;
                figure(1)
                plot(tcenter,'go')
    L1 = CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_Link_Length;
    L2 = CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_Link_Length;

    
    
    pos = data.(fldID).CMSpos;
    for jj = 1:length(pos)
        % Calculate pos in Cobra's frame
        
        cbraPos = pos(jj) - tcenter;
        
        % CW wrt X-axis
        J1Angle = i;
        J2Angle = J1Angle;
        cc = 0;
        
        
        while ~isreal(J1Angle) && ~isreal(J2Angle)
            [J1Angle J2Angle] = cobraArmAngles(cbraPos, L1, L2, data.(fldID).leftArm);
            linkMod = 1.05;
            L2 = linkMod*L2;
            cc = cc+1;
            if cc>1
                disp(sprintf('increasing link2 length for pid%d',ii))
            end
        end
        
        J1Angle = 2*pi-J1Angle;
        J2Angle = 2*pi-J2Angle;
        
        data.(fldID).J1 = [data.(fldID).J1 J1Angle];
        data.(fldID).J2 = [data.(fldID).J2 J2Angle];
%         data.(fldID).psi = [data.(fldID).psi angle(cbraPos)];
%         if data.(fldID).psi < 0
%             data.(fldID).psi = data.(fldID).psi + 2*pi;
%         end
%         % Convert to MSIM coord sys
%         data.(fldID).psi = 2*pi - data.(fldID).psi;
    end
    
    figure(1)
    tx = real(tcenter) +(L1+L2)*cos(J1Angle);
    ty = imag(tcenter) -(L1+L2)*sin(J1Angle);
    plot([tcenter, tx+ty*i],'k--')
end

output = data;

return;