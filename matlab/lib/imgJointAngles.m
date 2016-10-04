function output = imgJointAngles(imgdata,CobraConfig,flipimg,thetaIn)

% There should only be one fiber image per positioner per image for this
% function to work.

if ~exist('flipimg','var'), flipimg=false, end;
if ~exist('thetaIn','var'), thetaIn=false, end;

for ii = imgdata.activeCobras
    
    % Generate strings for struct fields
    fldID = strcat('pId',num2str(ii));
    cnfgID = sprintf('ARM_DATA_%d',ii);
    
    % Initialize Joint arrays
    imgdata.(fldID).J1 = [];
    imgdata.(fldID).J2 = [];
    
    % Assume not left-arm to start
    imgdata.(fldID).leftArm = false;
    
    % Get positioner center from Config
    tcenter = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_x.Text) + ...
                    str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_y.Text)*i;
                
    % Plot center
    existingFigure('Orientation Reference img')
    plot(tcenter,'go')
    hold on
    
    % Get Link lengths
    L1 = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_Link_Length.Text);
    L2 = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_Link_Length.Text);

    
    % Get centroid position 
    if flipimg
%         tcenter = real(tcenter) + (2048-imag(tcenter))*i;
        pos = imgdata.(fldID).CMSpos;
        pos = real(pos) + (2048-imag(pos))*i;
    else
        pos = imgdata.(fldID).CMSpos;
    end
    
    for jj = 1:length(pos)
        % Calculate pos in Cobra's frame
        cbraPos = pos(jj) - tcenter;
        
        % CW wrt X-axis
        J1Angle = i;
        J2Angle = J1Angle;
        cc = 0;
        
        % Calculate arm angles using law of cosines and increase arm
        % lengths slightly if solution is not found
        while ~isreal(J1Angle) && ~isreal(J2Angle)
            [J1Angle J2Angle] = cobraArmAngles(cbraPos, L1, L2, imgdata.(fldID).leftArm, false);
            linkMod = 1.05;
            L2 = linkMod*L2;
            cc = cc+1;
            if cc>1
                disp(sprintf('increasing link2 length for pid%d',ii))
            end
        end
       
        
        
        % Write the angles to structure joint arrays
        imgdata.(fldID).J1 = [imgdata.(fldID).J1 J1Angle];
        imgdata.(fldID).J2 = [imgdata.(fldID).J2 J2Angle];

    end
    
    if length(pos) == 1
        existingFigure('Orientation Reference img')
        % Plot a line from the positioner center representing theta angles
        tx = real(tcenter) +(L1)*cos(J1Angle);
        ty = imag(tcenter) +(L1)*sin(J1Angle);
        if thetaIn
            l1spec = 'r--';
            l2spec = 'g--';
        else
            l1spec = 'g--';
            l2spec = 'r--';
        end
        plot([tcenter, tx+ty*i],l1spec)
        text(tx,ty,sprintf('%0.2f',J1Angle*180/pi),'color',[1 1 1]);

        % Plot a line from elbow to fiber
        plot([tx+ty*i,pos],l2spec)
        text(real(pos),imag(pos),sprintf('%0.2f',J2Angle*180/pi),'color',[1 1 1]);
        
    end
end

output = imgdata;

return;