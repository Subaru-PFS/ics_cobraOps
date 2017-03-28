function output = getTargetList(fieldid)
% Retreive the targets for a given field from test data (mat
% files?)
%
% usage: output = getTargetList(fieldID)
% where fieldID counts from 0.
    
    fieldid = fieldid + 1;
    
    c = loadCfgXml();
    s = loadTestData('',c,1);
    
    ncobras = length(c.ARM_DATA.ARM_DATA_CONTAINER);
    targets = nan(ncobras,1) + i * nan(ncobras,1);
    
    for jj=1:ncobras
        pid = str2num(c.ARM_DATA.ARM_DATA_CONTAINER{jj}.DATA_HEADER.Positioner_Id.Text);
        if ~isempty(s(pid).xtg)
            targets(jj) = s(pid).xtg(fieldid) + 1i *s(pid).ytg(fieldid);
        else
            % use the same-sense home position as the target location.
            ck = c.ARM_DATA.ARM_DATA_CONTAINER{jj}.KINEMATICS;
            cen  = str2num(ck.Global_base_pos_x.Text) + ...
                   str2num(ck.Global_base_pos_y.Text) * i ;
            tht0 = str2num(ck.CCW_Global_base_ori_z.Text)*pi/180;
            phi0 = max( str2num(ck.Joint2_CCW_limit_angle.Text)*pi/180 - pi, ...
                        -pi ) + 0.1;
            L1   = str2num(ck.Link1_Link_Length.Text);
            L2   = str2num(ck.Link2_Link_Length.Text);
            targets(jj) = cen + + L1 .* exp(i*tht0) + L2 .* exp(i*(tht0 + phi0));

        end
    end
    
    output = targets;

