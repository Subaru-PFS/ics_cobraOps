% load config xml
clear all
close all

 cfgPath = '/Users/johannes/Dropbox/PFS_EM/TEST_RESULTS/TestDataAnalysis';
CobraConfig = loadCfgXml(cfgPath);

Npid = length(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER);

% GROUP C

%activeCobras = [3,5,11,15,17,23,27,33,35,39,41,47,53];
%highStill = ([1:Npid/3+1]-1)*6+3;
%lowStill = ([1:Npid/3+1]-1)*6+5;

% %GROUP B
% activeCobras = [5,11,17,23,29,35,41,47,53];
% highStill = ([1:Npid/3+1]-1)*6+1;
% lowStill = ([1:Npid/3+1]-1)*6+3;

% %GROUP A
 activeCobras = [3,9,15,21,27,33,39,45,51,57];
 highStill = [1:Npid/3+1]*6-1;
 lowStill = [1:Npid/3+1]*6+1;

phiAngle = 45;
thtAngleStill = [10,350,60,300,110,240];
thtAngleActv = 0;


% Create new target list files
for ii=1:Npid
%     Get PID
    pid = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii}.DATA_HEADER.Positioner_Id.Text);
%     Create files for non-active cobras
	if sum(activeCobras == pid) ~= 1
        tfile = fopen(sprintf('TargetList_mId_%d_pId_%d.txt',1,pid),'w');
        fclose(tfile);
    end
    
    
end

% Generate targets for each sharp edge angle
for jj=1:length(thtAngleStill)
%     Start new figure for current sharp edge angle
    figure('Name',['Sharp edge @ ',num2str(thtAngleStill(jj))]);
    
%     Pull out centers from database and plot them with PID labels
    for ii=1:Npid
        xx = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii}.KINEMATICS.Global_base_pos_x.Text);
        yy = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii}.KINEMATICS.Global_base_pos_y.Text);
        pid(ii) = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii}.DATA_HEADER.Positioner_Id.Text);

        centers(ii) = xx + yy*i;

        plot(xx + yy*i,'o');
        text(xx,yy,num2str(pid(ii)))
        hold on
    end
    
%     For only the stationary cobras, generate sharp edge targets
    for ii=1:Npid-1
        if ii > 1
            % Starboard stationary cobras
            if sum(pid(ii) == highStill) == 1

                % Get link length
                l1length = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii-1}.KINEMATICS.Link1_Link_Length.Text);
                % Find vector to neighboring active cobra and rotate by theta
                % angle
                l1vec = (centers(ii-1) - centers(ii))*exp(i*(thtAngleStill(jj))*pi/180);
                refDir = [centers(ii),centers(ii-1)];
                l1 = l1vec/abs(l1vec)*l1length;
                l1line = [centers(ii),centers(ii)+l1];

                % Get link 2 length
                l2length = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii-1}.KINEMATICS.Link2_Link_Length.Text);
                % Rotate link1 by phi angle to make link2 direction
                l2dir = l1*exp(i*((phiAngle-180)*pi/180));
                % Normalize and scale to link2 length
                l2 = l2dir/abs(l2dir)*l2length;
                % Create target at end of link2
                target = l1line(2)+l2;
                l2line = [l1line(2),target];

                tfile = fopen(sprintf('TargetList_mId_%d_pId_%d.txt',1,pid(ii)),'a+');
                fprintf(tfile,'%4.2f,%4.2f,0,0\n',real(target),imag(target));
                fclose(tfile);

                plot(l1line,'-r')
                plot(l2line,'-b')
                plot(refDir,'--r')
                plot(centers(ii),'+r')

            % Port Stationary Cobras
            elseif sum(pid(ii) == lowStill) == 1

                l1length = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii}.KINEMATICS.Link1_Link_Length.Text);
                l1vec = (centers(ii+1) - centers(ii))*exp(i*(thtAngleStill(jj))*pi/180);
                refDir = [centers(ii),centers(ii+1)];
                l1 = l1vec/abs(l1vec)*l1length;
                l1line = [centers(ii),centers(ii)+l1];

                l2length = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii}.KINEMATICS.Link2_Link_Length.Text);
                l2dir = l1*exp(i*((phiAngle-180)*pi/180));
                l2 = l2dir/abs(l2dir)*l2length;
                target = l1line(2)+l2;
                l2line = [l1line(2),target];

                tfile = fopen(sprintf('TargetList_mId_%d_pId_%d.txt',1,pid(ii)),'a+');
                fprintf(tfile,'%4.2f,%4.2f,0,0\n',real(target),imag(target));
                fclose(tfile);


                plot(l1line,'-g')
                plot(l2line,'-b')
                plot(refDir,'--g')
                plot(centers(ii),'+g')
                plot(real(target),imag(target),'+')
            end
        end


        axis equal
    end
end
    

%   pids = [data.tht.pid];
%   
%   if pidORpnum < 0
%     pnum = -pidORpnum;
%   else
%     pnum = find(pids == pidORpnum);
%   end
%   %% THETA
% 
%     ZZ = data.tht(pnum).centroids - data.tht(pnum).center;
%     angles = angle(ZZ);
%     radii  = abs(ZZ);
%     ok = data.tht(pnum).ok;
% 
%     figure(10)
%     axis([0 2000 0 2000]);
%     daspect([1,1,1]);
%     hold on;
%     plot([data.tht.center],'rx')
%     plot([data.tht.centroids], 'bo');
%     plot(data.tht(pnum).centroids,'gx');
%     plot([data.phi.center],'rx')
%     plot([data.phi.centroids], 'bo');
%     plot(data.phi(pnum).centroids,'gx');
%     hold off;
%     xlabel('X [pix]');
%     ylabel('Y [pix]');
%     title(sprintf('\\Theta and \\Phi(pId %d)',data.tht(pnum).pid));
% 
%     figure(11) 
%     plot(angles/(2*pi),radii, 'bx');
%     hold on;
%     plot(angles(ok)/(2*pi), radii(ok),'go');%,'MarkerFace','g');
%     hold off;
%     xlabel('angle/\tau');
%     ylabel('radius [pix]');
%     title(sprintf('\\Theta PID %d',data.tht(pnum).pid));
% 
%     figure(12)
%     hist(radii(ok),50);
%     xlabel('radius [pix]');
%     title(sprintf('\\Theta PID %d',data.tht(pnum).pid));
%     
%     %% PHI    
%     ZZ = data.phi(pnum).centroids - data.phi(pnum).center;
%     angles = angle(ZZ);
%     radii  = abs(ZZ);
%     ok = data.phi(pnum).ok;
% 
% % $$$     figure(20)
% % $$$     axis([0 2000 0 2000]);
% % $$$     daspect([1,1,1]);
% % $$$     hold on;
% % $$$     plot([data.phi.center],'rx')
% % $$$     plot([data.phi.centroids], 'bo');
% % $$$     plot(data.phi(pnum).centroids,'gx');
% % $$$     hold off;
% % $$$     xlabel('X [pix]');
% % $$$     ylabel('Y [pix]');
% % $$$     title(sprintf('\\Phi (pId %d)',data.phi(pnum).pid));
%     
%     figure(21) 
%     plot(angles/(2*pi),radii, 'bx');
%     hold on;
%     plot(angles(ok)/(2*pi), radii(ok),'go');%,'MarkerFace','g');
%     hold off;
%     xlabel('angle/\tau');
%     ylabel('radius [pix]');
%     title(sprintf('\\Phi PID %d',data.phi(pnum).pid));
% 
%     figure(22)
%     hist(radii(ok),50);
%     xlabel('radius [pix]');
%     title(sprintf('\\Phi Pid %d',data.tht(pnum).pid));
%     return;
%     
