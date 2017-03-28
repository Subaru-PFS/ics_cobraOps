function [output strArr] = processTargetLog(targetLogFile,pIdSTR,CobraConfig,CCT,saveFigs,figPrefix)


%% Target Log File Format Reference
% Table 4 - Target Mode Output Log Files
% Column	Name	Description
% 1	Target #	Starts with 0
% 2	Iteration #	Starts with 0 as the first image of initialized Cobra before any moves
% 3	mId	module ID
% 4	pId	positioner ID
% 5	X_t	X coordinate of the target position (in pixels)
% 6	Y_t	Y coordinate of the target position (in pixels)
% 7	J1_t	Joint 1 target angle (in degrees)
% 8	J2_t	Joint 2 target angle (in degrees)
% 9	-	Zero
% 10	-	Zero
% 11	X_c	X coordinate of the current position (in pixels)
% 12	Y_c	Y coordinate of the current position (in pixels)
% 13	J1_c	Joint 1 current angle (in degrees)
% 14	J2_c	Joint 2 current angle (in degrees)
% 15	-	Zero
% 16	-	Zero
% 17	Distance	Distance to target (in microns)
% 18	J1_Steps	Joint 1 steps remaining – calculated based on motor calibration matrix
% 19	J2_Steps	Joint 2 steps remaining – calculated based on motor calibration matrix
% 20	-	Zero
% 21	-	Zero
% 22	Status	1=converged, 0=not converged

% maximum number of iterations per convergence
maxIters = 11;
 

%% EXECUTION CODE
% Read in target log file and config file
numIt = 10;%Number of iterations viewed.

TL = dlmread(targetLogFile,',');
savePath = targetLogFile(1:find(targetLogFile=='\',1,'last'));

% Initialize data structure
data.targetIds = [];
data.cvrgIters = [];

% iteration counter
lastIterN = -10;

tt = 0;
lastId = -1;
% Populate data structure from file contents
for ii = 1:length(TL(:,1))
    
    % Get target ID from file
    targId = TL(ii,1);
    
    % Get iteration number
    iterN = TL(ii,2);
    
    % Check if the target has changed
    if targId ~= lastId || iterN~=lastIterN+1
        % Incriment target counter
        tt = tt+1;
        
        % Initialize target arrays 
        targStr = ['target' num2str(tt)];
        data.targetIds = [data.targetIds tt];
        data.(targStr).pid = [TL(ii,4)];
        data.(targStr).iter = [];
        data.(targStr).status = [];
        data.(targStr).dist = [];
        data.(targStr).targCmplx = [TL(ii,5) + TL(ii,6)*i];
        data.(targStr).J1 = [];
        data.(targStr).J2 = [];
        data.(targStr).J1_t = []; %[TL(ii,7)];
        data.(targStr).J2_t = []; %[TL(ii+1,8)];
%         data.(targStr).J1_tx = CobraConfig.(pIdSTR).L1*cos(2*pi-TL(ii,7)*pi/180);
%         data.(targStr).J1_ty = CobraConfig.(pIdSTR).L1*sin(2*pi-TL(ii,7)*pi/180);
%         data.(targStr).J2_tx = CobraConfig.(pIdSTR).L2*cos( 2*pi-TL(ii,7)*pi/180 + (pi - TL(ii,8)*pi/180) );
%         data.(targStr).J2_ty = CobraConfig.(pIdSTR).L2*sin( 2*pi-TL(ii,7)*pi/180 + (pi - TL(ii,8)*pi/180) );
%         data.(targStr).Jx = data.(targStr).J1_tx + data.(targStr).J2_tx;
%         data.(targStr).Jy = data.(targStr).J1_ty + data.(targStr).J2_ty;
        data.(targStr).J1rgn = [];
        data.(targStr).J2rgn = [];
        data.(targStr).J1err = [];
        data.(targStr).J2err = [];
        data.(targStr).J1_S = [];
        data.(targStr).J2_S = [];
        data.(targStr).curPos = []; % Current position in cobra frame
%         data.(targStr).RdTh = []; % Vector from Cobra center to target
    end
    
    lastIterN = iterN;
    
    % Update last target tracker
    lastId = targId;
    
    % Add this data to target 
    targStr = ['target' num2str(tt)];  
    data.(targStr).iter = [data.(targStr).iter; (TL(ii,2)+1)];
%     data.(targStr).iter = [data.(targStr).iter; (TL(ii,2))];
    data.(targStr).status = [data.(targStr).status; TL(ii,22)];
    data.(targStr).dist = [data.(targStr).dist; TL(ii,17)];
    data.(targStr).J1 = [data.(targStr).J1; TL(ii,13)];
    data.(targStr).J2 = [data.(targStr).J2; TL(ii,14)];
    
    data.(targStr).J1_t = [data.(targStr).J1_t; TL(ii,7)];
    data.(targStr).J2_t = [data.(targStr).J2_t; TL(ii,8)];
    
    data.(targStr).J1_S = [data.(targStr).J1_S; TL(ii,18)];
    data.(targStr).J2_S = [data.(targStr).J2_S; TL(ii,19)];
 %keyboard
    if(iterN == 0)
    data.(targStr).J1err = [data.(targStr).J1err; data.(targStr).J1_t(iterN+1) - TL(ii,13)];
    data.(targStr).J2err = [data.(targStr).J2err; data.(targStr).J2_t(iterN+1) - TL(ii,14)];
    else 
    data.(targStr).J1err = [data.(targStr).J1err; data.(targStr).J1_t(iterN) - TL(ii,13)];
    data.(targStr).J2err = [data.(targStr).J2err; data.(targStr).J2_t(iterN) - TL(ii,14)];
    end
    data.(targStr).curPos = [data.(targStr).curPos; TL(ii,11) + TL(ii,12)*i]; % Current position in cobra frame
%     data.(targStr).RdTh = [data.(targStr).RdTh abs(data.(targStr).curPos-CobraConfig.(pIdSTR).s1Center) * (data.(targStr).J1_t - TL(ii,13)) ];
  
end

% Which fields are in degs?
degFlds = {'J1','J2','J1err','J2err','J1_t','J2_t'};

%% Analyze Data structure
CW = colormap(lines(length(data.targetIds)));
nCW=1;
Z = [];
unZ = [];
ConvN = 0;
ConvI = [];
unConvN = 0;
unConvI = [];
FinErr = [];
FinJ1 = [];
FinJ2 = [];
ConvIMG = zeros(maxIters,length(data.targetIds));
J1_T = [];
J2_T = [];
J1_z = [];
J2_z = [];
unJ1_z = [];
unJ2_z = [];
CCTI = [];
JointsIMG = zeros(200,400);
hh1 = figure('Name',strcat(pIdSTR,' Convergence Analysis'));
% hh5 = figure;

for ii=data.targetIds
    targStr = ['target' num2str(ii)];
    data.(targStr).cvrgIter = data.(targStr).iter(find(data.(targStr).status,1));
    data.cvrgIters = [data.cvrgIters data.(targStr).cvrgIter];
    
    % Pull out the final distance to target at last iteration
    FinErr = [FinErr data.(targStr).dist(end)];
    FinJ1 = [FinJ1 data.(targStr).J1(end)];
    FinJ2 = [FinJ2 data.(targStr).J2(end)];
    
    % Add target joint angles to master arrays
    J1_T = [J1_T data.(targStr).J1_t];
    J2_T = [J2_T data.(targStr).J2_t];
    
    % Declare custom convergence threshold iteration to zero
    data.(targStr).CCTI = 0;
    
    % Look for CCT to be met
    for uu=1:length(data.(targStr).dist)
        % Add value to convergence image
        ConvIMG(uu,ii) = data.(targStr).dist(uu);
        
        % Check if CCT is met and update CCTI if so
        if data.(targStr).CCTI==0 && data.(targStr).dist(uu) < CCT
            data.(targStr).CCTI = data.(targStr).iter(uu);
        end
        
    end
    
    % Add this target's CCTI to master array
    CCTI = [CCTI data.(targStr).CCTI];
    
%     % Plot convergence traces in RdTh space
%     figure(hh5)
%     plot(data.(targStr).RdTh)
%     hold on

    % Plot convergence trace
    X_vector = [0:length(data.(targStr).dist)-1];
    figure(hh1)
    subplot(2,3,1)
    semilogy(X_vector,data.(targStr).dist,'Color',CW(nCW,:));
    hold on
    
    figure(hh1)
    subplot(2,3,2)
    plot(X_vector,(data.(targStr).J1err),'Color',CW(nCW,:));
    hold on
    subplot(2,3,5)
    plot(X_vector,(data.(targStr).J2err),'Color',CW(nCW,:));
    hold on
   
    % Increment color wheel counter
    nCW = nCW+1;
    
    % Add targets to Z matrix
    if data.(targStr).CCTI > 0
        Z = [Z;
            real(data.(targStr).targCmplx),imag(data.(targStr).targCmplx),data.(targStr).CCTI];
        J2_z = [J2_z data.(targStr).J2_t];
        J1_z = [J1_z data.(targStr).J1_t];
        ConvN = ConvN+1;
        ConvI = [ConvI data.(targStr).CCTI];
    else
        unZ = [unZ;
            real(data.(targStr).targCmplx),imag(data.(targStr).targCmplx),data.(targStr).CCTI];
        unJ2_z = [unJ2_z data.(targStr).J2_t];
        unJ1_z = [unJ1_z data.(targStr).J1_t];
        unConvN = unConvN+1;
        unConvI = [unConvI data.(targStr).CCTI];
    end
    
    
    % Add data to structure array
    flds = fields(data.(targStr));
   
    for kk=1:length(flds)
%         disp(flds(kk))
        fldlength = length(data.(targStr).(flds{kk}));
        
        if fldlength > 1
          try
          strArr(ii).(flds{kk}) = padarray(data.(targStr).(flds{kk}),[maxIters-fldlength, 0],nan,'post');
          catch
            disp(['at Line 227 of processTargetLog.m.  Probably need to ' ...
                  'increase maxIters. Or you did not get a image toolbox license!']);
            keyboard
          end
        else
          strArr(ii).(flds{kk}) = data.(targStr).(flds{kk});
        end

        % Convert degree fields to radians
        if sum(strcmp(flds{kk},degFlds))
          strArr(ii).(flds{kk}) = strArr(ii).(flds{kk}) * pi / 180.;
        end
        
    end 
end

%% Get Motor Map Moves 

%keyboard;




%% Create Figures

figure(hh1)
subplot(2,3,2)
title('Theta Error')
subplot(2,3,5)
title('Phi Error')

%% Convergence Traces
figure(hh1)
subplot(2,3,1)
hl1 = semilogy(ones(1,numIt)*CCT,'--g','LineWidth',3);
legend(hl1,strcat(num2str(CCT),'um Convergence'))
title('Convergence to target')
ylabel('Dist to Target (um)')
xlabel('Iteration')

figure(hh1)
subplot(2,3,2)
xlabel('Iterations')
ylabel('J1 Error (deg)')
subplot(2,3,5)
xlabel('Iterations')
ylabel('J2 Error (deg)')



%% Cummulative Convergence Plot
ConvP = [];
for itr = 1:numIt
    if itr == 1
        ConvP(itr) = sum(ConvI == itr)/(ConvN+unConvN)*100;
    else
        ConvP(itr) = ConvP(itr-1) + sum(ConvI == itr)/(ConvN+unConvN)*100;
    end
end
X_vector3 = [0:length(ConvP)-1];
figure(hh1)
subplot(2,3,4)
l1 = plot(X_vector3,ConvP,'-bo');
hold on
l2 = plot(95*ones(1,numIt),'--g');
title(strcat(num2str(CCT),'um Cumulative Convergence'))
ylabel('Percentage Convergence (%)')
xlabel('Iterations')
set(gca,'YTick',linspace(0,100,11))
set(gca,'XTick',linspace(1,numIt,numIt))
grid on
ylim([0 100])
legend([l1 l2], '% Converged', '95%','location','SE')

data.convP = ConvP;



%% Iterations vs Joint Angles
% figure(hh1)
% 
% subplot(2,3,3)
% plot(J1_z,ConvI,'o')
% if ~isempty(J1_z)
%     
%     keyboard
%     J1_xx = linspace(min(J1_z),max(J1_z),100);
%     J1_PP = polyfit(J1_z,ConvI,3);
%     J1_yy = polyval(J1_PP,J1_xx);
% 
%     hold on
%     plot(J1_xx,J1_yy,'--r')
%     % plot(J1_fit,'--k')
%     xlabel('J1 Angle (deg)')
%     ylabel('Iterations')
%     title(strcat('Iterations to reach < ',num2str(CCT),'um'))
% end
% 
% 
% subplot(2,3,6)
% plot(J2_z,ConvI,'o')
% if ~isempty(J2_z)
%     J2_xx = linspace(min(J2_z),max(J2_z),100);
%     J2_PP = polyfit(J2_z,ConvI,3);
%     J2_yy = polyval(J2_PP,J2_xx);
%     hold on
%     plot(J2_xx,J2_yy,'--r')
%     xlabel('J2 Angle (deg)')
%     ylabel('Iterations')
% end

if saveFigs
    saveas(hh1,char(strcat(savePath,figPrefix,pIdSTR,'_anlys.fig')))
end


%% Convergence Image
% hh3 = figure('Name',strcat(pIdSTR,' Convergence Image'));
% 
% % Define color scale upper limit
% clrLim = 2000; %microns
% logLim = 10; % log(um)
% 
% % Plot the Convergence Image matrix
% % imagesc(ConvIMG,[0,clrLim])
% imagesc(log(abs(ConvIMG)),[0,logLim])
% colormap(pink)
% 
% % Create colorbar and update tick marks to data values
% cb1 = colorbar;
% cbt = get(cb1, 'YTick');
% set(cb1, 'YTickLabel', sprintf('%0.1f|', exp(cbt)))
% ylabel(cb1,'microns')
% 
% % Label plot
% title('Target Distance')
% xlabel('Target #')
% ylabel('Iteration')
% 
% % Flip the y-axis of the image
% set(gca,'YDir','normal')
% 
% data.convIMG = ConvIMG;
% 
% if saveFigs
%     saveas(hh3,char(strcat(savePath,figPrefix,pIdSTR,'_cnvImg.fig')))
% end



% %% Final Distance vs Joint Angle
% hh4 = figure('Name',strcat(pIdSTR,' Final Distance vs Joint Angle'));
% subplot(2,1,1)
% plot(FinJ1,FinErr,'ro','MarkerFaceColor','r')
% grid on
% title('Stage 1')
% xlabel('Final Joint Angle (Deg)')
% ylabel('Final Error (um)')
% subplot(2,1,2)
% plot(FinJ2,FinErr,'ro','MarkerFaceColor','r')
% grid on
% title('Stage 2')
% xlabel('Final Joint Angle (Deg)')
% ylabel('Final Error (um)')
% 
% if saveFigs
%     saveas(hh4,char(strcat(savePath,figPrefix,pIdSTR,'_jDist.fig')))
% end



% %% Patrol Region Difficulty
% hh5 = figure('Name',strcat(pIdSTR,' Patrol Region Difficulty'));
% 
% pId = regexp(pIdSTR,'\d*','match');
% cnfgID = ['ARM_DATA_' pId{1}];
% s1Center = getARMval(CobraConfig,pId{1}, 1,'Global_base_pos_x')+getARMval(CobraConfig,pId{1}, 1,'Global_base_pos_y')*i;
% 
% % s1Center = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_x.Text) + ...
% %                     str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_pos_y.Text)*i;
%   
% R_ptrl = getARMval(CobraConfig,pId{1}, 1,'Link1_Link_Length')+getARMval(CobraConfig,pId{1}, 1,'Link2_Link_Length');
% 
% % R_ptrl = str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_Link_Length.Text) + ...
% %     str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_Link_Length.Text);
% 
% orientation = getARMval(CobraConfig,pId{1}, 1,'Global_base_ori_z');
% 
% %str2num(CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_ori_z.Text))
% orient = (360-orientation)*pi/180; %in matlab angle rad
% circle(real(s1Center),imag(s1Center),R_ptrl,'b');
% hsVect = R_ptrl*(cos(orient)+sin(orient)*i);
% hold on
% plot([s1Center,s1Center+hsVect],'k--')
% if ~isempty(Z)
%     scatter(Z(:,1),Z(:,2),40,Z(:,3),'filled');%
% end
% if ~isempty(unZ)
%     scatter(unZ(:,1),unZ(:,2),40,unZ(:,3));
% end
% t = colorbar;
% set(get(t,'ylabel'),'String', 'Iterations');
% caxis([0 10]);
% axis equal
% title(strcat('Iterations to reach < ',num2str(CCT),' um from target'))
% xlabel('px')
% ylabel('px')
% if saveFigs
%     saveas(hh5,char(strcat(savePath,figPrefix,pIdSTR,'_ptrl.fig')))
% end
%keyboard
output = data;


return;








        
