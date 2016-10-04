function output=motorMap(data,nStep,stage,CobraConfig,direction)


%% Prerequisites: 
    % Data centroids must been processed and run through cobraCMS function
    % Data must have been run through the cobraCircle function
    % Data must have been taken with first points coming off the hardstop
    

%% Motor map generation
HM = figure('Name','Centroids');
plotFibers(data,'bx');
hold on;
c=1;
for ii=data.activeCobras
    
    %Clear temp vars
    clear I tpos cpos J HS tempPos tempcPos
    
    % Create field ID string
    fldID=sprintf('pId%d',ii);
    
    % Get HS vector from initPos array
    % For stage 1 this is position vector wrt cobra center
    % For stage 2 this is position vector wrt elbow position which should
    % not move during the test since it is not commanded.
    
    
    switch stage
        case 1
            clocking = CobraConfig.(fldID).orientation*pi/180; % CW from +X
            HS = cos(clocking) + sin(clocking)*i; % Unity Hardstop vector
            start = 1;
            maxRom = 400;
        case 2
            start = 2;
            P1 = data.(fldID).CCDpos(1,1) - data.(fldID).cfit0.c;
            clocking = angle(P1);
            clocking(clocking<0) = clocking(clocking<0) + 2*pi;
            HS = cos(clocking) + sin(clocking)*i; % Unity Hardstop vector
            maxRom = 200;
            plot(data.(fldID).cfit0.c + P1,'ro')
    end
    
    Thetas = [];
    Phis = [];
    CPs = [];
    Ps = [];
    imgN = [];
%     fh8 = figure('Name','Theta testing');
%     switch stage
%         case 1
            % Declare maximum range of motion for theta
            
            
            % For each centroid
            for ii=start:length(data.(fldID).CCDpos)
                % Calculate local cobra position vector to centroid
                P1 = data.(fldID).CCDpos(1,ii) - data.(fldID).cfit0.c;
                
                if stage==1
                    % Calculate arm angles. Note theta is CW wrt +X
                    [theta1, phi, alpha] = cobraArmAngles(P1, CobraConfig.(fldID).L1,CobraConfig.(fldID).L2);
                    % Rotate P CW by alpha to line up with theta arm
                    P1 = P1*exp(i*alpha);
                end
                
                % Normalize P
                P = P1/abs(P1);
                        
                
                figure(HM);
                plot([data.(fldID).cfit0.c,data.(fldID).cfit0.c+90*HS])
                
                theta = acos(dot([real(P),imag(P)],[real(HS),imag(HS)]));
                imgN = [imgN data.(fldID).CCDpos(2,ii)];
                CP = cross([real(P),imag(P),0],[real(HS),imag(HS),0]);
                
                % If this is the first image
                if data.(fldID).CCDpos(2,ii) == 1
                    plot(data.(fldID).cfit0.c + P1,'ro')
                    
                % If this is the second image
                elseif data.(fldID).CCDpos(2,ii) == 2
                    plot(data.(fldID).cfit0.c + P1,'go')

                end
                    
                Ps = [Ps data.(fldID).cfit0.c + P1];
                CPs = [CPs CP(3)];
                Thetas = [Thetas theta];
%                 Phis = [Phis phi];
            end

            % Theta matrix
            thetaM1 = [Thetas',imgN',CPs',Ps.'];
            
            
            switch stage
                case 1
                    % Subtract from 360deg theta angles which are from 2nd image
                    % and pos cross product
                    thetaM2 = thetaM1;
                    thetaM2(thetaM2(:,2)==2 & thetaM2(:,3)>0,1) = 2*pi - thetaM2(thetaM2(:,2)==2 & thetaM2(:,3)>0,1);

                    % Subtract from 360deg IMG1 theta between 45deg and 180deg with
                    % pos cross product
                    thetaM2(thetaM2(:,2)==1 & thetaM2(:,3)>0 & thetaM2(:,1)>pi/4,1) = 2*pi - thetaM2(thetaM2(:,2)==1 & thetaM2(:,3)>0 & thetaM2(:,1)>pi/4,1);

                    % Add 360 to IMG2 theta between 0 and 90deg with neg xproduct
                    thetaM2(thetaM2(:,2)==2 & thetaM2(:,3)<0 & thetaM2(:,1)<pi/2,1) = 2*pi + thetaM2(thetaM2(:,2)==2 & thetaM2(:,3)<0 & thetaM2(:,1)<pi/2,1);
                    
                
                case 2
                    thetaM2 = thetaM1;
                    % Subtract from 360deg IMG1 theta between 90deg and 180deg with
                    % pos cross product
                    thetaM2(thetaM2(:,2)==1 & thetaM2(:,3)>0 & thetaM2(:,1)>pi/2,1) = 2*pi - thetaM2(thetaM2(:,2)==1 & thetaM2(:,3)>0 & thetaM2(:,1)>pi/2,1);
            end
                    
            thetaM3 = sortrows(thetaM2,1);
            
            for ii=1:length(Ps)
%                 cmplx(@text,thetaM3(ii,4),num2str(round(thetaM3(ii,1)*180/pi)));
                cmplx(@text,thetaM3(ii,4),num2str(ii));
            end
            
            
            % Get first image centroids and joint angles
%             [J tpos] = jointsForMotorMap(1,data,fldID,stage,armLengths(c,:),direction,clocking)
%             [JResid tposResid] = jointsForMotorMap(2,data2,fldID,stage,armLengths(c,:),direction,clocking)
            
%             switch direction
%                 case 'FWD'
%                     J = [J JResid];
%                     tpos = [tpos tposResid];
%                 case 'REV'
%                     J = [JResid J];
%                     tpos = [tposResid tpos];
%             end
            J = thetaM3(:,1);
            data.(fldID).J1 = J;
                    
%             
%         case 2
%             % Negative hardstop CW wrt upper stage home position
% %             [unused clocking] = cobraArmAngles(HS, armLengths(c,1), armLengths(c,2));
% %             maxRom = 220;
%     end

   
    
    
    %% Create Motor Map

%     angles = angle(tpos - data.(fldID).cfit0.c)
%     dPsi = [diff(unwrap((angle(tpos - data.(fldID).cfit0.c))))]
    dJ = diff(J);
    
    
    % Initialize motor map array
    mMap = [];
    % Initialize join angle tracker. This is the current angle from hardstop
    Jt = 0;
    % Initialize region angle tracker. Integrates J during loop until
    % reset by new region condition
    cJ = 0;
    % Initialize region counter
    cRegion = 1;
    % Initialize step counter
    cStep = 0;
    % How many moves form a region?
    moveN = 1;
    % Iterate across all motions
    for ii=1:length(dJ)
        % Incriment counters
        cJ = cJ+dJ(ii);
        Jt = Jt+dJ(ii);
        cStep = cStep+nStep(c);
        % Check if new region threshold has been exceeded
        if mod(ii,moveN)==0 || ii==length(dJ)
            % Store current region motor map data
            mMap(cRegion,:) = [abs(Jt*180/pi()),abs(cJ*180/pi()/cStep)];
            % Increment region and reset counters
            cRegion = cRegion + 1;
            cStep = 0;
            cJ = 0;
        end
    end
    
    % Extend last motor map region to end of ROM
    mMap(cRegion,:) = [maxRom,mMap(cRegion-1,2)];
    
    %% Store motor map in data structure
    data.(fldID).mMap = mMap';
    data.(fldID).negHS = clocking;
%     data.(fldID).posHS = anglePosHS;
    
    c=c+1;
   
end

output = data;
return