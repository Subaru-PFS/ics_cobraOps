function output = createMotorMap(motorDirection, oldMap) 
    % keyboard;
    osa = [0 oldMap(1,1:end-1)]*180/pi;% old start angles
    ofa = oldMap(1,:)*180/pi; % old finish angles
    oms = diff([0 oldMap(1,:)])*180/pi; % [old move sizes]
    omv = oldMap(2,1:end);  % [old map values]
    angularLength  = 400;
    bw = 3.6; % bin width in degree
    numBins = (angularLength/bw);
    mma = [0:ceil(numBins)]*bw;
    %mma = linspace(0,angularLength + bw,numBins+2); % [motor map angles]
    msa = motorDirection.startAngle; % [motor start angles]
    msa = vertcat(msa, osa'); % Append motor map values
    sam = bsxfun(@minus, mma, msa); % [start angle matrix]
    
    sel1 = sam>0 & sam < bw; %get just the angular portion in the first bin
    sam2 = sam; % Make a copy to delete all but the first bin values 
    sam2(~sel1)= 0;
    sam2(:,1) =[];   % Use the values to your right (those are the right ones).
    mla = motorDirection.moveSizes; % [move length angles]
    mla = vertcat(mla, oms');

    mfa = motorDirection.finishAngle; % [move final angles]
    mfa = vertcat(mfa, ofa'); % Append motor map values
    fam = bsxfun(@minus, mma, mfa); % [final angle matrix]
    self = fam>0 & fam < bw;  % Select the ending angles
    fam(~self) = 0; % [final angles matrix]
    fam2 = fam(:,2:end);       % Use the values to your right (those are the right ones) 
    msm = sam2; %  [move size matrix]
    
    for ii = 1:size(sam2,1) % Calculate the move size within the bin.
        istart = find(sam2(ii,:) > 0, 1, 'first');
        ifinsh = find(fam2(ii,:) > 0, 1, 'last');
        if(mla(ii) > 0)
            msm(ii,istart+1:ifinsh) = bw;
        else
            msm(ii,ifinsh+1:istart) = bw;
        end
    end
    msm = msm - fam2;
    msm(msm>bw) = msm(msm>bw)-bw; % in the 1-bin case we are overcounting by bin-width (because start fraction and end fraction include the move and the parts to the end on both sides) 
    
    %% Weigthing the moves
    msma = msm.^2; % [move size matrix area] Area of the move 
    mlab = abs(mla).*bw; % [move length angles bin-width] Move Length angles times bin width is the area to divide by 
    msmf = bsxfun(@rdivide, msma, mlab); % [move size matrix fractions] Divide the fractions by the width of the bin and the length of the total move.
    msmf(isnan(msmf)) = 0;
    mmv =  motorDirection.mmap; % [motor map values] 
    mmv = vertcat(mmv, omv'); % add the old motor map
    
    %motorMap = linspace(bw,angularLength, numBins+1);
    motorMap = [1:numBins+1]*bw;
    wmv = bsxfun(@times, msmf, mmv); % [weighted motormap values] Fractions times motormap values
    smmv = sum(wmv,1);  % [sum motor map values] per bin
    sfrac = sum(msmf,1); % [sum fractions] per bin
    nmv = smmv./sfrac; % [new motor-map values] 
    motorMap(2,:) = nmv;
    

    %% Correct Motor Maps:
  %  keyboard;
    %First filter out NANS in case there was no value available in the data
    %or the old map
    nans = find(isnan(motorMap(2,:)));
    mapMean =  mean(motorMap(2,~isnan(motorMap(2,:))));    
    for kk = nans
        if kk>1
            motorMap(2,kk) = motorMap(2,kk-1);    
        else
            motorMap(2,kk) = mean(motorMap(2,~isnan(motorMap(2,:))));    
        end
    end
    
    % Second, get rid of too small values < 0.001 
    motorMap(2,motorMap(2,:)<0.01) = mapMean;


    %% Calculate Variance of the bin:
    % the variance = E(X^2) - (E(X))^2
    % keyboard;
    mmvsq = mmv.^2; % [motor map value squared]
    wmvsq = bsxfun(@times, msmf, mmvsq); % [weighted motor map values (=mean) of squared values (= x^2)]
    smmvsq = sum(wmvsq, 1); % [
    nmvsq = smmvsq./sfrac; % [new motormap values squared]
    nmvsq(isnan(nmvsq)) = 0;
    stdmm = (abs(nmvsq - motorMap(2,:).^2)).^0.5;  % CAUTION abs is used because negative values can appear
    motorMap(3,:) = stdmm;   
    
% figure(23)
% plot(stdmm + nmv)
% hold on
% plot(nmv, 'r')
% plot(nmv - stdmm)    
    
    output = motorMap;
end

function investigateColumn(col, msmb, wmv)
msmb(:,find(msmb(:,13)~=0))
wmv(find(msmb(:,15)~=0),15)
ksc = msmb(find(msmb(:,15)~=0),15);
motr = wmv(find(msmb(:,15)~=0),15);
msmb(msmb(:,13)>0,13);
find(msmb(:,68)==max(msmb(:,68)));

end
    %% Plot Map
%     figure(1)
%     clf;
%     hold on;
%     plot(msa,mmv,'*r') % Move Start Angles
%     plot(mfa,mmv,'*r') % Move Finsh Angles
%     plot([msa, mfa]',[mmv, mmv]','m') % Just the line
%     plot(motorMap(1,:), nmv,'bo-', 'linewidth',3); % new map
%     plot(oldMap(1,:)*180/pi, omv, 'go-', 'linewidth', 3)  % old map
%     legend('new map','old map');
%     axis([0 360 0 0.3]);
%     keyboard;
     