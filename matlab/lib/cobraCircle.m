function output=cobraCircle(data,debugplots)

CW = colormap(lines(length(data.activeCobras)*2));
nCbr = 1;
for ii=data.activeCobras

    
    fldID=sprintf('pId%d',ii);
    
    
    % Check if there is CMS data and use it if so, else use raw
    % Store whichever position array exists in tpos array to be used hereon
    if isfield(data.(fldID),'CMSpos')
        tpos = data.(fldID).CMSpos;
    elseif isfield(data.(fldID),'pos')
        tpos = data.(fldID).pos;
    elseif isfield(data.(fldID),'RawPos')
        tpos = data.(fldID).RawPos;
    elseif isfield(data.(fldID),'CCDpos')
        tpos = data.(fldID).CCDpos(1,:);
    end

    
    

    
    
    %% first pass circle fit on raw points
    cfit0 = circfit(tpos);
    
    %% First level analysis of points
    % Find X,Y points in local cobra frame based on first pass fit
    lev1.xy = tpos - cfit0.c;
    % Calculate radii and smooth them out
    lev1.R  = smooth(abs(lev1.xy),5);
    % Calculate theta angles at each point
    lev1.theta = angle(lev1.xy);
    % Sort theta in ascending order
    [lev1.theta, I] = sort(lev1.theta,'ascend');
    % Sort radii vector accordingly
    Rtemp = [];
    for ll = 1:length(lev1.theta)
        Rtemp(ll) = lev1.R((I(ll)));
    end
    lev1.R = Rtemp;
    % Unwrap theta and cat it with itself to encompass > 2*pi range. 
    lev1.theta = unwrap(lev1.theta);
    lev1.theta = [lev1.theta lev1.theta+2*pi];
    % Cat radii vector with itself to match theta vector
    lev1.R = [Rtemp Rtemp];
    
    %% Second level analysis with interpolation
    % Decide what the interval shall be for interpolated points
    %     dtheta = abs(mean(diff(angle(lev1.xy))));
    dtheta = 10*pi/180;
    % Create theta vector for interpolation
	lev2.theta = subarray(linspace(min(lev1.theta),min(lev1.theta)+2*pi,ceil(2*pi/dtheta)+1),'1:end-1');
    % Select interpolation method    
    interpmethod = 'pchip';    
      
	% Interpolate radii for lev2.theta
    lev2.R = interp1(lev1.theta, lev1.R, lev2.theta ,interpmethod, 'extrap');
    % Calculate lev2 xy points in cobra local frame
    lev2.xy = lev2.R .* exp(i*lev2.theta);
    % Fit circle to new points translated back to CCD frame   
    cfit2 = circfit(lev2.xy + cfit0.c);
    

    %% Calculate residuals of circle fits
    resid0 = abs(tpos - cfit0.c) - cfit0.R;
    resid2 = abs(tpos - cfit2.c) - cfit2.R;
    
    
    %% Debug plots
    if debugplots
        try
%         figID = 20000+ii;
%         
%         while ishandle(figID)
%             figID= figID+1;
%         end
%         
%         if ishandle(figID)
%             figure(figID)
%             close(gcf) 
%         end
%         figure(figID)
%         title(data.Name,'Interpreter','none')
%         hold on
        
        color1 = CW(nCbr,:);
        color2 = CW(nCbr+length(CW(:,1))/2,:);
        nCbr = nCbr+1;
        
        figure(500+nCbr)
%         title(data.Name,'Interpreter','none')
        hold on
        plot(lev1.theta,lev1.R,'Color',color1);
        plot(lev2.theta,lev2.R,'Color',color2);
%         legend(lh,'Raw','iner')
%         xlabel('theta (deg)')
%         ylabel('R (px)')
        
        
%         hold on
%         plot(tpos-data.origin,'rx');
%         plot(tpos,'bo');
%         circle(cfit0.xc,cfit0.yc,cfit0.R,'r')
%         circle(cfit1.xc,cfit1.yc,cfit1.R,'b')
%         hold off
%         axis equal
%         hgsave(strcat(data.CntrImgsDir,data.Name,'_',num2str(figID),'.fig'));
%         
%         figID = 30000+ii;
%         if ishandle(figID)
%             figure(figID)
%             close(gcf) 
%         end
%         figure(figID)
%         hold on
%         plot(unwrap(angle(tpos-cfit1.c)),resid1,'r')
%         plot(unwrap(angle(tpos(ok)-cfit2.c)),resid2,'b')
%         hold off
%         hgsave(strcat(data.CntrImgsDir,data.Name,'_',num2str(figID),'.fig'));
        catch err
            disp(err)
        end
    end
    
    % Store results in data structure
    data.(fldID).cfit0 = cfit0;
    data.(fldID).cfit2 = cfit2;
    data.(fldID).resid0 = resid0;
    data.(fldID).resid2 = resid2;
end

output = data;

return 
