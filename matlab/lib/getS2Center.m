function output=getS2Center(data,debugplots)

CW = colormap(lines(length(data.activeCobras)*2));
nCbr = 1;
for ii=data.activeCobras

    
    fldID=sprintf('pId%d',ii);
    CCDpos = data.(fldID).CCDpos;
        if(length(data.(fldID).CCDpos) > 2) % @Johannes

    %% first pass circle fit on raw points
    cfit0 = circfit(CCDpos);
 
	%% First level analysis of points
    % Find X,Y points in local cobra frame based on first pass fit
    lev1.xy = CCDpos - cfit0.c;
    % Calculate radii and smooth them out
    lev1.R  = smooth(abs(lev1.xy),5);
    % Calculate theta angles at each point
    lev1.theta = unwrap(angle(lev1.xy));
    % Sort theta in ascending order
    [lev1.theta, I] = sort(lev1.theta,'ascend');
    % Sort radii vector accordingly
    Rtemp = [];
    for ll = 1:length(lev1.theta)
        Rtemp(ll) = lev1.R((I(ll)));
    end
    lev1.R = Rtemp;


    %% Calculate residuals of circle fits
    resid0 = abs(CCDpos - cfit0.c) - cfit0.R;
    
    
    
    
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
        
        
%         figure(500)
%         subplot(2,3,nCbr)
% %         title(data.Name,'Interpreter','none')
%         plot(lev1.theta/pi,lev1.R,'.:','Color',color1);
%         hold on
%         title(data.activeCobras(nCbr))
%         xlabel('angle (n\pi)')
%         ylabel('Link2 (px)')
        
        nCbr = nCbr+1;
%         legend(lh,'Raw','iner')
%         xlabel('theta (deg)')
%         ylabel('R (px)')
        
        
%         hold on
%         plot(CCDpos-data.origin,'rx');
%         plot(CCDpos,'bo');
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
%         plot(unwrap(angle(CCDpos-cfit1.c)),resid1,'r')
%         plot(unwrap(angle(CCDpos(ok)-cfit2.c)),resid2,'b')
%         hold off
%         hgsave(strcat(data.CntrImgsDir,data.Name,'_',num2str(figID),'.fig'));
        catch err
            disp(err)
        end
    end
    
    % Store results in data structure
    data.(fldID).cfit0 = cfit0;
    data.(fldID).resid0 = resid0;
end
end
output = data;

return 
