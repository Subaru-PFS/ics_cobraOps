function output = S2CentersAnalysis_CIT(imgsFile, dataTitle,CobraConfig,mID)

% Load the image structures file
load(imgsFile);

% Create color wheel
CW = colormap(lines);

% Determine if a current set of plots exists from another instance of this
% function
try
    figure(findobj('type','figure','name','Fiducial Spreads'))
    % count number entries in legend already
    [LEGH,OBJH,OUTH,OUTM] = legend;
    Ntraces = length(OUTM);
    disp(OUTM)
    
    % Set next color for this instance of function
    dataColor = CW(Ntraces+1,:);
    
catch err
    dataColor = CW(1,:);
end

% fixDataStruct;
cfit = 'cfit0';


%% Process data structures
% Query all variables in the current workspace
S = whos;
centImgs = regexp(cellstr([S.name]),'centerms2fw_ImageId_\d*','match');
centImgs = sortOnRegexp(centImgs,'(?<=ImageId_)\d*');
refImgs = regexp(cellstr([S.name]),'refs2_ImageId_\d*','match');
refImgs = sortOnRegexp(refImgs,'(?<=ImageId_)\d*');
fidPosRaw = [];
centerPosCCD = [];
centersCMS = [];
Link1s = [];
c=1;
dnames = {};
refpos = cellfun(@str2num,regexp(CobraConfig.ARM_DATA.refpos.Text,' ','split'));
allFidCMS = [];

hf1 = existingFigure('Center Spreads');
clf;
hf2 = existingFigure('Fiducial Spreads');

% Extract positions from data structures
for vv = 1:length(centImgs{1})
    centimgname = char(centImgs{1}(vv));
    thiscentimg = eval(centimgname);
    refimgname = char(refImgs{1}(vv));
    thisrefimg = eval(refimgname);
    %  try
    disp(centimgname)
    disp(refimgname)
    
    dnames{c} = centimgname;
    
    % Check that active cobra arrays have multiple centroids
    data.activeCobras = thiscentimg.activeCobras;
    for ii=thiscentimg.activeCobras
        fldID=sprintf('pId%d',ii);
        if length(thiscentimg.(fldID).CCDpos) == 1
            
        end
    end
    
    % Fit circles and find centers of positioners
    thiscentimg = getS2Center(thiscentimg,true);
    
    % FIDUCIAL STUFF
    %Nfid = length(thisrefimg.fixedCobras);
    Nfid = length(thisrefimg.fiducials)
    activeCobras = thisrefimg.activeCobras;
   
    for ii=1:length(thisrefimg.fiducials)
        fldID=sprintf('pId%d',thisrefimg.fiducials(ii));
        if isfield(thisrefimg,fldID) % @Johannes added for running it with less than all fixed cobras. j
            
            if length(thisrefimg.(fldID).CCDpos)>0
                % keyboard;
                fldID;
                fidPosRaw(c,ii) = thisrefimg.(fldID).CCDpos(1);
            else
                disp(strcat('No fiducal data on : ', fldID));
                 fidPosRaw(c,ii) = 0 + 0i;
            end
        end
    end
     
    if(~ismember(0, fidPosRaw(c,:))) % Do not use the images if refimage doesnt have fiducal data. Assuming fiducials are never at 0.
        %Calculate horn transform
        hh = horn87(fidPosRaw(c,:), refpos);
        maxNum = length(thisrefimg.fiducials);
        % Apply horns to fiducials
        for ii=1:maxNum
            figure(hf2)
            subplot(2,ceil(maxNum/2),ii);
            
            fldID=sprintf('pId%d',thisrefimg.fiducials(ii));
            thisrefimg.(fldID).CMSpos = thisrefimg.(fldID).CCDpos * exp(i*hh.R) * hh.S + hh.T;
            plot(thisrefimg.(fldID).CMSpos,'x','MarkerEdgeColor',dataColor);
            hold on
            %cmplx(@text,thisrefimg.(fldID).CMSpos,num2str(c));
            %keyboard;
            if(~isempty(thisrefimg.(fldID).CMSpos))
                allFidCMS(c,ii) = thisrefimg.(fldID).CMSpos;
            end
            axis equal
        end
        
        % Apply horns to centers
        cf=1;
        CW3 = colormap(hsv(length(thiscentimg.activeCobras)));
        
        maxNum = length(thiscentimg.activeCobras);
        for ii= thiscentimg.activeCobras
            figure(hf1)
            subplot(3,ceil(maxNum/3),cf)
            hold on
            fldID=sprintf('pId%d',ii);
            tcenter = getARMval(CobraConfig,ii,mID,'Global_base_pos_x') + ...
                getARMval(CobraConfig,ii,mID,'Global_base_pos_y')*i;
            if(isfield(thiscentimg.(fldID), 'cfit0'))   % @Johannes this check was added since  there must not be a center for every pid anymore.
                %       keyboard
                thiscentimg.(fldID).(cfit).CMSc = thiscentimg.(fldID).(cfit).c * exp(i*hh.R) * hh.S + hh.T;
                lh = plot(thiscentimg.(fldID).(cfit).CMSc,'x','MarkerEdgeColor',dataColor);
                centersCMS(c,cf) = thiscentimg.(fldID).(cfit).CMSc;
                Link2s(c,cf) = thiscentimg.(fldID).(cfit).R;
                Link1s(c,cf) = abs(thiscentimg.(fldID).(cfit).CMSc - tcenter);
                hold on
                axis equal
            else
                centersCMS(c,cf) = 0 + 0 * 1i;
            end
            
            cf=cf+1;
        end
    end
    c= c+1;
end

fidBar = [];
for ii=1:Nfid
    figure(hf2)
    subplot(2,ceil(Nfid/2),ii)
    fidBar(ii) = mean(allFidCMS(:,ii));
    fidSTD(ii) = std(allFidCMS(:,ii));
    title(['fiducial ' num2str(ii)])
    ch1 = cmplx(@circle,fidBar(ii),3*fidSTD(ii),'Color',dataColor);
    appendLegend(ch1,strcat(dataTitle{1},' 3\sigma circle (\sigma=',num2str(fidSTD(ii)),')'))
end

Nact = length(activeCobras);
xwidth = zeros(1,Nact);
ywidth = zeros(1,Nact);
figure(hf1)
for ii=1:Nact
    subplot(3,ceil(Nact/3),ii)
    
     centersCMSlocal = centersCMS(:,ii);
    centersCMSlocal(centersCMSlocal == 0) = [];
   if(~isempty(centersCMSlocal))
    cBar(ii) = mean(centersCMSlocal);
    cSTD(ii) = std(centersCMSlocal)./sqrt(2);
    
    
    pIdstr = ['pid' num2str(activeCobras(ii))];
    
    
    title(sprintf('PID%d',activeCobras(ii)))
    
    
    title(pIdstr)
    if (cSTD(ii) ~= 0)
        ch2 = cmplx(@circle,cBar(ii),3*cSTD(ii),'Color',dataColor);
        nm4lgnd = char(regexp(strrep(dataTitle{1},'_','-'),'.*(?=\.mat)','match'));
        hl = appendLegend(ch2,strcat(nm4lgnd,' 3\sigma circle (\sigma=',num2str(cSTD(ii)),')'));
        
        % Get axis limits and scale
        v = axis;
        xwidth(ii) = v(2)-v(1);
        ywidth(ii) = v(4)-v(3);
    end
   end
end

%     figure(hf1)
%     for ii = 1:Nact
%         subplot(3,ceil(Nact/3),ii)
%         v = axis;
%         xW = max(xwidth);
%         yW = max(ywidth);
%         xM = v(1)+(v(2)-v(1))/2;
%         yM = v(3)+(v(4)-v(3))/2;
%         xlim([xM-xW/2,xM+xW/2])
%         ylim([yM-yW/2,yM+yW/2])
%     end

output.dnames = dnames;
output.cBar = cBar;
output.cSTD = cSTD;
output.fidBar = fidBar;
output.fidSTD = fidSTD;
output.Link1s = Link1s;
output.Link2s = Link2s;
output.activeCobras = activeCobras;



%% Plot Link Lengths
hf3 = existingFigure('Link Lengths');
for ii=1:size(Link1s,2)
    
    figure(hf3)
    % Link 1
    Link1DataPos = Link1s(:,ii);
    Link1DataPos(Link1DataPos==0)=[];
    if(~isempty(Link1DataPos))
    subplot(Nact,2,2*(ii-1)+1)
    %     plot(Link1s(:,ii),'*','Color',dataColor)
    
    hist(Link1DataPos,50);
    hold on
    ylabel('px')
    % Find std and mean
    stdStr = sprintf('pId%d_stdL1',activeCobras(ii));
    assignin('caller',stdStr,std(Link1DataPos))
    meanStr = sprintf('pId%d_meanL1',activeCobras(ii));
    assignin('caller',meanStr,mean(Link1DataPos))
    
    title(['pId' num2str(activeCobras(ii)) ' Link1']);
    nm4lgnd = char(regexp(strrep(dataTitle{1},'_','-'),'.*(?=\.mat)','match'));
    end
      Link2DataPos = Link2s(:,ii);
    Link2DataPos(Link2DataPos==0)=[];
    % Link 2
     if(~isempty(Link2DataPos))
    subplot(Nact,2,2*ii)
    hist(Link2DataPos,50);
    hold on
    ylabel('px')
    % Find std and mean
    stdStr = sprintf('pId%d_stdL2',activeCobras(ii));
    assignin('caller',stdStr,std(Link2DataPos))
    meanStr = sprintf('pId%d_meanL2',activeCobras(ii));
    assignin('caller',meanStr,mean(Link2DataPos))
    
    title(['pId' num2str(activeCobras(ii)) ' Link2']);
    
    nm4lgnd = char(regexp(strrep(dataTitle{1},'_','-'),'.*(?=\.mat)','match'));
     end
end


end
function output = sortOnRegexp(arr, expStr)
IDs = regexp(arr{1},expStr,'match');
for ii=1:length(IDs)
    Num(ii) = str2num(IDs{ii}{1});
end
[sorted Isort] = sort(Num);
arr{1} = arr{1}(Isort);
output = arr;
end

