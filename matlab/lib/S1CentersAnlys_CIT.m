function output = S1CentersAnalysis_CIT_(imgsFile, dataTitle, CobraConfig, cobraMap)

% Load the image structures file
load(imgsFile);

% Create color wheel
CW = colormap(lines);

% Pixel Scale
PxSc = 90.5;

% Determine if a current set of plots exists from another instance of this
% function
try
    figure(findobj('type','figure','name','Center Spreads'))
    % count number entries in legend already
    [LEGH,OBJH,OUTH,OUTM] = legend;
    Ntraces = length(OUTM);
    
    % Set next color for this instance of function
    dataColor = CW(Ntraces+1,:);
    
catch err
    dataColor = CW(1,:);
end

cfit = 'cfit2';



%% Process data structures
% Query all variables in the current workspace
S = whos;

centImgs = regexp(cellstr([S.name]),'centerms1_ImageId_\d*','match');
centImgs = sortOnRegexp(centImgs,'(?<=ImageId_)\d*');
refImgs = regexp(cellstr([S.name]),'refs1_ImageId_\d*','match');
refImgs = sortOnRegexp(refImgs,'(?<=ImageId_)\d*');
fidPosRaw = [];
centerPosCCD = [];
centersCMS = [];
c=1;
dnames = {};
allFidCMS = [];
refpos = cellfun(@str2num,regexp(CobraConfig.ARM_DATA.refpos.Text,' ','split'));

hf1 = existingFigure('Center Spreads');
hf2 = existingFigure('Fiducial Spreads');

% Extract positions from data structures
for vv = 1:length(centImgs{1})
    centimgname = char(centImgs{1}(vv));
    thiscentimg = eval(centimgname);
    refimgname = char(refImgs{1}(vv));
    thisrefimg = eval(refimgname);
    %     try
    disp(centimgname)
    disp(refimgname)
    
    dnames{c} = centimgname;
    
    % Check that active cobra arrays have multiple centroids 
    data.activeCobras = thiscentimg.activeCobras;
    for ii=thiscentimg.activeCobras
        fldID=sprintf('pId%d',ii);
        if length(thiscentimg.(fldID).CCDpos) == 1
            disp(strcat('Not enough data in: ', thiscentimg.Name));
            disp(strcat('For Cobra: ',fldID));
        end
    end
 
    % Fit circles and find centers of positioners 
    thiscentimg = getS1Center(thiscentimg,true);
    
    % FIDUCIAL STUFF  if isfield(Ptemp,pfld)
    %Nfid = length(thisrefimg.fixedCobras)
    Nfid = length(thisrefimg.fiducials)
    for ii=1:Nfid
        fldID=sprintf('pId%d',thisrefimg.fiducials(ii));
        if isfield(thisrefimg,fldID) % @Johannes added for running it with less than all fixed cobras.
            fidPosRaw(c,ii) = thisrefimg.(fldID).CCDpos;
        else
            disp(strcat('No data on this? ', fldID));
        end
    end
    
    %Calculate horn transform
    hh = horn87(fidPosRaw(c,:), refpos);
     
    % Apply horns to fiducials
    maxNum = length(thisrefimg.fiducials);
    for ii=1:maxNum
        figure(hf2)
        subplot(2,ceil(maxNum/2),ii);
        
        fldID=sprintf('pId%d',thisrefimg.fiducials(ii));
        thisrefimg.(fldID).CMSpos = thisrefimg.(fldID).CCDpos * exp(i*hh.R) * hh.S + hh.T;
        plot(thisrefimg.(fldID).CMSpos,'o','MarkerEdgeColor',dataColor);
        plot(thisrefimg.(fldID).CCDpos,'rx')
        hold on
        %             cmplx(@text,thisrefimg.(fldID).CMSpos,num2str(c));
        allFidCMS(c,ii) = thisrefimg.(fldID).CMSpos;
        axis equal
    end
     
    %  keyboard;
    % Apply horns to centers
    cf=1;
    maxNum = length(thiscentimg.activeCobras);
    for ii= thiscentimg.activeCobras
        figure(hf1)
        subplot(3,ceil(maxNum/3),cf) 
        hold on
        fldID=sprintf('pId%d',ii);
        %                 lh = plot(thisvar.(fldID).(cfit).c,'bx');
        if(isfield(thiscentimg.(fldID), 'cfit0'))   % @Johannes this check was added since  there must not be a center for every pid anymore.
            thiscentimg.(fldID).(cfit).CMSc = thiscentimg.(fldID).(cfit).c * exp(i*hh.R) * hh.S + hh.T;
            lh = plot(thiscentimg.(fldID).(cfit).CMSc,'x','MarkerEdgeColor',dataColor);
            %             cmplx(@text,thiscentimg.(fldID).(cfit).CMSc,num2str(c));
            centersCMS(c,cf) = thiscentimg.(fldID).(cfit).CMSc;
            axis equal
        end
        cf=cf+1; 
    end 
    c=c+1;
    
  
%     figure(222)
%     
%     I=fitsread(fullfile(thiscentimg.ImgDir,thiscentimg.imageName));
%     imagesc(I)
%     hold on
%      for ii= thiscentimg.activeCobras
%           fldID=sprintf('pId%d',ii);
%            if(isfield(thiscentimg.(fldID), 'cfit2')) 
%          plot(thiscentimg.(fldID).(cfit).CMSc,'rx')
%            end
%      end
end

fidBar = [];
for ii=1:Nfid
    figure(hf2)
    subplot(2,ceil(Nfid/2),ii)
    fidBar(ii) = mean(allFidCMS(:,ii));
    fidSTD(ii) = std(allFidCMS(:,ii))./sqrt(2);
    title(['fiducial ' num2str(ii)])
    ch1 = cmplx(@circle,fidBar(ii),3*fidSTD(ii),'Color',dataColor);
    nm4lgnd = char(regexp(strrep(dataTitle{1},'_','-'),'.*(?=\.mat)','match'));
    hl = appendLegend(ch1,sprintf('3\\sigma circle (\\sigma=%4.2fum)',PxSc*fidSTD(ii)));   
end


 

Nact = size(centersCMS,2);
xwidth = zeros(1,Nact);
ywidth = zeros(1,Nact);
figure(hf1)
for ii=1:Nact % Plot
    subplot(3,ceil(Nact/3),ii)
    centersCMSlocal = centersCMS(:,ii);
    centersCMSlocal(centersCMSlocal == 0) = [];
   if(~isempty(centersCMSlocal))
    cBar(ii) = mean(centersCMSlocal);
    cSTD(ii) = std(centersCMSlocal)./sqrt(2);
    
    if exist('cobraMap','var')
        title(sprintf('SN%d Center Measurements',cobraMap(data.activeCobras(ii))))
    else
        title(sprintf('PID%d',data.activeCobras(ii)))
    end
    % keyboard;
    
    if (cSTD(ii) ~= 0)
        ch2 = cmplx(@circle,cBar(ii),3*cSTD(ii),'Color',dataColor);
        nm4lgnd = char(regexp(strrep(dataTitle{1},'_','-'),'.*(?=\.mat)','match'));
        hl = appendLegend(ch2,sprintf('3\\sigma circle (\\sigma=%4.2fum)',PxSc*cSTD(ii)));
        
        % Get axis limits and scale
        v = axis;
        xwidth(ii) = v(2)-v(1);
        ywidth(ii) = v(4)-v(3);
    else
    end
   end
    
end

% figure(hf1)
% for ii = 1:Nact
%     subplot(3,ceil(Nact/3),ii)
%     v = axis;
%     xW = max(xwidth);
%     yW = max(ywidth);
%     xM = v(1)+(v(2)-v(1))/2;
%     yM = v(3)+(v(4)-v(3))/2;
%     xlim([xM-xW/2,xM+xW/2])
%     ylim([yM-yW/2,yM+yW/2])
% end

output.dnames = dnames;
output.cBar = cBar;
output.cSTD = cSTD;
output.fidBar = fidBar;
output.fidSTD = fidSTD;
output.activeCobras = data.activeCobras;
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

