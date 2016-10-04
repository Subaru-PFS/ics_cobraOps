function output = S1CentersAnalysis(imgsFile, dataTitle,configFile)

% Load the image structures file
load(imgsFile);

% Create color wheel
CW = colormap(lines);

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


% fixDataStruct;
load(configFile)
cfit = 'cfit2';



%% Process data structures
% Query all variables in the current workspace
S = whos;
A = regexp(cellstr([S.name]),'centerms1_ImageId_\d*','match');
fidPosRaw = [];
centerPosCCD = [];
centersCMS = [];
c=1;
dnames = {};
refpos = CobraConfig.refpos;
allFidCMS = [];


hf1 = existingFigure('Center Spreads');
hf2 = existingFigure('Fiducial Spreads');

% Extract positions from data structures
for vv = 1:length(A{1})
    vname = char(A{1}(vv));
    thisvar = eval(vname);
    try
        disp(vname)
        dnames{c} = vname;

        % Scrub fiducial arrays for erroneous centroids
%             thisvar = fixFalseFiducials(thisvar,refpos);

        % Check that active cobra arrays have multiple centroids
        stop = 0;
        data.activeCobras = thisvar.activeCobras;
        for ii=thisvar.activeCobras
            fldID=sprintf('pId%d',ii);
            if length(thisvar.(fldID).CCDpos) == 1
                stop = 1;
            end
        end
        if stop == 1
            continue;
        end


        thisvar = getS1Center(thisvar,true);

        Nfid = length(thisvar.fixedCobras);

        for ii=1:length(thisvar.fixedCobras)
            fldID=sprintf('pId%d',thisvar.fixedCobras(ii));
            fidPosRaw(c,ii) = thisvar.(fldID).CCDpos;
        end

        %Calculate horn transform
        hh = horn87(fidPosRaw(c,:), refpos);
        
%         keyboard

        % Apply horns to fiducials
        for ii=1:length(thisvar.fixedCobras)
            figure(hf2)
            subplot(2,2,ii)

            fldID=sprintf('pId%d',thisvar.fixedCobras(ii));
            thisvar.(fldID).CMSpos = thisvar.(fldID).CCDpos * exp(i*hh.R) * hh.S + hh.T;
            plot(thisvar.(fldID).CMSpos,'o','MarkerEdgeColor',dataColor);
            plot(thisvar.(fldID).CCDpos,'rx')
            hold on
            cmplx(@text,thisvar.(fldID).CMSpos,num2str(c));
            allFidCMS(c,ii) = thisvar.(fldID).CMSpos;
            axis equal
        end

        % Apply horns to centers
        cf=1;
        for ii=thisvar.activeCobras
            figure(hf1)
            subplot(3,3,cf)

            hold on
            fldID=sprintf('pId%d',ii);
%                 lh = plot(thisvar.(fldID).(cfit).c,'bx');
            thisvar.(fldID).(cfit).CMSc = thisvar.(fldID).(cfit).c * exp(i*hh.R) * hh.S + hh.T;
            lh = plot(thisvar.(fldID).(cfit).CMSc,'x','MarkerEdgeColor',dataColor);
            cmplx(@text,thisvar.(fldID).(cfit).CMSc,num2str(c))
            centersCMS(c,cf) = thisvar.(fldID).(cfit).CMSc;
            axis equal
            cf=cf+1;
        end


        c=c+1
    catch err
        disp(err)
        err.stack.line
    end
        
%     end
end


fidBar = [];
for ii=1:Nfid
    figure(hf2)
    subplot(2,2,ii)
    fidBar(ii) = mean(allFidCMS(:,ii));
    fidSTD(ii) = std(allFidCMS(:,ii));
%     cmplx(@text,fidBar(ii)+2*fidSTD(ii),['sigma=' num2str(fidSTD(ii))])
%     title(['fiducial ' num2str(ii) ' 1sigma=' num2str(fidSTD(ii))])
    title(['fiducial ' num2str(ii)])
    ch1 = cmplx(@circle,fidBar(ii),3*fidSTD(ii),'Color',dataColor);
    nm4lgnd = char(regexp(strrep(dataTitle,'_','-'),'.*(?=\.mat)','match'));
    hl = appendLegend(ch1,strcat(nm4lgnd,' 3\sigma circle (\sigma=',num2str(fidSTD(ii)),')'));
    
end

Nact = length(thisvar.activeCobras);
xwidth = zeros(1,Nact);
ywidth = zeros(1,Nact);
figure(hf1)
for ii=1:length(data.activeCobras)
    subplot(3,3,ii)
    cBar(ii) = mean(centersCMS(:,ii));
    cSTD(ii) = std(centersCMS(:,ii));
    pIdstr = ['pid' num2str(data.activeCobras(ii))];
    title(pIdstr)
    ch2 = cmplx(@circle,cBar(ii),3*cSTD(ii),'Color',dataColor);
    nm4lgnd = char(regexp(strrep(dataTitle,'_','-'),'.*(?=\.mat)','match'));
    hl = appendLegend(ch2,strcat(nm4lgnd,' 3\sigma circle (\sigma=',num2str(cSTD(ii)),')'));
    
    % Get axis limits and scale
    v = axis;
    xwidth(ii) = v(2)-v(1);
    ywidth(ii) = v(4)-v(3);

%     legend(ch2,strcat('3\sigma circle (\sigma=',num2str(cSTD(ii)),')'))
%     mtit('Red circle is 3sigma')
end

figure(hf1)
for ii = 1:Nact
    subplot(3,3,ii)
    v = axis;
    xW = max(xwidth);
    yW = max(ywidth);
    xM = v(1)+(v(2)-v(1))/2;
    yM = v(3)+(v(4)-v(3))/2;
    xlim([xM-xW/2,xM+xW/2])
    ylim([yM-yW/2,yM+yW/2])
end

output.dnames = dnames;
output.cBar = cBar;
output.cSTD = cSTD;
output.fidBar = fidBar;
output.fidSTD = fidSTD;
output.activeCobras = data.activeCobras;

% RSLT(1,:) = {'Cobra','MeanCenter','Sigma (px)'};
% 
% 
% RSLT{2,1} = 'Fiducial 1.3';
% RSLT{3,1} = 'Fiducial 1.6';
% RSLT{4,1} = 'Fiducial 1.7';
% RSLT{5,1} = 'GEN3 2.5';
% RSLT{6,1} = 'GEN2 2.1';
% RSLT{2,3} = std(allFidCMS(:,1));
% RSLT{3,3} = std(allFidCMS(:,2));
% RSLT{4,3} = std(allFidCMS(:,3));
% RSLT{5,3} = std(centersCMS(:,2));
% RSLT{6,3} = std(centersCMS(:,1));
% RSLT{2,2} = mean(allFidCMS(:,1));
% RSLT{3,2} = mean(allFidCMS(:,2));
% RSLT{4,2} = mean(allFidCMS(:,3));
% RSLT{5,2} = mean(centersCMS(:,2));
% RSLT{6,2} = mean(centersCMS(:,1));
