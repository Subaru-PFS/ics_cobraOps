function output = S2CentersAnalysis(imgsFile, dataTitle,CobraConfig)

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
A = regexp(cellstr([S.name]),'centerms2fw_ImageId_\d*','match');
fidPosRaw = [];
centerPosCCD = [];
centersCMS = [];
Link1s = [];
c=1;
dnames = {};
refpos = CobraConfig.refpos;
allFidCMS = [];

hf1 = existingFigure('Center Spreads');
clf;
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
            for ii=thisvar.activeCobras
                fldID=sprintf('pId%d',ii);
                if length(thisvar.(fldID).CCDpos) == 1
                    stop = 1;
                end
            end
            if stop == 1
                continue;
            end
            
            thisvar = getS2Center(thisvar,true);

            Nfid = length(thisvar.fixedCobras);
            activeCobras = thisvar.activeCobras;

%             disp(thisvar.fixedCobras)
            for ii=1:length(thisvar.fixedCobras)
                fldID=sprintf('pId%d',thisvar.fixedCobras(ii));
                fidPosRaw(c,ii) = thisvar.(fldID).CCDpos(1);
            end
            
            %Calculate horn transform
            hh = horn87(fidPosRaw(c,:), refpos);
            
            % Apply horns to fiducials
            for ii=1:length(thisvar.fixedCobras)
                figure(hf2)
                subplot(2,2,ii)
                fldID=sprintf('pId%d',thisvar.fixedCobras(ii));
                thisvar.(fldID).CMSpos = thisvar.(fldID).CCDpos * exp(i*hh.R) * hh.S + hh.T;
                plot(thisvar.(fldID).CMSpos,'x','MarkerEdgeColor',dataColor);
                hold on
                cmplx(@text,thisvar.(fldID).CMSpos,num2str(c));
                allFidCMS(c,ii) = thisvar.(fldID).CMSpos;
                axis equal
            end
            
            % Apply horns to centers
            cf=1;
            CW3 = colormap(hsv(length(thisvar.activeCobras)));
            for ii=thisvar.activeCobras
                figure(hf1)
                subplot(3,3,cf)
                hold on
                fldID=sprintf('pId%d',ii);
                tcenter = CobraConfig.ARM_DATA.(['ARM_DATA_' num2str(ii)]).KINEMATICS.Global_base_pos_x + ...
                    CobraConfig.ARM_DATA.(['ARM_DATA_' num2str(ii)]).KINEMATICS.Global_base_pos_y*i;
                thisvar.(fldID).(cfit).CMSc = thisvar.(fldID).(cfit).c * exp(i*hh.R) * hh.S + hh.T;
                lh = plot(thisvar.(fldID).(cfit).CMSc,'x','MarkerEdgeColor',dataColor);
                centersCMS(c,cf) = thisvar.(fldID).(cfit).CMSc;
                Link2s(c,cf) = thisvar.(fldID).(cfit).R;
                Link1s(c,cf) = abs(thisvar.(fldID).(cfit).CMSc - tcenter);
                hold on
                axis equal
                cf=cf+1;
            end
            

            c=c+1;
        catch err
            disp(err)
            err.stack.line
        end
        
    
end


fidBar = [];
for ii=1:Nfid
    figure(hf2)
    subplot(2,2,ii)
    fidBar(ii) = mean(allFidCMS(:,ii));
    fidSTD(ii) = std(allFidCMS(:,ii));
    title(['fiducial ' num2str(ii)])
    ch1 = cmplx(@circle,fidBar(ii),3*fidSTD(ii),'Color',dataColor);
    appendLegend(ch1,strcat(dataTitle,' 3\sigma circle (\sigma=',num2str(fidSTD(ii)),')'))
end

Nact = length(activeCobras);
xwidth = zeros(1,Nact);
ywidth = zeros(1,Nact);
figure(hf1)
for ii=1:Nact
    subplot(3,3,ii)
    cBar(ii) = mean(centersCMS(:,ii));
    cSTD(ii) = std(centersCMS(:,ii));
    pIdstr = ['pid' num2str(activeCobras(ii))];
    title(pIdstr)
    ch2 = cmplx(@circle,cBar(ii),3*cSTD(ii),'Color',dataColor);
    nm4lgnd = char(regexp(strrep(dataTitle,'_','-'),'.*(?=\.mat)','match'));
    hl = appendLegend(ch2,strcat(nm4lgnd,' 3\sigma circle (\sigma=',num2str(cSTD(ii)),')'));
    
    % Get axis limits and scale
    v = axis;
    xwidth(ii) = v(2)-v(1);
    ywidth(ii) = v(4)-v(3);
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
output.Link1s = Link1s;
output.Link2s = Link2s;
output.activeCobras = activeCobras;



%% Plot Link Lengths
hf3 = existingFigure('Link Lengths');
for ii=1:Nact

    figure(hf3)
    % Link 1
    subplot(9,2,2*(ii-1)+1)
    plot(Link1s(:,ii),'Color',dataColor)
    hold on
    ylabel('px')
    % Find std and mean
    stdStr = sprintf('pId%d_stdL1',activeCobras(ii));
    assignin('caller',stdStr,std(Link1s(:,ii)))
    meanStr = sprintf('pId%d_meanL1',activeCobras(ii));
    assignin('caller',meanStr,mean(Link1s(:,ii)))
    % Plot mean line
    mh1 = plot(evalin('caller',meanStr)*ones(length(Link1s(:,ii))),'--','Color',dataColor);
    title(['pId' num2str(activeCobras(ii)) ' \sigma=' num2str(evalin('caller',stdStr))]);
    nm4lgnd = char(regexp(strrep(dataTitle,'_','-'),'.*(?=\.mat)','match'));
    hl = appendLegend(mh1,strcat(nm4lgnd,' Mean=',num2str(evalin('caller',meanStr))));
    
    % Link 2
    subplot(9,2,2*ii)
    plot(Link2s(:,ii),'Color',dataColor)
    hold on
    ylabel('px')
    % Find std and mean
    stdStr = sprintf('pId%d_stdL2',activeCobras(ii));
    assignin('caller',stdStr,std(Link2s(:,ii)))
    meanStr = sprintf('pId%d_meanL2',activeCobras(ii));
    assignin('caller',meanStr,mean(Link2s(:,ii)))
    % Plot mean line
    mh2 = plot(evalin('caller',meanStr)*ones(length(Link2s(:,ii))),'--','Color',dataColor);
    % Title and legend
    title(['pId' num2str(activeCobras(ii)) ' \sigma=' num2str(evalin('caller',stdStr))]);
    nm4lgnd = char(regexp(strrep(dataTitle,'_','-'),'.*(?=\.mat)','match'));
    hl = appendLegend(mh2,strcat(nm4lgnd,' Mean=',num2str(evalin('caller',meanStr))));
%     legend('Link2',strcat('Mean=',num2str(evalin('caller',meanStr))))
%     figure(500)
%     subplot(2,3,ii)
%     XL = xlim;
%     dummyx = linspace(XL(1),XL(2),10);
%     HH = plot(dummyx,evalin('caller',meanStr)*ones(length(dummyx)),'--','Color',dataColor);
% %     legend(HH(1),strcat('Mean Radius=',num2str(evalin('caller',meanStr))))
%         nm4lgnd = char(regexp(strrep(dataTitle,'_','-'),'.*(?=\.mat)','match'));
%     hl = appendLegend(mh2,strcat(nm4lgnd,' Mean=',num2str(evalin('caller',meanStr))));
end

