%% FIDUCIAL MEASUREMENTS AND REFERENCE POSITION SETUP
% This script will process several msim centroid logs containing fiducial
% centroids and find the average and standard error. It will write these to
% a Cobra Config XML
clear all
close all

%% INPUTS

% Provide paths to centroid logs with good fiducial centroids
logFiles = {
    'D:\PfsTests\02_03_15_14_04_00\Log\Centroid_log.txt'
    };
    
%     'D:\PfsTests\06_24_14_17_56_21_TargetRun\Log\Centroid_log.txt'};

% Provide the starting config xml (usually at this step you should choose
% the template xml)
baseCfg = loadCfgXml; % Uncomment this to have dialog popup
% baseCfg = loadCfgXml('.','CobraConfigTemplate.xml'); 


%% EXECUTION (BEWARE OF MESSING WITH WHATS BELOW)
doHorns = true;

for ii=1:length(logFiles)
    logFile = logFiles{ii};
    data = processCentroidLog(logFile);
    varName = strcat('FID_',char(regexp(logFile,'(\d\d_)*','match')));
    assignin('base',varName,data)
    Nfid = data.nFid;
end




% load('C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\FiducialMeasurements\cntrds_03_31_14_15_16_44.mat');
% tdata = cntrds_03_31_14_15_16_44;
% tdata.fixedCobras = [6,7,8];
% tdata.activeCobras = [];
% for ii=tdata.fixedCobras
%     fldID = ['pId',num2str(ii)];
%     tdata.(fldID).CCDpos = tdata.(fldID).CCDpos + (-1-1*i);
% end
% assignin('base','FID_03_31_14_15_16_44_ML',tdata)



clear data
S = whos;
filterSTR = 'FID_';
tdata.name = 'FF_Measurements';

figure(1)
hold on
CW = colormap(lines(length(logFiles)));
c=1;
legendStr = {};
plotHandles = [];
refall = [];

for vv = 1:length(S)
    thisvar = eval(S(vv).name);
    if regexp(S(vv).name,filterSTR)
        thisvar.CCDrefpos = [];
        nFid = thisvar.nFid
        
        for ii=1:nFid
           
            fldID = sprintf('fId%d',ii);
            
            figure(1)
            subplot(floor(nFid/2),ceil(nFid/2),ii)
            hold on
            title(sprintf('Fiducial #%d',ii))
            
            
            if isfield(thisvar.(fldID),'RawPos')
                posStr = 'RawPos';
            elseif isfield(thisvar.(fldID),'CCDpos')
                posStr = 'CCDpos';
            elseif isfield(thisvar.(fldID),'CLpos')
                posStr = 'CLpos';
            end
            
%             %%%%%%TESTING%%%%%%%%%%%%%%%%%%%
%             if regexp(S(vv).name,'MSIM') & thisvar.fixedCobras(ii)==7
%                 lh = plot(thisvar.(fldID).(posStr),'x','MarkerEdgeColor',CW(c,:));
%                 figure(10)
%                 IRght = find(real(thisvar.(fldID).(posStr))>1056.7)
%                 ILft = find(real(thisvar.(fldID).(posStr))<=1056.7)
%                 lh = plot(thisvar.(fldID).(posStr),'x','MarkerEdgeColor',CW(c,:));
%                 hold on
%                 axis equal
%                 figure(1)
%             elseif regexp(S(vv).name,'MATLAB') & ii==thisvar.fixedCobras(2)
%                 lh = plot(thisvar.(fldID).(posStr),'x','MarkerEdgeColor',CW(c,:));
%                 figure(10)
%                 lh = plot(thisvar.(fldID).(posStr)(IRght),'ro');
%                 lh = plot(thisvar.(fldID).(posStr)(ILft),'ko');
%                 hold on
%                 axis equal
%                 figure(1)
%             else
                lh = plot(thisvar.(fldID).(posStr),'o','MarkerEdgeColor',CW(c,:));
%                 for kk = 1:length(thisvar.(fldID).(posStr))
%                     cmplx(@text,thisvar.(fldID).(posStr)(kk),num2str(kk));
%                 end
%             end
            %%%%%%TESTING%%%%%%%%%%%%%%%%%%%
            
            
               
            
            % Filter outliers that are greater than 5 sigma from mean
%             pos = thisvar.(fldID).(posStr);
%             thisvar.(fldID).(posStr) = pos(abs(pos-Pbar) < 7*Pstd)
            
            Pbar = mean(thisvar.(fldID).(posStr));
            Pstd = std(thisvar.(fldID).(posStr))/sqrt(2);
            
            
            axis equal
            
%             if ~isfield(data,fldID)
%                 data.(fldID).CCDpos = [];
%             end
%                 
%             data.(fldID).CCDpos = [data.(fldID).CCDpos thisvar.(fldID).(posStr)];
            
            thisvar.CCDrefpos = horzcat(thisvar.CCDrefpos, thisvar.(fldID).(posStr).');
            

    %         cmplx(@circle,Pbar,Pstd)
            hc1 = circle(real(Pbar),imag(Pbar),3*Pstd,'Color',CW(c,:));
            hl1 = appendLegend(hc1,'3\sigma Circle');
%             disp(S(vv).name)
            C = regexp(S(vv).name,'_','split');
            nameF = sprintf('%s ',C{:});
            hl1 = appendLegend(lh,strcat(nameF,' - \sigma=',num2str(Pstd)));
%             set(hl1,'Interpreter','latex')
%             cmplx(@text,Pbar,strcat('\sigma=',num2str(Pstd)))
                            
                    
        end
        
%         disp(thisvar.CCDrefpos)
        refall = vertcat(refall,thisvar.CCDrefpos);

        legendStr{c} = S(vv).name;
        plotHandles(c) = lh;
        c=c+1;
        
        assignin('base',S(vv).name,thisvar)
    end
end
for ii=1:Nfid
    subplot(2,3,ii)
        Pbar = mean(refall(:,ii));
        Pstd = std(refall(:,ii))/sqrt(2);
        hc1 = circle(real(Pbar),imag(Pbar),3*Pstd,'k');
        appendLegend(hc1,strcat('Overall \sigma=',num2str(Pstd)));
end


% legend(plotHandles,legendStr,'Interpreter','none')
% mtit('Raw Fiducial Positions')


if doHorns
    %% Mean Fiducial positions
%     CobraConfigFile = '../Metrology/config_01-17-14.mat';
%     load(CobraConfigFile)
%     refpos = CobraConfig.refpos;
    % Calculate mean position of all fiducials wrt CCD coord
    refpos = mean(refall);

    %% Establish Origin

    % % This block is to be used if fiducial avg is origin
    % % Calculate new origin defined by mean of fiducial positions
    % origin = mean(refpos);
    % % Calculate new reference positions of fiducials wrt common origin
    % refpos = refpos - origin;

    % Use this block for CCD as origin
    origin = 0;


    figure(2)
    set(2,'Name','Horns Adjusted Fiducials')
    hold on
    c=1;
    FFrefall = [];
    Ntot = 0;
    for vv = 1:length(S)
        thisvar = eval(S(vv).name);
        if regexp(S(vv).name,filterSTR)

            disp(S(vv).name)
            % For each image
            for ii=1:length(thisvar.CCDrefpos)

                % Derive horn's transfomation
                hh = horn87(thisvar.CCDrefpos(ii,:),refpos);

                % Apply horn's transformation to each Cobra
                for jj = 1:thisvar.nFid
                    
                    % Get field id
                    fldID = sprintf('fId%d',jj);

                    % Check if data is from matlab or msim log
                    if isfield(thisvar.(fldID),'CCDpos')
                        posStr = 'CCDpos';
                    elseif isfield(thisvar.(fldID),'CLpos')
                        thisvar.(fldID).FFpos(ii) = thisvar.(fldID).CLpos(ii);
                    end

                    % Apply horn
                    thisvar.(fldID).FFpos(ii) =  thisvar.(fldID).(posStr)(ii)* exp(i*hh.R) * hh.S + hh.T;
                end
            end

            % For each fiducial
            temp = [];
            for ii=1:thisvar.nFid
                    % Get field id
                    fldID = sprintf('fId%d',ii);

                    % Plot fiducial positions
                    figure(2)
                    subplot(2,3,ii)
                    hold on
                    title(sprintf('Fiducial #%d',ii))
%                     lh = plot(thisvar.(fldID).FFpos,'o','MarkerEdgeColor',CW(c,:));
lh = plot(thisvar.(fldID).FFpos,'o','color',CW(c,:));
                    FFstd = std(thisvar.(fldID).FFpos)/sqrt(2);
                    C = regexp(S(vv).name,'_','split');
                    nameF = sprintf('%s ',C{:});
                    appendLegend(lh,strcat(nameF,' - \sigma=',num2str(FFstd)))
                    if ii==1
                        Ntot = Ntot + length(thisvar.(fldID).FFpos);
                    end
                    plot(refpos(ii),'x')
                    axis equal
                    temp = horzcat(temp,thisvar.(fldID).FFpos.');
            end

            thisvar.allFids = temp;
            FFrefall = vertcat(FFrefall,temp);
            c=c+1;

            assignin('base',S(vv).name,thisvar)
        end

    end

%     Rows = floor(Nfid/2)+1;
%     Cols = ceil(Nfid/2);
    PosM = {};
    for ii=1:Nfid
        figure(2)
        subplot(2,3,ii)
        hold on
        fldID = sprintf('fId%d',ii);
        disp(fldID)
        Pbar = mean(FFrefall(:,ii));
        Pstd = std(FFrefall(:,ii))/sqrt(2);
        hc1 = circle(real(Pbar),imag(Pbar),3*Pstd,'k');
        appendLegend(hc1,strcat('Overall \sigma=',num2str(Pstd)));
        refpos(ii) = Pbar;
        PosM{1,1} = 'Fiducial';
        PosM{2,1} = 'Position';
        PosM{3,1} = 'Standard Deviation';
        PosM{1,ii+1} = fldID;
        PosM{2,ii+1} = Pbar;
        PosM{3,ii+1} = Pstd;
    end
%     mtit('Horn-adjusted Positions')

newCfg = baseCfg;
newCfg.ARM_DATA.refpos = refpos;

[xmlfile, xmlfilepath] = uiputfile('*.xml','Save new CobraConfig XML file with reference positions');
cobraCfg2xml(newCfg,fullfile(xmlfilepath,xmlfile));  
end
    

