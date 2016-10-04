function [output, CobraConfig] = loadTestData(dataDirs, CobraConfig, useLocal)

% Specify prefix of files to screen for.
loadFilter = '5um';
 
isLoaded = 0;

if(isstruct(dataDirs))
   nTests = 1;
   isLoaded = 1;
else
   nTests = length(dataDirs);
end
if exist('useLocal') 
if(useLocal)
nTests = 1;
end
else
   useLocal = 0; 
end


% Load cobra config xml
if exist('CobraConfig') 
   
else
    %cfgFile = dir2cell([dataDirs{end} '/*.xml']);
    CobraConfig = loadCfgXml();
end


numPos = 57;
pos = [];
s1 = [];
s2 = [];

for id=1:numPos % 9 being the number of positioners used.
    pos(id).s1.steps = [];
    pos(id).s1.moveSizes = [];
    pos(id).s1.startAngle = [];
    pos(id).s1.finishAngle = [];
    pos(id).s1.targetAngle = [];
    pos(id).s1.mmap = [];
    pos(id).s1.iteration = [];
    pos(id).s1.targetNo = [];
    pos(id).xst = [];% x location of start of movement
    pos(id).yst = [];
    pos(id).xtg = [];% x location of target
    pos(id).ytg = [];
    pos(id).xfn = [];% x location of end of movement
    pos(id).yfn = [];
    pos(id).hs = []; % 1 = cw hardstop; 2 = ccw hardstop
    
   
    
    pos(id).s2.steps = [];
    pos(id).s2.moveSizes = [];
    pos(id).s2.startAngle = [];
    pos(id).s2.finishAngle = [];
    pos(id).s2.targetAngle = [];
    pos(id).s2.mmap = [];
    pos(id).s2.iteration = [];
    pos(id).s2.targetNo = [];
end
  
for ii=1:nTests;
    if(~isLoaded)
        if(useLocal)
            QQ = loadmats([loadFilter '_mId*.mat']);
        else
            dataDir = dataDirs{ii};
            QQ = loadmats([loadFilter '_mId*.mat'],dataDir);
        end
    eval(unpackstruct(QQ));
    S = whos;
    Names = {S.name};
    posFileInd = find(~cellfun(@isempty,regexp(Names,'^mId')));
    end
    for vv = posFileInd % for all positioners
        if(isLoaded)
            vname = '_str';
            keyboard;
        else 
        vname='';
        vname = char(Names(vv));
        id = str2num(char(regexp(vname,'(?<=pId_)\d*','match')));
        pos(id).name = vname;
        s = eval(vname);
        end
        if ~isempty(regexp(vname,'_str'))
             
            angles = [s.J1].*180/pi;
            targetAngles= [s.J1_t].*180/pi;
          
            tt  = [s.J1_S];
            xp = real([s.curPos]);
            yp = imag([s.curPos]);
            xt = real([s.targCmplx]);
            yt = imag([s.targCmplx]);
            ms = diff(angles);
            iter = [s.iter];   
            pos(id).s1.moveSizes =     horzcat(pos(id).s1.moveSizes ,ms);
            pos(id).s1.steps =         horzcat(pos(id).s1.steps, tt(1:end-1,:));
            pos(id).s1.startAngle =    horzcat(pos(id).s1.startAngle, angles(1:end-1,:));
            pos(id).s1.finishAngle =   horzcat(pos(id).s1.finishAngle, angles(2:end,:));
            pos(id).s1.mmap =          horzcat(pos(id).s1.mmap, ms ./ tt(1:end-1,:));
            pos(id).s1.iteration =     horzcat(pos(id).s1.iteration, iter(1:end-1,:));
            pos(id).xst =              horzcat(pos(id).xst, xp(1:end-1,:)); % x location of start of movement
            pos(id).yst =              horzcat(pos(id).yst, yp(1:end-1,:));% y location of start of movement
            pos(id).xtg =              horzcat(pos(id).xtg, xt(1:end,:));
            pos(id).ytg =              horzcat(pos(id).ytg, yt(1:end,:));
            pos(id).xfn =              horzcat(pos(id).xfn, xp(2:end,:));
            pos(id).yfn =              horzcat(pos(id).yfn, yp(2:end,:));
            pos(id).hs =               horzcat(pos(id).hs, [s.hardstop]); % Hardstops
            pos(id).s1.targetNo =      horzcat(pos(id).s1.targetNo, meshgrid(1:size(ms,2), 1:size(ms,1)));
            pos(id).s1.targetAngle =   horzcat(pos(id).s1.targetAngle, targetAngles(1:end-1,:));
            
            angles = [s.J2].*180/pi;
            targetAngles= [s.J2_t].*180/pi;
            tt = [s.J2_S];
            ms = diff(angles);
            pos(id).s2.moveSizes =     horzcat(pos(id).s2.moveSizes ,ms);
            pos(id).s2.steps =         horzcat(pos(id).s2.steps, tt(1:end-1,:));
            pos(id).s2.startAngle =    horzcat(pos(id).s2.startAngle, angles(1:end-1,:));
            pos(id).s2.finishAngle =   horzcat(pos(id).s2.finishAngle, angles(2:end,:));
            pos(id).s2.targetAngle =   horzcat(pos(id).s2.targetAngle, targetAngles(1:end-1,:));
            pos(id).s2.mmap =          horzcat(pos(id).s2.mmap, ms ./ tt(1:end-1,:));
            pos(id).s2.iteration =     horzcat(pos(id).s2.iteration, iter(1:end-1,:));
            pos(id).s2.targetNo =      horzcat(pos(id).s2.targetNo, meshgrid(1:size(ms,2), 1:size(ms,1)));
        end
    end
end

for id=1:numPos

  fwdmoves = pos(id).s1.steps > 0;
  fwdmoves = fwdmoves & ~isnan(pos(id).s1.mmap);
  rvsmoves = pos(id).s1.steps < 0;
  rvsmoves = rvsmoves & ~isnan(pos(id).s1.mmap);
  nonmoves = pos(id).s1.steps == 0;
  

  % Filter forward Moves
  pos(id).s1f.moveSizes        = pos(id).s1.moveSizes(fwdmoves);
  pos(id).s1f.steps            = pos(id).s1.steps(fwdmoves);
  pos(id).s1f.startAngle       = pos(id).s1.startAngle(fwdmoves);
  pos(id).s1f.finishAngle      = pos(id).s1.finishAngle(fwdmoves);
  pos(id).s1f.mmap             = pos(id).s1.mmap(fwdmoves);
  pos(id).s1f.iteration        = pos(id).s1.iteration(fwdmoves);
  pos(id).s1f.targetNo        = pos(id).s1.targetNo(fwdmoves);
  
  % Filter reverse Moves       
  pos(id).s1r.moveSizes        = pos(id).s1.moveSizes(  rvsmoves);
  pos(id).s1r.steps            = pos(id).s1.steps(      rvsmoves);
  pos(id).s1r.startAngle       = pos(id).s1.startAngle( rvsmoves);
  pos(id).s1r.finishAngle      = pos(id).s1.finishAngle(rvsmoves);
  pos(id).s1r.mmap             = pos(id).s1.mmap(       rvsmoves);
  pos(id).s1r.iteration        = pos(id).s1.iteration(  rvsmoves);
  pos(id).s1r.targetNo        = pos(id).s1.targetNo(  rvsmoves);
  % Filter non-moves
  pos(id).s1noMove.moveSizes   = pos(id).s1.moveSizes(nonmoves);
  pos(id).s1noMove.steps       = pos(id).s1.steps(nonmoves);
  pos(id).s1noMove.startAngle  = pos(id).s1.startAngle(nonmoves);
  pos(id).s1noMove.finishAngle = pos(id).s1.finishAngle(nonmoves);
  pos(id).s1noMove.mmap        = pos(id).s1.mmap(nonmoves);
  pos(id).s1noMove.iteration   = pos(id).s1.iteration(nonmoves);
  pos(id).s1noMove.targetNo   = pos(id).s1.targetNo(nonmoves);
  
  
  fwdmoves = pos(id).s2.steps > 0;
  fwdmoves = fwdmoves & ~isnan(pos(id).s2.mmap);
  rvsmoves = pos(id).s2.steps < 0;
  rvsmoves = rvsmoves & ~isnan(pos(id).s2.mmap);
  nonmoves = pos(id).s2.steps == 0;


  % Filter forward Moves
  pos(id).s2f.moveSizes        = pos(id).s2.moveSizes(fwdmoves);
  pos(id).s2f.steps            = pos(id).s2.steps(fwdmoves);
  pos(id).s2f.startAngle       = pos(id).s2.startAngle(fwdmoves);
  pos(id).s2f.finishAngle      = pos(id).s2.finishAngle(fwdmoves);
  pos(id).s2f.mmap             = pos(id).s2.mmap(fwdmoves);
  pos(id).s2f.iteration        = pos(id).s2.iteration(fwdmoves);
  pos(id).s2f.targetNo        = pos(id).s2.targetNo(fwdmoves);
  % Filter reverse Moves       
  pos(id).s2r.moveSizes        = pos(id).s2.moveSizes(rvsmoves);
  pos(id).s2r.steps            = pos(id).s2.steps(rvsmoves);
  pos(id).s2r.startAngle       = pos(id).s2.startAngle(rvsmoves);
  pos(id).s2r.finishAngle      = pos(id).s2.finishAngle(rvsmoves);
  pos(id).s2r.mmap             = pos(id).s2.mmap(rvsmoves);
  pos(id).s2r.iteration        = pos(id).s2.iteration(rvsmoves);
  pos(id).s2r.targetNo        = pos(id).s2.targetNo(rvsmoves);
  % Filter non-moves
  pos(id).s2noMove.moveSizes   = pos(id).s2.moveSizes(nonmoves);
  pos(id).s2noMove.steps       = pos(id).s2.steps(nonmoves);
  pos(id).s2noMove.startAngle  = pos(id).s2.startAngle(nonmoves);
  pos(id).s2noMove.finishAngle = pos(id).s2.finishAngle(nonmoves);
  pos(id).s2noMove.mmap        = pos(id).s2.mmap(nonmoves);
  pos(id).s2noMove.iteration   = pos(id).s2.iteration(nonmoves);
  pos(id).s2noMove.targetNo   = pos(id).s2.targetNo(nonmoves);
end

%% Get old Motor Map values too
if(~useLocal)
pos = getOldMotorMapValues(CobraConfig, pos);
end
% file = strsplit(CobraConfig.cfgFile, '\');
% file = file(1);
% pos(1).folderPath = file;
output = pos;
end