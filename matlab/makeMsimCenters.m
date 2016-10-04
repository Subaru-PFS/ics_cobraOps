function logData = makeMsimCenters(stage,direction,logPrefix,Nsteps,moveMin,movePersist)


% Find the positioner centroid logs in current directory and store pids
posFiles = dir2cell([logPrefix '_mId_*pId_*.txt']);

mid_cells = regexp(posFiles,'(?<=mId_)\d*','match');
pid_cells = regexp(posFiles,'(?<=pId_)\d*','match');
for jj = 1:length(posFiles)
  mId(jj)  = str2num(mid_cells{jj}{1}); %% argghh! need to
                                        %% dereference cells!!
  mpId(jj) = str2num(pid_cells{jj}{1});
end

pId = (mId-1)*57 + mpId;

% For each log file
for ii=1:length(posFiles)

    logData(ii).mId = mId(ii);
    logData(ii).pId = pId(ii); 
    logData(ii).rawData = dlmread(posFiles{ii},','); 
     
end  
return

