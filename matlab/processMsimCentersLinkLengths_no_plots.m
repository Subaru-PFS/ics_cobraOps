function output = processMsimCentersLinkLengths_no_plots();
% processMsimMaps
 
baseCfg = loadCfgXml;

%% EXECUTION (Don't edit below here)
% Get theta dots.
thetaFw = makeMsimCenters(1,1,'thetaFwMap');

% Get phi dots.
phiFw = makeMsimCenters(2,1,'PhiFwMap');

ctr=1;
for pid = [thetaFw.pId]
    pidname = sprintf('pid%02d', pid);
    
   for kk = 1:sum(regexpcmp(fieldnames(baseCfg.ARM_DATA),'ARM_DATA_*')) % Count all struct that start with the name ARM_DATA_
        ARM_DATA_n = sprintf('ARM_DATA_%d', kk);
        if(pid ==  str2num(baseCfg.ARM_DATA.(ARM_DATA_n).DATA_HEADER.Positioner_Id.Text))
            armid = ARM_DATA_n;
            break;
        end
   end
  
   thetaCentroids = thetaFw(ctr).rawData(:,1:2) * [1;1i];
   cfit0     = circfit(thetaCentroids);
    
    
   thetaCentroidsCM = thetaCentroids - cfit0.c;
   thetaRadii = abs(thetaCentroidsCM);
   thetaAngles = angle(thetaCentroidsCM);
   thetaUangles = unwrap(thetaAngles); % use a continuous set of angles for calcs
   
   % filter out hard stops.  need to deal with the fact that diff
   % reduces the length of the vector by 1.
   %%%%  "ok" logic example
   % X      : 1 2 3 4 4 5 6 7 8 8 8 9 9
   % dX     :  1 1 1 0 1 1 1 1 0 0 1 0
   % [dX 1] : 1 1 1 0 1 1 1 1 0 0 1 0 1
   % [1 dX] : 1 1 1 1 0 1 1 1 1 0 0 1 0
   % AND    : 1 1 1 0 0 1 1 1 0 0 0 0 0
   
   ok = (abs(diff(thetaUangles)) > .01);
   ok = [ok ; 1] & [0 ; ok]; % always toss the first data point.
   % refit the data    
   cfit1 = circfit(thetaCentroids(ok));
   thetaCentroidsCM = thetaCentroids - cfit1.c;
   thetaRadii = abs(thetaCentroidsCM);
   thetaAngles = angle(thetaCentroidsCM);

   %% collect outputs
   theta(ctr).centroids = thetaCentroids;
   theta(ctr).center    = cfit1.c;
   theta(ctr).ok        = ok;
   theta(ctr).pid       = pid;
 
   phiIndex = find([phiFw.pId] == pid);
   phiCentroids =  phiFw(phiIndex).rawData(:,1:2) * [1;1i];
   cfit0 = circfit(phiCentroids);
   
 
   phiCentroidsCM = phiCentroids - cfit0.c;
   phiRadii = abs(phiCentroidsCM);
   phiAngles = angle(phiCentroidsCM);
   phiUangles = unwrap(phiAngles); % use a continuous set of angles for calcs
   
   % filter out hard stops.  need to deal with the fact that diff
   % reduces the length of the vector by 1.
   ok = (abs(diff(phiUangles)) > .01);
   ok = [ok ; 1] & [0 ; ok];
   % refit the data    
   cfit1 = circfit(phiCentroids(ok));
   phiCentroidsCM = phiCentroids - cfit1.c;
   phiRadii = abs(phiCentroidsCM);
   phiAngles = angle(phiCentroidsCM);

 
   phi(ctr).centroids = phiCentroids;
   phi(ctr).center    = cfit1.c;
   phi(ctr).ok        = ok;
   phi(ctr).pid       = pid;
   
% $$$     keyboard;
   
   
   %% Now the calculate the link lengths:
 
   ctr = ctr + 1; % increment counter
end

%% I want jjsrt so that the data structures will be in increasing
%% pid order
[xxsrt jjsrt] = sort([theta.pid]);

output.tht = theta(jjsrt);
output.phi = phi(jjsrt);

