function output = processMsimCenters(baseCfg);
% processMsimMaps
 
% Get theta dots.
thetaFw = makeMsimCenters(1,1,'thetaFwMap');

% Get phi dots.
phiFw = makeMsimCenters(2,1,'PhiFwMap');


data = [];
ctr=1;
for pid = [thetaFw.pId]
    pidname = sprintf('pid%02d', pid);

   %% THETA , ie, J1
   thetaCentroids = thetaFw(ctr).rawData(:,1:2) * [1;1i];
   cfit0     = circfit(thetaCentroids);
    
    
   thetaCentroidsCM = thetaCentroids - cfit0.c; 
   thetaAngles = angle(thetaCentroidsCM);
   thetaUangles = unwrap(thetaAngles); % use a continuous set of angles for calcs
   
   % filter out hard stops.  need to deal with the fact that diff
   % reduces the length of the vector by 1. 
    tok = (abs(diff(thetaUangles)) > .01);
   %   ok = ok & [0; ok(1:end-1)];
   tok = [tok ; 1] & [0 ; tok]; % always toss the first data point.
   
   % refit the data    
   cfit1 = circfit(thetaCentroids(tok));
   thetaCentroidsCM = thetaCentroids - cfit1.c;
   thetaRadii = abs(thetaCentroidsCM);
 
   % collect outputs
   j1centroids = thetaCentroids;
   j1center    = cfit1.c;
   j1ok        = tok;
   
   %% PHI, ie ,J2
   phiIndex = find([phiFw.pId] == pid);
   phiCentroids =  phiFw(phiIndex).rawData(:,1:2) * [1;1i];
   cfit0 = circfit(phiCentroids);
    
   phiCentroidsCM = phiCentroids - cfit0.c; 
   phiAngles = angle(phiCentroidsCM);
   phiUangles = unwrap(phiAngles); % use a continuous set of angles for calcs
   
   % filter out hard stops.  need to deal with the fact that diff
   % reduces the length of the vector by 1.
   pok = (abs(diff(phiUangles)) > .01);
   pok = pok & [0; pok(1:end-1)];
   % refit the data  
   if(sum(pok)>3)
       cfit1 = circfit(phiCentroids(pok));
       phiCentroidsCM = phiCentroids - cfit1.c;
       phiRadii = abs(phiCentroidsCM);
       j2center    = cfit1.c;
   else % if finding a center failed due to no valid points, create fake center. 
       j2center    = 0 + 1i;
   end

   j2centroids = phiCentroids;
   j2ok        = pok;
   
   output(ctr) = packstruct(pid,j1centroids,j1center,j1ok,...
                            j2centroids,j2center,j2ok);
   
   
   ctr = ctr + 1; % increment counter
end
