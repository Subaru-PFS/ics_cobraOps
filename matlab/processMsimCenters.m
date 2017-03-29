function output = processMsimCenters(baseCfg);
% processMsimMaps
 
% Get theta dots.
thetaFw = makeMsimCenters(1,1,'thetaFwMap');

% Get phi dots.
phiFw = makeMsimCenters(2,1,'PhiFwMap');

% Get phi dots.
try 
cwthtmap = makeMsimCenters(2,1,'CWThtMap');
catch
end

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
       %sort out radial outliers.
       pok = pok & phiRadii(1:end-1) < (mean(phiRadii) + std(phiRadii));
       pok = pok & phiRadii(1:end-1) > (mean(phiRadii) - std(phiRadii)); 
       
        cfit1 = circfit(phiCentroids(pok));
       phiCentroidsCM = phiCentroids - cfit1.c;
       phiRadii = abs(phiCentroidsCM);
       j2center    = cfit1.c;
   else % if finding a center failed due to no valid points, create fake center. 
       j2center    = 0 + 1i;
   end
   
   j2centroids = phiCentroids;
   j2ok        = pok;
      %% PHI, ie ,J2 on the pos hardstop:
     try
   phiIndex2 = find([cwthtmap.pId] == pid);
   phiCentroids2 =  cwthtmap(phiIndex2).rawData(:,1:2) * [1;1i];
   cfit0 = circfit(phiCentroids2);
    
   phiCentroids2CM = phiCentroids2 - cfit0.c; 
   phiAngles2 = angle(phiCentroids2CM);
   phiUangles2 = unwrap(phiAngles2); % use a continuous set of angles for calcs
   
   % filter out hard stops.  need to deal with the fact that diff
   % reduces the length of the vector by 1.
   pok2 = (abs(diff(phiUangles2)) > .01);
   pok2 = pok2 & [0; pok2(1:end-1)];
   % refit the data  
   if(sum(pok2)>3)
       cfit2 = circfit(phiCentroids2(pok2));
       phiCentroids2CM = phiCentroids2 - cfit2.c;
       phiRadii = abs(phiCentroids2CM);
       j2center2    = cfit2.c;
   else % if finding a center failed due to no valid points, create fake center. 
       j2center2    = 0 + 1i;
   end

   j2centroids2 = phiCentroids2;
   j2ok2        = pok2;
     catch
         % this is the case for old folders
         % create empty lists if required.
         j2centroids2 =[];
         j2ok2 = [];
              j2center2    = 0 + 1i;

     end
     
   output(ctr) = packstruct(pid,j1centroids,j1center,j1ok,...
                            j2centroids,j2center,j2ok,j2centroids2, j2center2, j2ok2);
   
   ctr = ctr + 1; % increment counter
end
