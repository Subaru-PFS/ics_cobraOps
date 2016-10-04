%TARGET GENERATOR
clear all
close all
addpath 'C:/Users/sage/Desktop/Dropbox/PFS_EM/SVN/MATLAB/collisions/'

%% INPUTS
% Get geometry from xml
empt = [];
bench = defineBenchGeometry(empt, 'true', 'true', 1);

%activeCobras = [1,5,9,13,17,21,25,29,33,37,41,45,49,53]; % Same purpose as movingCobras but more intuitive to use pId
%activeCobras = [9,15,31,37,45]; % Same purpose as movingCobras but more intuitive to use pId
%activeCobras = [3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51];
%activeCobras = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53];
numPos = length(bench.center);
% Specify number of random targets to generate
numtrg = 200; 
minDist = 2000/90; % 2mm in pixel
LF = [10*ones(1,numtrg-1) 48]; % quick hack to remove \n after last target.



%% EXECUTION

 
% for ii = 1:length(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER);
%     cnfgID = sprintf('ARM_DATA_%d',ii);
%     
%   
%     center = [center; getARMval(CobraConfig,ii, moduleID,'Global_base_pos_x') + ...
%         getARMval(CobraConfig,ii, moduleID, 'Global_base_pos_y')*i];
%     L1 = [L1 getARMval(CobraConfig,ii, moduleID,'Link1_Link_Length')];
%     L2 = [L2 getARMval(CobraConfig,ii, moduleID,'Link2_Link_Length')];
%     if(~isempty(getARMval(CobraConfig,ii, moduleID,'phiMax')))
%         phiMax = [phiMax getARMval(CobraConfig,ii, moduleID,'phiMax')];
%     else
%         phiMax = [phiMax 160];
%     end
%     if(~isempty(getARMval(CobraConfig,ii, moduleID,'phiMin')))
%         phiMin = [phiMin getARMval(CobraConfig,ii, moduleID,'phiMin')];
%     else
%         phiMin = [phiMin 5];
%     end 
% end

  
 
KeepOutAngle = 0.1;
%% Rule #0 stick to the geometry (inner and outer radii)
% Inner keepout radius
Rmin = abs(bench.L1 + bench.L2 .* exp(1i.*bench.phiIn));

% Outer reach
Rmax = abs(bench.L1 + bench.L2 .* exp(1i.*bench.phiOut));

% angular coordinate of target
THT = rand(numPos,numtrg)*2*pi;

% radial coordinate of target
% radial coordinate of target
dA = 1 ./ ( (Rmax./Rmin).^2 - 1 ); %fraction of the keepout area.
RDS = sqrt(bsxfun(@times, (bsxfun(@plus, rand(size(THT)), dA)), (Rmax.^2 - Rmin.^2)));


% RDS = -ones(numPos,numtrg);
% for ii = 1:numPos
%   while(true)
%     % "indices of radii to assign" 
%     jj_R_assign = find(RDS(ii,:) < Rmin(ii));
%     n_assign = length(jj_R_assign);
%     % asssign radii, choosing uniformly over a disk of radius Rmax
%     % or end if there's nothing to do.
%     if n_assign > 0
%       RDS(ii,jj_R_assign) = Rmax(ii) * sqrt(rand(1,n_assign));
%     else
%       break;
%     end
%   end
% end
 
%% Rule 1 No targets closer than 2mm to any line in its target location. 
disp('Number of bad targets from rule #2 (< 2mm distance)')
targets = bsxfun(@plus, RDS.*exp(1i*THT), bench.center);

targetsTP  = XY2TP(bsxfun(@minus, targets, bench.center), bench.L1, bench.L2); 
numcoll = 0;
%numcoll = zeros(1,numtrg);
% 
% figure(1)
% plot(targets,'b.');
% axis equal;
% hold on;
%  
% for kk = 1:numtrg
%     firstRun = true;
%     while(true)
%         distances = CalcDistanceMatrix(targets(:,kk), bench);
%         distances = symmetrize(distances); 
%         if(sum(distances.dst<bench.minDist) == 0) 
%             break;
%         else
%             if(firstRun)
%                 firstRun = false;  
%                 numcoll = numcoll + sum(sum(distances.dmatrix<bench.minDist));  
%                 %numcoll(kk) = sum(sum(distances.dmatrix < minDist));
%                 [row col] = find(distances.dmatrix<minDist);
%                 plot(targets(row,kk),'r.','markersize', 15);
%             end
%             mcc = find((sum(distances.dmatrix<minDist,2) == max(sum(distances.dmatrix<bench.minDist,2)))); % returns the one who collides with most of the others.
%             
%             RDSnew = sqrt( ( rand(1) + dA(mcc(1)) ) * (Rmax(mcc(1))^2 - Rmin(mcc(1))^2) );
%             THTnew = rand(1) * 2 * pi;
%             targets(mcc(1),kk) =  bench.center(mcc(1)) + RDSnew * exp(1i * THTnew);       
%         end
%     end
% end
% $$$ for kk = 1:numtrg
% $$$   distances = CalcDistanceMatrix(targets(:,kk), center, L1, L2);
% $$$   distances = symmetrize(distances); 
% $$$   numcoll(kk) = sum(sum(distances.dmatrix < minDist));
% $$$   [row col] = find(distances.dmatrix < minDist);
% $$$   plot(targets(row,kk),'r.','markersize', 15);
% $$$   if numcoll(kk) > 10, keyboard; end;
% $$$ end
% collisions = numcoll/2;
% %collisions = sum(numcoll)/2;
% percentColl = 100* collisions/(numtrg * numPos)
% 
% cmplx(@plotcircle,bench.center,Rmax,'k');


  jj= 0;
  
fh1 = figure;
hold on;
 for ii = 1:length(bench.pids)
  jj = jj +1;
  tfile = fopen(sprintf('TargetList_mId_%d_pId_%d.txt',bench.mids(ii),bench.pids(ii)), 0,0, 'w');
  plot(targets(jj,:),'r.');
  % label targets with integers
  %cmplx(@text,targets(ii,:),strsplit(num2str(1:length(targets(ii,:))),' '));
  %drawnow;
  targetoutput = [simple(targets(jj,:).')' ; LF];
  fprintf(tfile,'%4.2f,%4.2f%c',targetoutput);
  fclose(tfile);
end

axis equal
hold off;
