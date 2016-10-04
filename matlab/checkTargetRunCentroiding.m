I=fitsread('D:\PfsTests\03_13_15_13_14_55_msimMaps\Images\thetaFwMapImageId_1_Target_0_Iteration_-1_loopId_-1_64.fits');
figure(123)
imagesc(I)
hold on
% the array indecies from pId array.
data.activeCobras = [1,5,9,13,17,21,25,29,33,37,41,45,49,53]; % Same purpose as movingCobras but more intuitive to use pId

% Positioners which should be treated as fiducials. @Johannes change for
% debug, is set to the max number of fiducials found. These are acutally
% the real fiducals not the cobras.
data.fixedCobras = [3,7,11,15,19,23,27,31,35,39,43,47,51,]; % Same purpose as Fiducials, but more intuitive to use pId

% Get the config xml for pulling centers
CobraConfig = loadCfgXml;
 data.centers = [];
c=1;
for pid=data.activeCobras
    data.centers(c) = getARMval(CobraConfig,pid, 1,'Global_base_pos_x')+getARMval(CobraConfig,pid, 1, 'Global_base_pos_y')*i;
    data.Rptrl(c) = 80;
    c=c+1;
end

for pid=data.fixedCobras
    data.centers(c) = getARMval(CobraConfig,pid, 1,'Global_base_pos_x')+getARMval(CobraConfig,pid, 1,'Global_base_pos_y')*i;
    data.Rptrl(c) = 10;
    c=c+1;
end

c=1;
for pid=data.activeCobras
    data.sumLinkLengths(c) = getARMval(CobraConfig,pid, 1,'Link1_Link_Length')+getARMval(CobraConfig,pid, 1, 'Link2_Link_Length'); 
    c=c+1;
end

for pid=data.fixedCobras
    data.sumLinkLengths(c) = getARMval(CobraConfig,pid, 1,'Link1_Link_Length')+getARMval(CobraConfig,pid, 1, 'Link2_Link_Length'); 
   
    c=c+1;
end



% from actual target script 
% cmd_setHornMethodFiducialCoordinate 72.71819, 362.4398, 10
% cmd_setHornMethodFiducialCoordinate 461.49845, 733.04543, 10
% cmd_setHornMethodFiducialCoordinate 850.61126, 1104.9811, 10
% cmd_setHornMethodFiducialCoordinate 1238.5277, 1476.343, 10
% cmd_setHornMethodFiducialCoordinate 1628.8039, 1849.4198, 10

% Add the fiducial positions to centers array
data.centers = [data.centers.';
72.7181901965125 + 362.439887112847i;
461.498446356842 + 733.045426660864i;
850.611264646654 + 1104.98112670259i;
1238.52772825481 + 1476.34296141023i;
1628.80387570267 + 1849.41980465647i];
keyboard;
A =  [data.activeCobras, data.fixedCobras, 101, 102, 103, 104, 105];
for j=1:length(A)
plot(data.centers(j), 'rx')
cmplx(@text, data.centers(j), num2str(A(j)));
if(j<=length(data.sumLinkLengths))
circle(real(data.centers(j)), imag(data.centers(j)), 1.07*data.sumLinkLengths(j));
end
end