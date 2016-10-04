%% Function for obtaining centroids of backlit fibers in a single image
%% HOW TO USE
%   Input 'data' structure must have the following fields populated:
%         - data.ImgDir
%         - data.centers
%         - data.pId

function output=getCntrdsFromImage(data,NSIGMA)
start_time=tic;

imagefile = data.Name;

%% Read Image and Centroid
I=fitsread(fullfile(data.ImgDir, imagefile));
disp(imagefile)

if exist('NSIGMA','var')
    Ptemp = struct(imgGetCntrds(I,data.centers,data.Rptrl,NSIGMA));
else
    Ptemp = struct(imgGetCntrds(I,data.centers,data.Rptrl));
end


xx = 0;
c = 1;
for ii=data.pId % Count used dots.
    pfld = sprintf('center%d',c);
     if isfield(Ptemp,pfld)
    xx = xx + length( Ptemp.(pfld));
     end
      c = c+1;
end
xx
%% Create Data Structure for Centroids
c = 1; %Start counter
%keyboard;
figure(1002)
clf;
hold on;

%Create list in the same order as the centers were given. 
allIds = [data.activeCobras, data.fixedCobras, data.fiducials];
for ii=allIds % CHANGED from data.pId @Johannes
    fldID=sprintf('pId%d',ii);  
    pfld = sprintf('center%d',c);
    intnsID = sprintf('intns%d',c);
    data.(fldID).CCDpos = [];
    if isfield(Ptemp,pfld)
        if isfield(data,fldID)
            plot(Ptemp.(pfld),'ro');
            text(real(Ptemp.(pfld)), imag(Ptemp.(pfld)),fldID);
            data.(fldID).CCDpos = [data.(fldID).CCDpos Ptemp.(pfld)];
%             data.(fldID).pxSum = [data.(fldID).pxSum Ptemp.(intnsID)];
        else
            data.(fldID).CCDpos = [Ptemp.(pfld)];
%             data.(fldID).pxSum = [Ptemp.(intnsID)];
        end
    end
%     keyboard;
    data.npos(c) = length(data.(fldID).CCDpos);
    if data.npos(c) == 0
        sprintf('No centroids found for pid%d',ii)
    else
        sprintf('%d Centroids found for pid%d',data.npos(c),ii)
    end
    c = c+1;
end 

 
output = data;

return;