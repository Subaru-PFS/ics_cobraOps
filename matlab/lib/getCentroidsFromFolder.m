%% Function for obtaining centroids of backlit fibers in a series of images
%% HOW TO USE
%   Input 'data' structure must have the following fields populated:
%         - data.ImgDir
%         - data.centers
%         - data.pId

function output=getCentroidsFromFolder(data)
start_time=tic;

% [data.ImgDir data.ImgFilter]
imagefile = dir2cell([data.ImgDir data.ImgFilter]);
Nimg = length(imagefile);
imgorder = {};

%% Read Images and Centroid

% if ~exist('Ptemp','var')
  for ii=1:Nimg;
%    I=fitsread([data.ImgDir imagefile{ii}]);
    Img=fitsread([imagefile{ii}]);
    disp(imagefile{ii}) 
    imgorder{ii} = imagefile{ii};
%     Multiple_MP_Centroiding(I,data.centers,false,false)
%     Test{ii} = Multiple_MP_Centroiding(I,data.centers,false,false)
    try
        Ptemp{ii}= struct(imgGetCntrds(Img,data.centers));
    catch
        keyboard;
    end
   end
% end
data.imgorder = imgorder;
% P=Ptemp;
% 
% 
%% Create Data Structure for Centroids
c = 1; %Start counter
for ii=data.pId
    fldID=sprintf('pId%d',ii);  
    pfld = sprintf('center%d',ii);
    intnsID = sprintf('intns%d',ii);
    TOGGLE = false;
    for jj = 1:length(Ptemp)
        if isfield(Ptemp{jj},pfld)
            TOGGLE = true;
            if isfield(data,fldID)
                data.(fldID).CCDpos = [data.(fldID).CCDpos Ptemp{jj}.(pfld)];
                data.(fldID).pxSum = [data.(fldID).pxSum Ptemp{jj}.(intnsID)];
            else
                data.(fldID).CCDpos = [Ptemp{jj}.(pfld)];
                data.(fldID).pxSum = [Ptemp{jj}.(intnsID)];
            end
        end
    end
    if TOGGLE
      data.npos(c) = length(data.(fldID).CCDpos);
      c = c+1;
    end
end

data.Nimg = Nimg;
% output = Ptemp;
output = data;
return;