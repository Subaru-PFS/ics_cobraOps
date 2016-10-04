%% Inputs and Outputs
%-Inputs-
% IMG: input image [tiff]
% CENTERS: Positioner centers in the form of imaginary numbers
% NSIGMA: Number of image standard deviations to add to image median for
% determing BW conversion threshold
% MINPX: Minimum connected pixels to be counted as a fiber image
% DIAG: True or False. Runs diagnostic plots

%-Outputs-
% raw_centroids= array of current positions, unordered [pixels]

function output=imgGetCntrds(IMG, CENTERS, RPTRL, NSIGMA, MINPX, DIAG)

%% Check for inputs and declare defaults if not given
if ~exist('NSIGMA','var'), NSIGMA = 5;end; % background threshold
if ~exist('MINPX','var'), MINPX = 5;end; % minimum pixels
if ~exist('DIAG','var'), DIAG = true;end; 
if ~exist('RPTRL','var'), RPTRL = 10^6;end; % If no patrol region radius matrix is given, raise it so catches everything 

global tbout;

%% Centroiding

% These are the centers given to start with
centers=CENTERS.';

% This is how far from given centers to look for centroids and associate
% with the positioner for that center
Rptrl = RPTRL.';

% Get gray threshold with Peter-math
imgsize = size(IMG);
imgstats = allstats(double(reshape(IMG,1,[])));
tbout.imgstats = imgstats;
grythrsh = imgstats.median + NSIGMA*imgstats.std;
tbout.grythrsh = grythrsh;

% Use gray threshold to make BW image
IMG_bw = IMG > grythrsh;
tbout.IMG_bw = IMG_bw;
% imagesc(IMG_bw)

% Use matlab function to find connected objects in BW image
objs=bwconncomp(IMG_bw);
tbout.objs = objs;

% Initialize centroid array
barycenters = [];

 
% For each connected object matlab found do:
for jj = 1:objs.NumObjects
    % Count how many pix are in the obj
    npix(jj) = length(objs.PixelIdxList{jj});
    % If obj pix are greater than MINPX
    if npix(jj) > MINPX
        % Find centroid
        [yy xx] = ind2sub(imgsize, objs.PixelIdxList{jj});
        zz = xx + i*yy;
        tbout.p{jj} = zz;
        ww = double(IMG(objs.PixelIdxList{jj}));
        tbout.I{jj} = ww;
%         intns(jj) =  sum(ww)/length(ww); %INTENSITY EXTRACTION
        thisbarycent = zz.' * ww / sum(ww);
        barycenters = [barycenters thisbarycent];
        
% w=25;
% by = round(imag(thisbarycent)-w);
% bx = round(real(thisbarycent)-w);
% window = IMG([by:by+2*w],[bx:bx+2*w]);
% % imagesc(window)
% intns(jj) = sum(sum(window));
intns(jj) = 1;
% if intns(jj)==0
%     keyboard;
% end
        
%     else
% %         barycent(jj) = NaN;
    end
end

tbout.centroids = barycenters;
%keyboard;
if(length(barycenters) > 0)
output = assignCntrds2Cobras(barycenters, centers, Rptrl);

[junk username] = system('whoami');
if strcmp(username(1:end-1), 'petermao')
  fldnames = fields(output);
  for jj=1:length(fldnames);
    newoutput(jj) = output.(fldnames{jj});
  end
  output = newoutput;
  return;
end


%INTENSITY EXTRACTION
c=1;
for thiscntr=barycenters
    % Find the distance to each center
    dist2centers = abs(centers - thiscntr);
    
    % Determine the index in center array of minimum distance
    [unused, I]=min(dist2centers);
    
    fldID = strcat('intns',num2str(I));
    
    output.(fldID) = intns(c);
    c=c+1;
end 
else
    output =[];
end

if DIAG
  figure(1001);
  imagesc(IMG);
  colorbar;
  hold on;
 plot(barycenters,'w.--');
  length(barycenters)
  
  
%   figure(1002)
%   plot(barycenters, 'ro');
%   hold on 
%   plot(output
%   plot(barycent.center2,'w.--');
%   plot(barycent.center3,'w.--');
%   plot(barycent.center4,'w.--');
%   plot(barycent.center5,'w.--');
 % keyboard;
end

return;

