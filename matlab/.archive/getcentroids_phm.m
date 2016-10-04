function output=getcentroids_phm(fitsfile,NSIGMA);
% example code to get centroids by Peter Mao

  if ~exist('NSIGMA','var'), NSIGMA = 5;end; % background threshold
  minNumPix = 5;
  
% $$$   img = fitsread(fitsfile);
  img = imread(fitsfile);
  
  imgsize = size(img);
  
  %% this is junk
  level = graythresh(img);
  bw_otsu = im2bw(img,level);
  
  figure(980)
  colormap('gray');
  imagesc(bw_otsu);
  axis equal;
  title('B/W image from thresholding on Otsu''s level');
  %%                                                     
  
  imgstats = allstats(double(vectorize(img)));
  
  mythresh = imgstats.median + NSIGMA*imgstats.std
  bw_mine  = img > mythresh;
  
  objs = bwconncomp(bw_mine);
  
  
  for jj = 1:objs.NumObjects
    npix(jj) = length(objs.PixelIdxList{jj});
    if npix(jj) > minNumPix
      [yy xx] = ind2sub(imgsize, objs.PixelIdxList{jj});
      zz = xx + i*yy;
      ww = double(img(objs.PixelIdxList{jj}));
      barycent(jj) = zz.' * ww / sum(ww);
    else
      barycent(jj) = NaN;
    end
  end
  
  output.bc = barycent;
  output.npix = npix;
 
  %% plot
  figure(981)
  imagesc(img); hold on;
  colormap('gray');
  axis equal;
  colorbar;
  
  plot(barycent,'ro');
  hold off
  title(sprintf('simple barycenter using mean + %d \\sigma threshold', NSIGMA));