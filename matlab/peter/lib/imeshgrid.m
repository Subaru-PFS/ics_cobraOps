function xy=imeshgrid(xgv,ygv)
% 2D meshgrid output as a 1d complex column vectoc

  nx = length(xgv);
  ny = length(ygv);

  [x y] = meshgrid(xgv,ygv);
  
  xy = x + i*y;
  
  xy = reshape(xy,nx*ny,1);
  