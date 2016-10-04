function [output ST]=zerneval(zfit,xy)
% evaluate at points (XY) the Zernike-based fit.
% 
% see ZERNFIT.M and Chunyu Zhao's ZernikeTools

  jST = zfit.jST;
  Ast = zfit.Ast;
  xy = (xy - zfit.Trans)/zfit.Scale;
  
  jS = jST(jST > 0);
  jT = -jST(jST < 0);

  for jj=1:length(jS)
    SS(:,jj) = S_xy(jST(jj),xy);
  end
  
  ST = [SS -i*SS(:,jT-1)];
 
  output = ST*Ast;

  return
  
  %% evaluates the orthogonality/normalization of S and T functions
  %% evaluated on the positions listed in the XY input.
  fig = gcf + 100;
  figure(fig+1);
  imagesc((abs(ST'*ST)));
  axis equal;
  colorbar;
  figure(fig+2);
  plot(diag(ST'*ST));
  
  
  figure(fig-100)