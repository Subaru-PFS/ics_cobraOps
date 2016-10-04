function output=zernfit(xy,dxy,orders)
% fit vector vector field DXY at locations XY to Zernike-derived
% vector basis functions to ORDERS
%
% This is suitable for vector fields on the unit circle or subset
% thereof.
% 
% Uses Chunyu Zhao's ZernikeTools modified by PHM to accept and
% return complex variables at the top level.
  
  xycenter = ( max(real(xy)) + min(real(xy))  + ...
               i*( max(imag(xy)) + min(imag(xy)) ) ) / 2;
  xy = xy - xycenter;
  scale  = max(abs(xy))*1.01; % give 1% padding so everything is
                              % inside the unit circle.
  xy = xy / scale;

  % force xy and dxy into column vectors
  xy = xy(:);
  dxy = dxy(:);
  
  %% determine the T functions that are independent of the S
  %% functions
  nn = 1:orders;
  % http://oeis.org/A117142
  laplacians = (2 * nn.^2 + 10 * nn + 3 + (-1).^nn .* (2 * nn - 3))/16;
  laplacians = laplacians(laplacians < orders);
  solenoidal = logical(ones(1,orders));
  solenoidal(laplacians) = 0;
  % indices of the S and T functions
  jT = find(solenoidal);
  jS = 2:orders;
  
  for jj = 2:orders
    SS(:,jj-1) = S_xy(jj,xy);
  end
  
  ST = [SS -i*SS(:,jT-1)]; % T = -i * S
  
  Ast = ipnorm(dxy, ST)';
  dxy_fit = ST*Ast;

  %% outputs
  output.jST = [jS -jT];
  output.Ast = Ast;
  output.Trans = xycenter;
  output.Scale = scale;
  output.err = std(dxy - dxy_fit)/sqrt(2); % 1-sigma standard error
                                           % in fit
end

function output=ipnorm(x, y)
% calculates the inner product of x and y, normalized by y

  output = real(x' * y) / real(y' * y);
end