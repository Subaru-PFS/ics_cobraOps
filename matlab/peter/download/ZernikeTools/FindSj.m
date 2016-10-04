function [Z_Sjx C_Sjx Z_Sjy C_Sjy] = FindSj(j)

% for a given index j, find the corresponding Zernike indices and
% coefficients for both x and y component of Sj

threshold = 0.00001; % coefficient whose magnitudes are neither 1/2 or 1/sqrt(2) are discarded

if j==1
    disp('j=1 is trivial case!');
    Z_Sjx = 1;
    C_Sjx = 0;
    Z_Sjy = 1;
    C_Sjy = 0;
    return;
end
    
[n m q] = FindNMQ(j);
 
if n==m
    % Sj = 1/sqrt(2*n*(n+1))Gradient(Zj)
    C = ones(j, 1)/sqrt(2*n*(n+1));
    C_Sjx = DerivX_Zern(j, j).*C;
    C_Sjy = DerivY_Zern(j, j).*C;
    Z_Sjx = find(abs(C_Sjx)>threshold);
    Z_Sjy = find(abs(C_Sjy)>threshold);
    C_Sjx = C_Sjx(Z_Sjx);
    C_Sjy = C_Sjy(Z_Sjy);
    
else
    % Sj = 1/sqrt(4*n*(n+1))Gradient(Zj) - 1/sqrt(4*n*(n-1))Gradient(Zj'(n'=n-2,m'=m))
    C1 = ones(j, 1)/sqrt(4*n*(n+1));
    C_Sjx1 = DerivX_Zern(j, j).*C1;
    C_Sjy1 = DerivY_Zern(j, j).*C1;
    
    j2 = FindJ(n-2, m);
    if length(j2)==1
        jp = j2;
    else
        jp = j2(find(~mod(j2-[j j], 2)));
    end
    C2 = ones(jp, 1)/sqrt(4*n*(n-1));
    C_Sjx2 = DerivX_Zern(jp, jp).*C2;
    C_Sjy2 = DerivY_Zern(jp, jp).*C2;
    C_Sjx2 = vertcat(C_Sjx2, zeros(j-jp, 1));
    C_Sjy2 = vertcat(C_Sjy2, zeros(j-jp, 1));
    
    C_Sjx = C_Sjx1 - C_Sjx2;
    C_Sjy = C_Sjy1 - C_Sjy2;
    
    Z_Sjx = find(abs(C_Sjx)>threshold);
    Z_Sjy = find(abs(C_Sjy)>threshold);
    C_Sjx = C_Sjx(Z_Sjx);
    C_Sjy = C_Sjy(Z_Sjy);
end