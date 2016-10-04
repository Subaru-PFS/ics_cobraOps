function output=ellipsefit(XY)
% fit data to an ellipse

cfit = circfit(XY);

XYC = XY - cfit.c;
RR = abs(XYC);
THT = angle(XYC);

ellipse_radius = @(a,b,p,x)a*b./sqrt((b*cos(x-p)).^2 + (a*sin(x- ...
                                                  p)).^2);

RR = RR(:);
THT = THT(:);

[temp j_maxR] = max(RR);


fitellipse = fit(THT,RR,ellipse_radius,...
                 'startpoint',[cfit.R, cfit.R, THT(j_maxR)]);
output.fit = fitellipse;
output.ecc = sqrt(1 - (fitellipse.b/fitellipse.a)^2);
output.X = THT;
output.Y = RR;
output.c = cfit.c;

plot(fitellipse,THT,RR,'resid');