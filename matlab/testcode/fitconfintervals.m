% test code to understand fit's confidence intervals

% fake data set
xx = rand(100,1) * 2 * pi;
yy = 0.2*randn*xx + 0.05 * randn(size(xx)).*sqrt(xx) +...
     ...%.01 * sin(xx + randn) + ...
     .001 * randn(size(xx));

wts = 1e-3 + .05*sqrt(xx);
%modelfn = fittype({'x','1'});

S   = sum(wts.^-2);
Sy  = sum(yy .* wts.^-2);
Sx  = sum(xx .* wts.^-2);
Sxx = sum(xx.^2 .* wts.^-2);
Sxy = sum(xx.* yy .* wts.^-2);
Dlt = S*Sxx - Sx^2;

nr.a = (S*Sxy - Sx*Sy)/Dlt;
nr.b = (Sxx*Sy - Sx*Sxy)/Dlt;
nr.sigma.a = sqrt(S/Dlt);
nr.sigma.b = sqrt(Sxx/Dlt);

nr1.a = Sxy/Sxx; % y=ax least squares fit
nr1.sig = sqrt(1/Sxx)
% $$$ b = 0
% $$$ ft = fittype(@(a,x) a*x + b,'Problem','b');
% $$$ [slopefit gof output] = fit(xx, yy, ft, 'Problem', 0)
[slopefit gof output] = fit(xx, yy, @(a,x) a*x,'Startpoint', [0], 'Weights', wts.^2);

fprintf(1, 'fit:   %.4f\n', chisq(yy,feval(slopefit,xx),wts)/gof.dfe);
fprintf(1, 'nr:    %.4f\n', chisq(yy,nr.a*xx + nr.b, wts)/gof.dfe);
fprintf(1, 'nr1:   %.4f\n', chisq(yy,nr1.a*xx,wts)/(output.numobs - 1));
plot(xx,yy,'.'); hold on; 
plot(xx,feval(slopefit,xx),'r'); 
plot(xx,nr.a*xx + nr.b,'g'); 
plot(xx,nr1.a*xx,'g--'); 
hold off;