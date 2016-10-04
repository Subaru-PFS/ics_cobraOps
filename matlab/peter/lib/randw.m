function w = randw(p0, p1, n)
% USAGE: w = randw(p0, p1, n)
% generate a 1xn vector with random numbers in a wedge shaped
% distribution on the interval [0,1]
%   p0 = relative probability at low end
%   p1 = relative probability at high end
%   n == number of rv's

alpha = 2/(p0 + p1);
a2 = alpha * (p1 - p0);
b  = alpha * p0;

% $$$ x = 0:.01:1;
% $$$ plot(x, a2/2*x.^2 + b*x); hold on;
% $$$ plot(x,(-b + sqrt(b^2+2*a2*x))/a2,'r');
% $$$ hold off

u = rand(1,n);

if (a2 ~= 0)
  w = (-b + sqrt(b^2 + 2*a2*(u)))/(a2);
else 
  w = u;
end

% $$$ [h,x ] = hist(w,100);
% $$$ bar(x,h / (n*mean(diff(x))),'hist')
% $$$ fprintf(1,'%f %f\n', alpha*p0, alpha*p1);
