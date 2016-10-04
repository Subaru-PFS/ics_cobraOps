function handle=plotellipse(x,y,a,b,phi)
% h=plotellipse(x,y,a,b,phi)
% plots a ellipse at x,y or at each x,y (if vectors)
%%%%%
unfinished code
if ~exist('psym','var')
  psym = '';
end
x=x(:);
y=y(:);
r=r(:);

if length(a) == 1
  a=a*ones(size(x));
end
if length(b) == 1
  b=b*ones(size(x));
end

theta = 2*pi*(0:1e-3:1); % row vector

circles = bsxfun(@plus, r*exp(i*theta), x + i*y);
circles = transpose(circles);

wasnotheld = ~ishold;
hold on;
handle = plot(real(circles),imag(circles),psym);
axis equal

if wasnotheld
  hold off;
end
