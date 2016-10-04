function handle=plotcircle(x,y,r,psym)
% h=plotcircle(x,y,r,psym)
% plots a circle at x,y or at each x,y (if vectors)

if ~exist('psym','var')
  psym = '';
end
x=x(:);
y=y(:);
r=r(:);

if length(r) == 1
  r=r*ones(size(x));
end

theta = 2*pi*(0:1e-2:1); % row vector

circles = bsxfun(@plus, r*exp(i*theta), x + i*y);
circles = transpose(circles);

wasnotheld = ~ishold;
hold on;
handle = plot(real(circles),imag(circles),psym);
axis equal

if wasnotheld
  hold off;
end
