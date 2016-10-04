function output = circle(varargin)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
varargin;
x = varargin{1};
y = varargin{2};
r = varargin{3};

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
if nargin>=4
    output = plot(x+xp,y+yp,varargin{4:end});
else
    output = plot(x+xp,y+yp);
end

end