function [x_out y_out]=rotatexy(x,y,angle);
% rotate position data x and y through angle [radians]

% force x and y into row vectors
x = x(:)';
y = y(:)';

M_rotate = [ cos(angle) -sin(angle);...
             sin(angle)  cos(angle)];

v_out = M_rotate*[x; y];
x_out = v_out(1,:);
y_out = v_out(2,:);
