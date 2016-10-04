function qout=qcomp(v,phi)
% compose quaternion from vector + angle

v=v(:)/norm(v);
qout = [v*sin(phi/2) ; cos(phi/2)];
