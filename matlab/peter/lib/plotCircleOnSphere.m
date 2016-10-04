function handle=plotCircleOnSphere(ra,dec,arcdist,psym)
% ra and dec in radians!
% arcdist in radians!
% assumes a unit sphere
% 
% requires QVROT.M
% http://www.mathworks.com/matlabcentral/fileexchange/1176-quaternion-toolbox
if ~exist('psym','var')
  psym = '';
end

ra = ra(:);
dec = dec(:);
arcdist = arcdist(:);

[cx cy cz] = sph2cart(ra,dec,1);
[rx ry rz] = sph2cart(ra,dec+arcdist,1); % matlab appears to be handling dec+arcdist>pi/2 properly

centers = reshape([cx cy cz]',[],1);
edges = [rx ry rz];
phi = (2*pi*(0:5e-2:1));
qxyz = centers * sin(phi/2);
qxyz = permute(reshape(qxyz,3,length(ra),length(phi)), [3 1 2]);
qcos = cos(phi(:)/2);

washeld = ishold;
hold on;

for jj=1:length(ra)
  xyzr = qvrot([qxyz(:,:,jj) qcos], edges(jj,:));
  plot3(xyzr(:,1), xyzr(:,2), xyzr(:,3),psym);
end

if ~washeld
  hold off;
end
