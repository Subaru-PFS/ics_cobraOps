function output = XY2TP(xy, L1, L2, verify)
% with center at origin, calculates the cobra arm THT, PHI pair
% from the cartesian coordinate. 
% 
% xy is a complex number of the form X + i*Y and may be an ARRAY of
% such numbers.

if ~exist('L1','var')
  L1 = 2.375;
end
if ~exist('L2','var')
  L2 = 2.375;
end
  
threshold = mean(L1)/22;

R = abs(xy);
T = angle(xy);

% check for inside/outside limits errors
dr_outside = bsxfun(@minus, R, L1 + L2);
dr_inside  = bsxfun(@minus, R, abs(L1 - L2));
dr_outside_violation = find(dr_outside > threshold);
dr_inside_violation  = find(dr_inside < -threshold);
for jj=1:length(dr_outside_violation)
    kk = dr_outside_violation(jj);
    [row, col] = ind2sub(size(R), kk);
    fprintf(1, ['XY2TP detected illegal location outside patrol area (%d,%d) ' ...
                'by %g\n'], row, col, dr_outside(kk));
end
for jj=1:length(dr_inside_violation)
    kk = dr_inside_violation(jj);
    [row, col] = ind2sub(size(R), kk);
    fprintf(1, ['XY2TP detected illegal location inside patrol area (%d,%d) ' ...
                'by %g\n'], row, col, dr_inside(kk));
end

if size(R,2) == 1
    phi = -acos((R.^2 - L1.^2 - L2.^2)./(2*L1.*L2));
    tht_in = acos((-L2.^2  + R.^2 + L1.^2)./(2*L1.*R));
else
    phi = -acos(bsxfun(@rdivide, ...
                       bsxfun(@minus, R.^2, L1.^2 + L2.^2), ...
                       2 * L1 .* L2));
    tht_in = acos(bsxfun(@rdivide,...
                         bsxfun(@plus, R.^2, L1.^2 - L2.^2),...
                         bsxfun(@times, R  , 2 * L1)));
end

IN_PATROL_AREA = dr_outside < 0 & dr_inside > 0;
tht_in(IN_PATROL_AREA) = real(tht_in(IN_PATROL_AREA));
tht_in(~IN_PATROL_AREA) = real(tht_in(~IN_PATROL_AREA)); % nan;

tht = mod(T - tht_in.*sign(real(phi)) + pi, 2*pi)-pi;

output.tht = tht;
output.phi = phi;

if exist('verify','var');
    [T tht_in tht]
    plot(xy,'ko','MarkerFace','k');
    axis equal
    plotcircle(0,0,L1+L2,'k:');
    hold on;
    plot([complex(0) L1* exp(i*tht) L1*exp(i*tht) + L2*exp(i*(tht+phi))], ...
         '-o');
    hold off;
end
