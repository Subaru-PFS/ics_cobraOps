function output=horn87(left,right);
% finds horn solution to translate left coordinates into right coordinates
% this is written for the case where all data lie in a plane.  to
% ensure this assumption, the code works exclusively on data
% presented as column vectors of complex numbers.
%
% Right ~= Left*exp(i*R)*S + T

if isrow(left)
  left = transpose(left);
end
if isrow(right)
  right = transpose(right);
end

lmean = mean(left);
rmean = mean(right);

lcm = left  - lmean;
rcm = right - rmean;
ndata = length(lcm);

% OPTIMIZE SCALE FACTOR Suv2xy (see BKP Horn)
Sr = sum( abs(rcm).^2 );
Sl = sum( abs(lcm).^2 );
output.S  = sqrt(Sr/Sl); % left2right

% OPTIMIZE ROTATION

% M for optimum rotation (Sij's) (sum of outer products)

% BPK's general solution (in 2D)
% r.l = [real(lcm)'; imag(lcm)'; zeros(1,ndata)];
% r.r = [real(rcm)'; imag(rcm)'; zeros(1,ndata)];
% M = r.r * r.l';
% N = [       trace(M),               0,                0, M(1,2) - M(2,1); ...
%                    0, M(1,1) - M(2,2),  M(1,2) + M(2,1),               0; ...
%                    0, M(1,2) + M(2,1), -M(1,1) + M(2,2),               0; ...
%      M(1,2) - M(2,1),               0,                0,       -trace(M)];
% [V, D] = eig(N);
% theta0 = acos(V(1,4))*2 - 2*pi;
% clear D;

% exp(i*R) rotates L into R
Rvec = sum(conj(lcm).*rcm);
output.R = angle(Rvec);


% OPTIMIZE TRANSLATION
output.T = rmean - output.S * exp(i*output.R) * lmean;

output.mean_err =  mean(abs(left*exp(i*output.R)*output.S + output.T - right));
output.stdev    =  std(left*exp(i*output.R)*output.S + output.T - right)/sqrt(2);
return

% for verification of the angle
%phi = -theta-.001:.0001:-theta+.001;
%for jj=1:length(phi)
%  D(jj) = sum(abs(rcm).*abs(lcm).*...
%              cos(angle(lcm) - angle(rcm*exp(i*phi(jj)))));
%end
%plot(phi,D,'x');
