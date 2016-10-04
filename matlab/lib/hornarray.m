function output = hornarray(left, right, subarray)
% runs horn87 on an mxn array (m measurements, n fibers)
% LEFT is the data, RIGHT is the reference

if ~exist('subarray','var'), subarray = 1:length(right), end;

for jj=1:size(left,1)
  hh = horn87(left(jj,subarray), right(subarray));
  data(jj,:) = left(jj,:) * exp(i*hh.R) * hh.S + hh.T;
  horn(jj) = hh;
end

output = packstruct(data,horn);
