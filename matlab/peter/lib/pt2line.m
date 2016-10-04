function [dist intersection]=pt2line(L1,L2,PT)
% calculate distance from PT to line defined by point L1 and L2
% in two dimensions

if (length(L1) > 1)
  L1 = L1(1) + i*L1(2);
end
if (length(L2) > 1)
  L2 = L2(1) + i*L2(2);
end
if (length(PT) > 1)
  PT = PT(1) + i*PT(2);
end


line1 = L2 - L1;
point = PT - L1;

rotated_point = point * exp(-i*angle(line1));
dist = abs(imag(rotated_point));
intersection = real(rotated_point) * exp(i*angle(line1)) + L1;