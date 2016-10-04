function [arr val] = popn(inarray, nn)
% pop val from the nnth position of inarray
% usage  [newarray value] = popn(oldarray, nn);
  
insize = size(inarray);
[lngth dim] = max(insize);
val = NaN;

if (nn > (lngth))
  arr = inarray;
  warning('position requested > length');
  return;
elseif (nn < 1)
  arr = inarray;
  warning('position requested < 1');
  return;
end

outsize = [1 1];
outsize(dim) = lngth - 1;

inarray = inarray(:);

val = inarray(nn);
arr = [inarray(1:nn-1); inarray(nn+1:end)];
arr = reshape(arr,outsize);

end