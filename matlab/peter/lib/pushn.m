function arr = pushn(inarray, nn, val)
% push val into the nnth position of inarray
% usage newarr=pushn(oldarr, nn, val);
  
insize = size(inarray);
[lngth dim] = max(insize);

if (nn > (lngth + 1))
  arr = inarray;
  warning('position requested > length + 1');
  return;
elseif (nn < 1)
  arr = inarray;
  warning('position requested < 1');
  return;
end

outsize = [1 1];
outsize(dim) = lngth + length(val);


inarray = inarray(:);

arr = [inarray(1:nn-1); val(:); inarray(nn:end)];
arr = reshape(arr,outsize);

end