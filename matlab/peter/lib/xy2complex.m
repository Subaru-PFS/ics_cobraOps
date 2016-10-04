function c=xy2complex(array, dim)
% convert mx2 or 2xn array into mx1 or 1xn complex numbers
% from first two columns or rows.v
% optional DIM specifies the dimension to keep (1,2...)
asize = size(array);

if (length(asize) >  2) 
  warning('this code was not written for arrays with > 2 dimensions');
end

if ~exist('dim','var');
  [minlength mindim] = min(asize);
  [maxlength maxdim] = max(asize);
  if (minlength >=2)
    dim = maxdim;
  else % minlength == 1
    dim = mindim;
  end
  % note in a 2x2 case, dim = 1;
end

if (dim==1)
  c = complex(array(:,1),array(:,2));
elseif (dim==2)
  c = complex(array(1,:),array(2,:));
else
  warning('dim should be 1 or 2');
  c = NaN;
end

