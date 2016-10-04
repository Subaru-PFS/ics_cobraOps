function newarr = poparr(oldarr, indices)
% remove elements from array.  no error checking yet
  
  keepers = ones(size(oldarr));
  keepers(indices) = 0;
  keepers = find(keepers);
  newarr = oldarr(keepers);
  
end
