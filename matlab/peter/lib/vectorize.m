function vec=vectorize(array,opt)
% turn ndim array into a row or column vector

  if ~exist('opt','var'), opt = 'r'; end
  
  vec = reshape(array, 1, numel(array));


  if ~strcmp('opt','r')
    vec = vec(:);
  end

 