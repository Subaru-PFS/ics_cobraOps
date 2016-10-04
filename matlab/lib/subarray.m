function array_out=subarray(array_in, varargin)
% Author: Peter H. Mao, Caltech
%
% simple program to access an element of an array when the array is
% the output of another function
%
% usage array_out=subarray(array_in,dim1_indices,dim2_indices,....);
% to use the ':' symbol or the 'end' keyword, pass the range as a string
%
% intent: I call this function when I have a function that returns an array
% and I want to access some subset or element of that array.  this saves me from
% having a temporary variable in my workspace.

  if isempty(array_in) % special case for empty matrices
    array_out = array_in;
    return;
  end
  
  dimlength = size(array_in);
  NonSingletonDimensions = find(dimlength > 1);
  n_dimensions = length(NonSingletonDimensions);
  if n_dimensions
    if n_dimensions < length(varargin)
      warning('you may have too many index ranges specified');
    end
    
    kk = NonSingletonDimensions(1);
    for jj=1:length(varargin)
      if ischar(varargin{jj}) & ~strcmp(varargin{jj},':')
        range = regexprep(varargin{jj},'end', 'dimlength(kk)');
        varargin{jj} = eval(range);
      end
      kk=kk+1;
    end
    
    array_out = array_in(varargin{:});
  else
    array_out = array_in(varargin{1}); % special case for scalars
  end
  
  return;
