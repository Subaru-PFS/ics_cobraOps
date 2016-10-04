function output=cmplx(funhandle,varargin)
% complex wrapper for functions that take the form F(x,y,...)

if isreal(varargin{1}) % normal case
  try  
    output = funhandle(varargin{:})
  catch
    funhandle(varargin{:})
  end
else
  try
    output = funhandle(real(varargin{1}), imag(varargin{1}), varargin{2:end});
  catch
    funhandle(real(varargin{1}), imag(varargin{1}), varargin{2:end});
  end
end
