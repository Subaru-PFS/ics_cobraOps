function output=loadxyz(file,varargin)
% fancier version of "load"
% usage : output=loadxyz('filename','field1','field2',...)
% default: output=loadxyz('filename');
%          output.x = first column
%          output.y = second, etc.
%          output.c = x+i*y
%          output.m = [x y z]

d = load(file);
ncol = size(d,2);
output.dat = d;

if (length(varargin) < ncol)
  fieldlabel = 'xyzabdefghijklnopqrstuvw';
  if ncol-length(varargin) > length(fieldlabel)
    warning('using automatic labels, but not enough field labels to go around!');
    return;
  end
  kk=1;
  for jj=length(varargin)+1:ncol
    varargin{jj} = fieldlabel(kk);
    kk=kk+1;
  end
end


for jj=1:length(varargin)
  if ~ischar(varargin{jj})
    warning(['argument ' num2str(jj+1) ' does not appear to be a string']);
    continue;
  end
  eval(['output.' varargin{jj} ' = d(:,' num2str(jj) ');']);
end

if ( isfield(output,'x') && isfield(output,'y') )
  if ~isfield(output,'c')
    output.c = output.x + i * output.y;
  else
    warning('field ''c'' is reserved -- change its name');
  end
  if isfield(output,'z')
    if ~isfield(output,'m')
      output.m = [output.x output.y output.z];
    else
      warning('field ''m'' is reserved -- change its name');
    end      
  end
end

