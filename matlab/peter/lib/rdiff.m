function output=rdiff(Matrix, iterations, dimension)
% USAGE: Y = rdiff(X, n, dim)
%
% ratio difference between elements in a Matrix along the specified
% dimension.  basically, this is "diff" with division instead of
% subraction.
%
% EXAMPLE:
% >> rdiff(1:7)                                               
% ans =                                                       
%     2.0000    1.5000    1.3333    1.2500    1.2000    1.1667
% >> rdiff(2.^(1:7))                                          
% ans =                                                       
%     2.0000    2.0000    2.0000    2.0000    2.0000    2.0000
%
% Peter H Mao, Caltech 2016-09-15 

  if ~exist('iterations','var'), iterations = 1; end;
  if ~exist('dimension','var'), dimension = 1; end;

  if iterations > 1
      disp('Warning (RDIFF): sign may not be correct for multiple diffs')
  end
  
  transposed = false;
  
  if isrow(Matrix) & dimension == 1
      Matrix = Matrix(:);
      transposed = true;
  end
  
  sgn = (diff(sign(Matrix),iterations,dimension) == 0) * 2 - 1;

  ratios = exp(diff(log(abs(Matrix)),iterations,dimension));

  
  output = ratios.*sgn;
  
  if transposed
      output = (ratios.*sgn).';
  end
