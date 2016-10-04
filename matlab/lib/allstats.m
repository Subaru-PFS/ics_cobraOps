function s=allstats(x,dim,flag)
% ALLSTATS All Common Statistics.
% S = ALLSTATS(X,DIM,Flag) returns common statistics of the data in X along
% the dimension DIM, with the variance type specified by Flag.
%
% X must be non-sparse, real, single or double, and 2D. NaNs contained in
% X are consdered missing data and are ignored in computations. Note that
% when X contains NaNs, this function will return different results than
% those returned by the standard MATLAB functions MEAN, MEDIAN, and MODE,
% STD and VAR because these function do not ignore NaNs.
%
% If DIM is [] or not given, DIM is the first non-singleton dimension of X.
%
% If Flag = 0 or is not given, the variance is normalized by N-1.
% If Flag = 1, the variance is normalized by N.
%
% The output S is a structure containing the statistics with field names
% describing the respective statistics:
%
% S.n = number of observations (excluding NaNs)
% S.nan = number of NaNs
% S.min = minimum
% S.q1 = 25th percentile
% S.median = median (50th percentile)
% S.q3 = 75th percentile
% S.max = maximum
% S.mode = mode
% S.mean = mean
% S.std = standard deviation
% S.var = variance
% S.skew = skewness
% S.kurt = kurtosis (= 3 for Gaussian or Normal Data)
%
% For vector data, outputs are scalars.
% For matrix data, outputs are row vectors for DIM=1
%                       and column vectors for DIM=2.
%
% Reference: http://www.itl.nist.gov/div898/handbook/
%
% Examples:
% ALLSTATS(X) for vector X computes statistics along the vector, for matrix
% X computes statistics for each column of X.
% ALLSTATS(X,1) for matrix X computes statistics along the row dimension,
% i.e., for each column of X.
% ALLSTATS(X,2) for matrix X computes statistics along the column
% dimension, i.e., for each row of X.
% ALLSTATS(X,[],1) computes statistics along the first non-singleton
% dimension of X, and normalizes the variance and standard deviation by N.
%
% See also MAX, MIN, MEAN, MEDIAN, MODE, STD, VAR

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-03-04
%
% The gracious assistance of John D'Errico and Urs Schwarz is hereby noted.

if nargin<1
   error('allstats:NotEnoughInputArguments',...
      'At Least One Input Argument Required.')
elseif ~isfloat(x)
   error('allstats:InvalidInput','Input X Must be Single or Double.')
elseif ndims(x)>2
   error('allstats:InvalidNumberofDimensions','X Must be 2D.')
elseif issparse(x)
   error('allstats:NotFull','X Must Not be Sparse.')
elseif ~isreal(x)
   error('allstats:NotReal','X Must be Real Valued.')
end
xsiz=size(x);
if nargin<3
   flag=0;
end
if nargin<2 || isempty(dim)
   dim=find(xsiz>1,1);
end
if isempty(dim) || numel(dim)~=1 || fix(dim)~=dim || dim<1 || dim>ndims(x)
   error('allstats:DimOutofRange','DIM Invalid.')
end
if isempty(flag) || numel(flag)~=1 || fix(flag)~=flag || (flag~=0&&flag~=1)
   error('allstats:InvalidInput','Flag Must be 0 or 1.')
end
      
N=xsiz(dim);
if N==1                 % one observation along chosed DIM, return defaults
   sv=zeros(xsiz,class(x));
   in=isnan(x);
   s=struct('n',double(~in),'nan',double(in),'min',x,'q1',x,'median',x,...
            'q3',x,'max',x,'mode',x,'mean',x,'std',sv,'var',sv,...
            'skew',sv,'kurt',sv);
   return
else                                             % create struct for output
   s=struct('n',[],'nan',[],'min',[],'q1',[],'median',[],'q3',[],'max',[],...
            'mode',[],'mean',[],'std',[],'var',[],'skew',[],'kurt',[]);
end

tflag=false;
if dim==2      % work along DIM=1, convert back later if needed
   x=x.';
   tflag=true; % flag to note transpose taken
end

s.max=max(x);                                                % built in max
s.min=min(x);                                                % built in min

% To avoid redundant error checking and therefore maximize speed, implement
% other statistical measures here rather than call their respective M-files.

s.mean=nanmean(x,1);                                                 % mean

xs=x-repmat(s.mean,N,1); % subtract mean from each column
xss=xs.*xs;
s.var=nanmean(xss,flag);                                              % var
s.std=sqrt(s.var);                                                    % std
s.skew=nanmean(xss.*xs,flag)./s.var.^(1.5);                          % skew
s.kurt=nanmean(xss.*xss,flag)./s.var.^2;                         % kurtosis

xs=sort(x);         % sort data down columns, this puts NaNs last
inan=isnan(xs);     % true for NaNs
s.nan=sum(inan);    % number of NaNs per column
s.n=N-s.nan;        % number of non-NaN observations per column
ic=s.nan==0;        % columns where there are no NaNs

tmp=nan+zeros(size(s.mean),class(x));   % template
s.q1=tmp;
s.median=tmp;
s.q3=tmp;

% process columns having no NaNs to find percentiles and median
if any(ic)
   q=[0, (0.5:N -0.5)/N, 1]';
   xx=[s.min(ic); xs(:,ic); s.max(ic)];
   xq=interp1q(q,xx,[1;2;3]/4); % fast interp
   s.q1(ic)=xq(1,:);                                                   % q1
   s.median(ic)=xq(2,:);                                           % median
   s.q3(ic)=xq(3,:);                                                   % q3
end
% now process columns having NaNs
for k=find(s.nan & s.n>3) % columns having fewer than 4 non-NaNs return NaN
      N=s.n(k);
      q=[0, (0.5:N -0.5)/N, 1]';
      xx=[s.min(k); xs(1:N,k); s.max(k)];
      xq=interp1q(q,xx,[1;2;3]/4);
      s.q1(k)=xq(1);                                                   % q1
      s.median(k)=xq(2);                                           % median
      s.q3(k)=xq(3);                                                   % q3
end

s.mode=tmp;% columns having fewer than 4 non-NaNs return NaN
cols=1:size(xs,2);
cols(s.n<4)=[];
for k=cols	% get mode column by column
   N=s.n(k);
	y=[1;diff(xs(1:N,k))]~=0;	% where distinct values begin
	m=diff([find(y);N+1]);     % counts 
   y=xs(y,k);                 % the unique values
   [idx,idx]=max(m);
   s.mode(k)=y(idx);                                                 % mode
end

if tflag && numel(s.min)~=1               % transpose data back for DIM = 2
   s.n=s.n.';
   s.nan=s.nan.';
   s.min=s.min.';
   s.q1=s.q1.';
   s.median=s.median.';
   s.q3=s.q3.';
   s.max=s.max.';
   s.mode=s.mode.';
   s.mean=s.mean.';
   s.var=s.var.';
   s.std=s.std.';
   s.skew=s.skew.';
   s.kurt=s.kurt.';
end
%--------------------------------------------------------------------------
function m=nanmean(x,flag)
% mean along row dimension ignoring NaN elements
% flag is zero to divide by N-1
% flag is one to divide by N

bnan=isnan(x);                         % find NaNs in each column
N=size(x,1)-sum(bnan);                 % number of non-NaNs in each column
x(bnan)=0;                             % set NaN elements to zero
m=nan+zeros(1,size(x,2),class(x));     % default NaN output
idx=N>1;                               % columns that are not all NaNs
m(idx)=sum(x(:,idx))./(N(idx)-1+flag); % results excluding NaNs
