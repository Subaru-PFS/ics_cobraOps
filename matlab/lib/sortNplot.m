function output = sortNplot(x,y, varargin)
    x = reshape(x,length(x),1);
    y = reshape(y,length(y),1);
    M = [x,y];
    M = sortrows(M,1);
    output = plot(M(:,1),M(:,2),varargin{:});
end