function output = appendLegend(plotHandle,legendEntry,varargin)
    
    % Get object handles
    lh1 = legend;

    if length(lh1)==0
%         output = legend(plotHandle,legendEntry,varargin{:});
        output = legend(legendEntry,varargin{:});
    else
        % Add object with new handle and new legend string to legend
        legendString = lh1.String;
        legendString{length(legendString)+1} = legendEntry;
%         output = legend([OUTH;plotHandle],OUTM{:},legendEntry,varargin{:});
        output = legend(legendString, varargin{:});
    end

    return
end