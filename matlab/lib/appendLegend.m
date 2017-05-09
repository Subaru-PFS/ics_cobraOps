function output = appendLegend(plotHandle,legendEntry,varargin)
% Note: this was written by Chaz Morantz for Matlab ~2014a.  At
% 2016b, "legend" updates to include new entries in the figure.
% At the moment, this does not reflect the change.
%
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