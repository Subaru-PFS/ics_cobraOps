%% Assign centroids to Positioner
% return centroids assigned to the nearest centers found for them

%% INPUTS
% cntrds is a 1xn array of centroids as imaginary numbers
% centers is a nx1 array of cobra centers as imaginary numbers

function output = assignCntrds2Cobras(cntrds,centers,Rptrl)
% For each centroid
 
if length(Rptrl) == 1
    Rptrl = repmat(Rptrl,size(cntrds));
end

for thiscntrd=cntrds
    
    % Find the distance to each center
    dist2centers = abs(centers - thiscntrd);

    % Determine the index in center array of minimum distance
    [unused, indxMIN]=min(dist2centers);
    
    % Sanity Check that centroid is within patrol region radius of cobra
    % This weeds out the neighbors
    if dist2centers(indxMIN) > Rptrl(indxMIN)
        continue
    end
    
    % Declare the centerID for this centroid
    centerID=sprintf('center%d',indxMIN);
    
    % Assign centroids to output structure by 
    if exist('output','var') & isfield(output,centerID)
        output.(centerID)=(horzcat(output.(centerID),thiscntrd)); 
    else
        output.(centerID)=(thiscntrd);
    end
end

return;
