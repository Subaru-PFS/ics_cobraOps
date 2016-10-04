function output = symmetrize(dist);
    dist.dmatrix = min(dist.dmatrix, dist.dmatrix');

    for kk = 1:length(dist.dst)
        dist.dst(kk) = dist.dmatrix(dist.rc(kk,1), dist.rc(kk,2));
    end

    output = dist;

end

% 12/4/15:  not currently used anywhere