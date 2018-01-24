function output=realizeTrajectory2(gt2,geom,useP,useL)
% realize a trajectory.  pick shorest route for each target
% acquisition, calculate the trajectory. plot takes output of
% generateTrajectory2 as input also requires the bench geometry
% 
% inputs:
% gt2    : output of generateTrajectory2
% geom   : output of defineBenchGeometry
% useP   : logical array indicating which cobras use the ss hard stop
%
% outputs:
% traj  : MxN array of global x+iy trajectories (early theta)
% trajL : MxN array of global x+iy trajectories (late theta)
% useP  : 1xM logical: theta moves out in positive direction (same
%         sense) for the primary trajectory
% useL  : 1xM logical: theta moves out late
% tht   : MxN local theta
% phi   : MxN local phi
% thtDT(L) : 1xM fwd bin shifts (0 == maximally early) >= 0
% phiDT : 1xM fwd bin shifts (0 == maximally late)  >= 0
% ltdiff: 1xM time difference over lmax for changing theta direction
    
    if ~exist('useP','var')
        % if not specified, use the shorter theta route
        useP = (gt2.nthtP < gt2.nthtN);
    end
    if ~exist('useL','var')
        useL = logical(zeros(size(useP)));
    end

    useN = ~useP;
    ltht = useP .* gt2.lthtP + useN .* gt2.lthtN;
    lphi = gt2.lphi;
    lmax = max([ltht lphi]);  % max # time bins needed

    ntht = useP .* gt2.nthtP + useN .* gt2.nthtN;
    nphi = gt2.nphi;
    nmax = max([ntht nphi]);  % max # steps needed


    nCobras = length(useP);
    
    % what is the penalty for using the longer theta path?
    ltdiff = (useN .* gt2.lthtP + useP .* gt2.lthtN) - lmax;
    
    Tht       = gt2.ThtP;
    Tht(useN) = gt2.ThtN(useN);
    
    ThtN = gt2.ThtN;
    Phi  = gt2.Phi;
    
    for jj = 1:nCobras
        % default trajectories: mixture of thtP and thtN, whichever is shorter
        tht(jj,:) = padarray(Tht{jj}, [0, lmax - ltht(jj)], 'replicate','post');
        phi(jj,:) = padarray(Phi{jj}, [0, lmax - lphi(jj)], 'replicate','pre');
        % it only makes sense to run tht late in the thtN case (L == late)
        if lmax >= gt2.lthtN(jj)
            thtL(jj,:) = padarray(ThtN{jj}, [0, lmax - gt2.lthtN(jj)], 'replicate','pre');
        else
            % also, if lthtN(jj) > lmax, then we need to lengthen lmax
            % this is a subset of ltdiff > 0.
            thtL(jj,:) = nan(1,lmax);
        end
    end
    
    % the early/late trajectory (ss or os, depending on useP)
    traj = (bsxfun(@times, geom.L1, exp(1i*tht)) + ...
            bsxfun(@times, geom.L2, exp(1i*(tht+phi)) ) );
    traj = bsxfun(@plus, geom.center, traj);

    % the late/late trajectory (os only)
    trajL = (bsxfun(@times, geom.L1, exp(1i*thtL)) + ...
            bsxfun(@times, geom.L2, exp(1i*(thtL+phi)) ) );
    trajL = bsxfun(@plus, geom.center, trajL);

    % use late theta trajectory when we are not using the primary strategy
    traj(useL,:) = trajL(useL,:);
    
    % DT's are # of steps to push back the start of this component
    % of the trajectory.
    thtDT = zeros(size(Tht));
    phiDT = nmax - gt2.nphi;
    thtDTL = max(nmax - gt2.nthtN,0);

    output = packstruct(traj, trajL, useP, useL, tht, phi, thtDT, ...
                        thtDTL, phiDT, ltdiff);
    
    return