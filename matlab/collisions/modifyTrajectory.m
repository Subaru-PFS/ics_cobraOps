function output=modifyTrajectory(Traj,preTraj,geom,phiDT,thtDT)
% 
% modify the theta and phi timings of the specified trajectories
%

    useP = Traj.useP;
    useN = ~useP;
    ltht = useP .* preTraj.lthtP + useN .* preTraj.lthtN;
    lphi = preTraj.lphiP;
    lmax = max([ltht lphi]);  % max # time bins needed

% $$$     Tht       = preTraj.ThtP;
% $$$     Tht(useN) = preTraj.ThtN(useN);
% $$$     Phi       = preTraj.PhiP;

    %% phi modifications -- shift to earlier time

    phiIndx = find(phiDT);
    for jj = 1:length(phiIndx)
        thisTrj = phiIndx(jj);
        dt = min(phiDT(thisTrj), lmax - lphi(thisTrj));
% $$$         figure(thisTrj);
% $$$         plot(Traj.phi(thisTrj,:)); hold on;
        Traj.phi(thisTrj,:) = padarray(Traj.phi(thisTrj, (dt+1):end), [0 dt],'replicate','post');
        Traj.phiDT(thisTrj) = dt;
% $$$         plot(Traj.phi(thisTrj,:)); hold off;

        
    end

    tht = Traj.tht;
    phi = Traj.phi;
    
    traj = (bsxfun(@times, geom.L1, exp(1i*tht)) + ...
            bsxfun(@times, geom.L2, exp(1i*(tht+phi)) ) );
    traj = bsxfun(@plus, geom.center, traj);

    Traj.traj = traj;
    
    output = Traj;