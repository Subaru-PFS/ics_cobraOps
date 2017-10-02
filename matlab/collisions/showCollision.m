function showCollision(simFunOut, pid1, pid2)
% show collision parameters graphically.
% if 2 ID's given, then interpret as pids
% if 1 ID given, interpret as a nearest-neighbor index.

    Coll = simFunOut.Coll;
    Traj = simFunOut.Traj;
    geom = simFunOut.bench;
    i0   = complex(0);

    if nargin == 2 % pid1 is really a NN index
        nn12 = pid1;
        pid1 = Coll.row(nn12);
        pid2 = Coll.col(nn12);
    end    
    %  first get the nearest-neighbor indices for the pid pair
    nn12 = Coll.rcindx(pid1,pid2);
    nn21 = Coll.rcindx(pid2,pid1);

    % find where collisions are happening
    isColliding = Coll.detected(nn12,:) | simFunOut.Coll.detected(nn21,:);
    if isempty(find(isColliding))
        disp('no collsions');
        return;
    end
    one_vec = ones(size(find(isColliding)));

    % list the types of collisions
    full([Coll.type(nn12,isColliding); ...
          Coll.type(nn21,isColliding)])
    
    %% from here out, we're only going to look at the first instance of interference.
    collisionIndx = find(isColliding,1); % the index of the first interference
    
    collisionType = (Coll.type(nn12,collisionIndx) + Coll.type(nn21,collisionIndx));
    % for each, 1: F-E, 3: F-F
    
    % default victim/perp identification
    victim = pid2;
    perp   = pid1;

    if Coll.type(nn21,collisionIndx) == 1 % 2 hits 1's elbow
        victim = pid1;
        perp   = pid2;
    end

    % need to calculate elbow position from trajectory and geometry
    toi = collisionIndx-1:collisionIndx; %"time of interest"
    perpTP   = XY2TP(Traj.traj(perp,toi)   - geom.center(perp),   geom.L1(perp),   geom.L2(perp));
    victimTP = XY2TP(Traj.traj(victim,toi) - geom.center(victim), geom.L1(victim), geom.L2(victim));

        
    ELBv = exp(i*victimTP.tht) * geom.L1(victim) + geom.center(victim);
    ELBp = exp(i*  perpTP.tht) * geom.L1(perp  ) + geom.center(perp);
    FIBv = Traj.traj(victim,toi);
    FIBp = Traj.traj(perp  ,toi);
    
    if collisionType == 1 % clean fiber-elbow hit
        CollAxis = exp(i*angle(ELBv(2) - FIBp(2)));
        vperp    = diff(FIBp)./CollAxis;
        vvictim  = diff(ELBv)./CollAxis;
        vmax = max(abs([vperp vvictim]));
       
        distance = abs(FIBp(2) - ELBv(2));

        plotcircle(-distance,0,1,'k:'); %perp fiber
        hold on; axis equal;
        plotcircle(0,0,1,'k:'); % victim elbow
        plot([i0 vperp/vmax]-distance,'r'); % perp velocity
        ph(2) = plot([i0 vvictim/vmax],'r'); % victim velocity
        
        ph(1) = plot(([geom.center(victim) ELBv(2) FIBv(2)]-ELBv(2))/CollAxis,'b'); % victim arms
        cmplx(@plotcircle, (FIBv(2) - ELBv(2)) / CollAxis, 1, 'k:'); % victim fiber
        plot(([geom.center(perp)   ELBp(2) FIBp(2)]-FIBp(2))/CollAxis - distance,'b'); % perp arms
    end
    
    if collisionType >= 4 % fiber-fiber collisions

        % calculated the interpolated position of impact by linear interpolation
        distance = abs(FIBv - FIBp); % 1x2 distances
% $$$         t_impact = (2 - distance(1))/diff(distance);
        
        % calculate interpolated position by fitting a quadratic
        tt = [0 .5 1];
        zv = FIBv(1) + diff(FIBv) * tt;
        zp = FIBp(1) + diff(FIBp) * tt;
        dsquared = abs(zv-zp).^2;
        t_impact = roots(polyfit(tt,dsquared - 4, 2));
        t_impact = t_impact(t_impact > 0 & t_impact <= 1);

        % velocities and positions in focal plane coords
        v_vict   = diff(FIBv);
        z_vict   = FIBv(1) + t_impact * v_vict;
        v_perp   = diff(FIBp);
        z_perp   = FIBp(1) + t_impact * v_perp;
        vmax     = max(abs([v_perp v_vict]));

        % find the elbow locations in the interpolated position
        tp_vict  = XY2TP(z_vict - geom.center(victim), geom.L1(victim), geom.L2(victim));
        elb_vict = geom.center(victim) + geom.L1(victim) * exp(i*tp_vict.tht);
        tp_perp  = XY2TP(z_perp - geom.center(perp), geom.L1(perp), geom.L2(perp));
        elb_perp = geom.center(perp) + geom.L1(perp) * exp(i*tp_perp.tht);

        % put velocities in collision axis coordinates
        CollAxis = exp(i*angle(z_vict - z_perp));
        v_perp    = v_perp./CollAxis;
        v_vict    = v_vict./CollAxis;

        % generate the figure
        plotcircle(-2,0,1,'k:'); % perp fiber
        hold on; axis equal;
        plotcircle( 0,0,1,'k:'); % victim fiber
        plot([i0 v_perp/vmax]-2,'r'); % perp velocity
        ph(2) = plot([i0 v_vict/vmax] ,'r'); % victim velocity
        
        ph(1) = plot(([geom.center(victim) elb_vict z_vict]-z_vict)/CollAxis,'b'); % victim arms
        plot(([geom.center(perp)   elb_perp z_perp]-z_perp)/CollAxis - 2,'b'); % perp arms

        cmplx(@plotcircle, (elb_vict-z_vict) / CollAxis, 1, 'k:'); % victim elbow
 
    end
    
    legend(ph,'cobra arms','velocity vectors');
    title(sprintf('pid %d vs. pid %d',perp,victim));
    hold off;

    if collisionType(1) < 4 & collisionType(1) > 1
        disp('not a fiber-fiber collision and not a clean fiber-elbow collision')
        return
    end
    
