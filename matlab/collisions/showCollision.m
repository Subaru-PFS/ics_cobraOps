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
        CollAxis = exp(i*angle(FIBv(2) - FIBp(2)));
        vperp    = diff(FIBp)./CollAxis;
        vvictim  = diff(FIBv)./CollAxis;
        vmax = max(abs([vperp vvictim]));
       
        distance = abs(FIBv(2) - FIBp(2));
        
        plotcircle(-distance,0,1,'k:'); %perp fiber
        hold on; axis equal;
        plotcircle(0,0,1,'k:'); % victim fiber
        plot([i0 vperp/vmax]-distance,'r'); % perp velocity
        ph(2) = plot([i0 vvictim/vmax],'r'); % victim velocity
        
        ph(1) = plot(([geom.center(victim) ELBv(2) FIBv(2)]-FIBv(2))/CollAxis,'b'); % victim arms
        cmplx(@plotcircle, (ELBv(2) - FIBv(2)) / CollAxis, 1, 'k:'); % victim elbow
        plot(([geom.center(perp)   ELBp(2) FIBp(2)]-FIBp(2))/CollAxis - distance,'b'); % perp arms
 
    end
    
    legend(ph,'cobra arms','velocity vectors');
    title(sprintf('pid %d vs. pid %d',perp,victim));
    hold off;

    if collisionType(1) < 4 & collisionType(1) > 1
        disp('not a fiber-fiber collision and not a clean fiber-elbow collision')
        return
    end
    
