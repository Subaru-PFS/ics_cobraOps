function output=showCollision(simFunOut, pid1, pid2)
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

% $$$     % list the types of collisions
% $$$     full([Coll.type(nn12,isColliding); ...
% $$$           Coll.type(nn21,isColliding)])
    
    %% from here out, we're only going to look at the first instance of interference.
    collisionIndx = find(isColliding,1); % the index of the first interference
    
    collisionType = (Coll.type(nn12,collisionIndx) + Coll.type(nn21,collisionIndx));
    % for each, 1: F-E, 3: F-F
    
    if collisionType(1) < 4 & collisionType(1) > 1
        disp('not a fiber-fiber collision and not a clean fiber-elbow collision')
        return
    end
    
    if collisionType(1) == 4 | collisionType(1) == 5
        disp('oddball')
        keyboard;
    end
    
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
    
    if collisionType == 1
        VICT = ELBv; % part that gets hit
    elseif collisionType > 4
        VICT = FIBv;
    end
    
    distance = abs(VICT - FIBp); % 1x2 distances
        
    % calculate interpolated position by fitting a quadratic
    tt = [0 .5 1];
    zv = VICT(1) + diff(VICT) * tt;
    zp = FIBp(1) + diff(FIBp) * tt;
    dsquared = abs(zv-zp).^2;
    t_impact = roots(polyfit(tt,dsquared - 4, 2));
    t_impact = t_impact(t_impact > 0 & t_impact <= 1);

    % velocities and positions in focal plane coords
    v_vict   = diff(VICT);
    z_vict   = VICT(1) + t_impact * v_vict;
    v_perp   = diff(FIBp);
    z_perp   = FIBp(1) + t_impact * v_perp;
    vscale   = 5; % for plotting, boost velocity by this factor

    % elbow of the perpetrator
    tp_perp  = XY2TP(z_perp - geom.center(perp), geom.L1(perp), geom.L2(perp));
    elb_perp = geom.center(perp) + geom.L1(perp) * exp(i*tp_perp.tht);

    if collisionType == 1
        % find the non-interacting end locations in the interpolated position
        fib_vict = FIBv(1) + diff(FIBv) * t_impact; % approximate fiber location
        elb_vict = z_vict;
    elseif collisionType > 4
        % find the elbow position of the victim
        temp     = XY2TP(z_vict - geom.center(victim), geom.L1(victim), geom.L2(victim));
        fib_vict = z_vict;
        elb_vict = geom.center(victim) + geom.L1(victim) * exp(i*temp.tht);
        clear temp;
    end

    % put velocities in collision axis coordinates
    CollAxis = exp(i*angle(z_vict - z_perp));
    v_perp    = v_perp./CollAxis;
    v_vict    = v_vict./CollAxis;

    impact_force = (v_perp - v_vict)/2; % on the victim.  imag part is tangential CW on both parts
    v_final_normal = real(v_perp + v_vict)/2;

    % the normal velocities suffer inelastic collision, while the tangential velocities are
    % assumed to be unaffected by the collision (ie, coeff_of_friction = 0)
    v_perp_final = v_final_normal + i*imag(v_perp);
    v_vict_final = v_final_normal + i*imag(v_vict);

    % one may aslo want to modify this by taking into account the components of velocity
    % (momentum) due to the two motors.
    
    % velocity ratios: log(abs(initial/final))
    alrv_perp = log(abs(v_perp / v_perp_final));
    alrv_vict = log(abs((v_vict + eps)/ v_vict_final));

    % unit normals of the arms
    L2_perp = (elb_perp - z_perp)/CollAxis/geom.L2(perp);
    L1_perp = (geom.center(perp) - elb_perp)/CollAxis/geom.L1(perp);

    L2_vict = (elb_vict - fib_vict)/CollAxis/geom.L2(victim);
    L1_vict = (geom.center(victim) - elb_vict)/CollAxis/geom.L1(victim);
    
    %% generate the figure with the collision at the origin
    plotcircle(-1,0,1,'b:','linewidth',3); %perp fiber
    hold on; axis equal;
    plotcircle( 1,0,1,'g:','linewidth',3); % victim struck point

    plot([i0 v_perp*vscale]-1,'r');         % perp initial velocity
    ph(2) = plot(complex([0 v_vict*vscale]+1),'r'); % victim initial velocity

    plot([i0 v_perp_final*vscale],'b','linewidth',3);      % perp final velocity
    plot([i0 v_vict_final*vscale],'g','linewidth',3);      % victim final velocity

    ph(1) = plot(([geom.center(perp)   elb_perp z_perp]-z_perp)/CollAxis - 1,'k'); % perp arms
    plot(([geom.center(victim) elb_vict fib_vict]-z_vict)/CollAxis + 1,'k'); % victim arms

    cmplx(@plotcircle, (elb_perp-z_perp)/CollAxis-1, 1, 'b:','linewidth',3); % perp elbow
    if collisionType == 1
        cmplx(@plotcircle, (fib_vict - z_vict) / CollAxis + 1, 1, 'g:','linewidth',3); % victim fiber
    elseif collisionType > 4
        cmplx(@plotcircle, (elb_vict-z_vict) / CollAxis + 1, 1, 'g:','linewidth',3); % victim elbow 
    end
    
    legend(ph,'cobra arms','velocity \times 5');
    title(sprintf('pid %d vs. pid %d',perp,victim));
    
    % bump or jam?
    text(-3,0,sprintf('%.2f',alrv_perp),'fontsize',15,'horiz','center');
    text( 3,0,sprintf('%.2f',alrv_vict),'fontsize',15,'horiz','center');
    if abs(alrv_perp) > .5 & abs(alrv_vict) > .5
        plot(alrv_perp,max(alrv_vict,-3),'rx','markersize',10,'linewidth',3,'DisplayName','JAM!');
        collResult = 1;
    else
        plot(alrv_perp,max(alrv_vict,-3),'go','markersize',10,'linewidth',3,'DisplayName','bump...');
        collResult = 0;
    end        
    
    % elbow jam
    if collisionType == 1 & v_final_normal > 0.05
        if abs(angle(L1_vict/v_vict_final)) < pi/6
% $$$             plot(([geom.center(victim) elb_vict]-z_vict)/CollAxis + 1,'r--','linewidth',3,...
% $$$                  'DisplayName','Elbow Jam'); % victim L1
            collResult = collResult + 2;
        end
    end

    hold off;

    %% outputs

    output = collResult;
    
    % assume that the tangential component has no friction, so the collision is elastic along the axis,
    % but velocities are unaffected in the perpendicular axis
    %
    % 2vf = v1_init + v2_init, so the change in velocity for each component of the collision is
    % vf - v1_init = v2_init - vf, ie, equal and opposite.

    normal_force = real(impact_force);
    mu_fn = 0;
    tangential_force = imag(impact_force) * mu_fn * normal_force;

    
    % forces
    L2_perp_forces = -normal_force * L2_perp;

    if collisionType == 1
        L1_vict_forces = normal_force * L1_vict;
    elseif collisionType > 4
        L2_vict_forces = normal_force * L2_vict;
    end
    
            
