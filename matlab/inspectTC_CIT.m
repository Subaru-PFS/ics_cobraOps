function output = inspectTC_CIT(pid,minsteps,targets)
% output = inspectTC(pid, [minsteps=1], [targets=1:length(data)]);
% pid: positioner ID of interest
% minsteps: ignore convergences (within 10 microns) taking fewer steps than minsteps
%           defaults to 1
%
% output: for now, it's an array of the number of steps per convergence.

% cal : structure containing theta0 [deg], center [pix], arm1 [pix],
%       arm2 [pix] for all motors (use loadCfgXml)
    fignum = 876;
    cal = loadCfgXml;
    data0 = loadmats;
    % move all the data structures with names ending in _str into a
    % numbered cell array for easy access
    kk = 0;
    for ff=sort(fields(data0)')
        if ~isempty(regexp(ff{1},'_str$'))
            kk = kk+1;
            pdata{kk} = data0.(ff{1});
            pids(kk) = pdata{kk}(1).pid;
        end
    end

    data = pdata{find(pid == pids)};

    if ~exist('minsteps','var'), minsteps = 1; end;
    if ~exist('targets','var'), targets=1:length(data); end;

    
    % Blame Johannes for the crazy specification of taking the fifth
    %  (ammendment?) -- no -- fifth element in the J2_t (joint 2
    %  target) vector.
    for jj = 1:length(data)
        data(jj).J1_t = data(jj).J1_t(5);
        data(jj).J2_t = data(jj).J2_t(5);
    end


    mid = 1; % hard coded to deal with new getARMval

    maxiter = length(data(1).J1);
    
    DtoR = pi/180;
    pID = data.pid;
    
    pix2um = getARMval(cal,pID,mid,'Pixel_scale'); % pixels
    theta0 = getARMval(cal,pID,mid,'CW_Global_base_ori_z') * DtoR;
    cobracenter = getARMval(cal,pID,mid,'Global_base_pos_x') + ...
        getARMval(cal,pID,mid,'Global_base_pos_y') * i;
    arm1 = getARMval(cal,pID,mid,'Link1_Link_Length') * pix2um;
    arm2 = getARMval(cal,pID,mid,'Link2_Link_Length') * pix2um;
    RR = @(phi) sqrt(cos(phi) * 2 * arm1 * arm2 + arm1^2 + arm2^2);
    J1err = [data.J1err];
    Rtgt = RR(pi - [data.J2_t]);
    RdJ1 = bsxfun(@times, J1err, Rtgt);
    RdJ2 = arm2 * [data.J2err];
    J2_t = [data.J2_t];
    RdJ1J2 = RdJ1 + i*RdJ2;

    for kk=1:length(targets)
        jj = targets(kk);
        e5  = CobraTPerror( 5,pi - data(jj).J2_t,arm1,arm2);
        e10 = CobraTPerror(10,pi - data(jj).J2_t,arm1,arm2);

        in5  = inpoly( e5, RdJ1J2(:,jj));
        in10 = inpoly(e10, RdJ1J2(:,jj));
        try
            nsteps(jj) = find(in10,1);
        catch
            nsteps(jj) = 16;
        end;
        if nsteps(jj) > minsteps
            figure(pid)
            clf;
            subplot(221)
            dataok = find(~isnan(data(jj).curPos) & data(jj).curPos ~= 0);
            % outer patrol limit
            plotcircle(0,0,arm1+arm2,'r');            hold on;
            % 10, 5 micron target circles
            cmplx(@plotcircle,(data(jj).targCmplx - cobracenter)*pix2um,10,'r');
            cmplx(@plotcircle,(data(jj).targCmplx - cobracenter)*pix2um, 5,'m');
            % hard stop angles for theta
            plot([0 (arm1+arm2)*exp(i*theta0)], 'k-.');%os
            plot([0 (arm1+arm2)*exp(i*theta1)], 'k-'); %ss
            % cobra arms at target
            plot([0                   arm1 * exp(i*THTtgt(jj)) ...
                  (data(jj).targCmplx - cobracenter)*pix2um],'ro--');
            % cobra arms in initial position
            XYinit = (data(jj).curPos(dataok(1)) - cobracenter)*pix2um;
            THTinit = angle(XYinit) + acos((XYinit*XYinit' + arm1.^2 - arm2.^2)./(2*arm1.*abs(XYinit)));
            plot([0, arm1 * exp(i*THTinit), XYinit],'go-');
            % cobra arms in final position
            XYlast = (data(jj).curPos(dataok(end)) - cobracenter)*pix2um;
            THTlast = angle(XYlast) + acos((XYlast*XYlast' + arm1.^2 - arm2.^2)./(2*arm1.*abs(XYlast)));
            plot([0, arm1 * exp(i*THTlast), XYlast],'bo-');
            % fiber trajectory from XY position
            plot((data(jj).curPos(dataok) - cobracenter)*pix2um,'k.:');
            % target location
            % plot((data(jj).targCmplx - cobracenter)*pix2um,'rx','MarkerSize',20);
            
            xlabel('X [um]');
            ylabel('Y [um]');
            title(sprintf('%d: XY coordinates',jj-1));
            hold off;
            
            subplot(223)
            plot(RdJ1J2(3:end,jj),'o:');axis equal; hold on;
            cmplx(@scatter,RdJ1J2(3:end,jj),30,(3:maxiter)','fill');
            try, caxis([3 find(data(jj).status,1)]); end;
            colormap('gray');
            cmplx(@text,RdJ1J2(3:end,jj)+1,num2cell(3:maxiter));
            plot(e10,'r');
            plot(e5,'m');
            xlabel('R\Delta\theta [um]');
            ylabel('r_2\Delta\phi [um]');
            title(sprintf('%d in \\Delta\\theta-\\Delta\\phi space',jj-1));
            %% labels of the target theta, phi angles.  not necessary, given 221 figure
% $$$             try
% $$$                 textbp(sprintf('\\theta_t = %.2f rad (%.0f deg)\n\\phi_t = %.2f rad(%.0f deg)',...
% $$$                                data(jj).J1_t, 180/pi*data(jj).J1_t, ...
% $$$                                data(jj).J2_t, 180/pi*data(jj).J2_t));
% $$$             catch
% $$$                 text(.2,.2,sprintf('\\theta_t = %.2f rad (%.0f deg)\n\\phi_t = %.2f rad(%.0f deg)',...
% $$$                                    data(jj).J1_t, 180/pi*data(jj).J1_t, ...
% $$$                                    data(jj).J2_t, 180/pi*data(jj).J2_t),'units','normal');
% $$$             end
            hold off;

            set(0,'DefaultLineMarkerSize',8);
            subplot(222)
            j1pos = find(J1err(:,jj) > 0);
            j2pos = find(data(jj).J2err > 0);
            j1neg = find(J1err(:,jj) < 0);
            j2neg = find(data(jj).J2err < 0);
            plot(j1pos,  J1err(j1pos,jj), 'b^'); hold on;
            plot(j1neg, -J1err(j1neg,jj), 'bv', 'MarkerFaceColor','none');
            plot(j2pos,  data(jj).J2err(j2pos), 'r^');
            plot(j2neg, -data(jj).J2err(j2neg), 'rv', 'MarkerFaceColor','none');
            plot(abs(J1err(:,jj)),'b-');
            plot(abs(data(jj).J2err),'r-');
% $$$             plot(abs(data(jj).J1_S)/1000,'k--');
            hold off;
            set(gca,'yscale','log');
            try, refline(find(in5,1),1,Inf,'r-'); end
            try, refline(find(in10,1),1,Inf,'k--'); end
            try, refline(find(data(jj).status,1),1,Inf,'g--'); end
            hold off;
            ylabel('RADIANS');
            xlabel('iteration #');
            title('angular error vs. iteration #')
            
% $$$       %% MICRONS vs ITERATION plot
% $$$       j1pos = find(RdJ1(:,jj) > 0);
% $$$       j2pos = find(RdJ2(:,jj) > 0);
% $$$       j1neg = find(RdJ1(:,jj) < 0);
% $$$       j2neg = find(RdJ2(:,jj) < 0);
% $$$       plot(j1pos,  RdJ1(j1pos,jj), '^'); hold on;
% $$$       plot(j1neg, -RdJ1(j1neg,jj), 'v', 'MarkerFaceColor','none');
% $$$       plot(j2pos,  RdJ2(j2pos,jj), 'r^');
% $$$       plot(j2neg, -RdJ2(j2neg,jj), 'rv', 'MarkerFaceColor','none');
% $$$       plot(abs(RdJ1(:,jj)),'-');
% $$$       plot(abs(RdJ2(:,jj)),'r-');
% $$$       hold off;
% $$$       set(gca,'yscale','log');
% $$$       try, refline(find(in5,1),1,Inf,'r-'); end
% $$$       try, refline(find(in10,1),1,Inf,'k--'); end
% $$$       try, refline(find(data(jj).status,1),1,Inf,'g--'); end
% $$$       hold off;
% $$$       ylabel('MICRONS');
% $$$       xlabel('iteration #');
% $$$       title('distance error vs. iteration #')
            
            subplot(224)
            %% STEPS vs ITERATIONS
            j1pos = find(data(jj).J1_S > 0);
            j2pos = find(data(jj).J2_S > 0);
            j1neg = find(data(jj).J1_S < 0);
            j2neg = find(data(jj).J2_S < 0);
            j1zero = find(data(jj).J1_S == 0);
            j2zero = find(data(jj).J2_S == 0);
            try, hh(1) = plot(j1pos,  data(jj).J1_S(j1pos), 'b^'); hold on;
            catch hh(1) = plot(0,NaN,'b^'); hold on; end
            try, hh(3) = plot(j1neg, -data(jj).J1_S(j1neg), 'bv');
            catch hh(3) = plot(0,NaN,'bv'); end
            try, hh(4) = plot(j2pos,  data(jj).J2_S(j2pos), 'r^');
            catch hh(4) = plot(0,NaN,'r^'); end
            try, hh(6) = plot(j2neg, -data(jj).J2_S(j2neg), 'rv');
            catch, hh(6) = plot(0,NaN,'rv'); end
            plot(abs(data(jj).J1_S),'b-');
            plot(abs(data(jj).J2_S),'r-');
            % set zeros to 0.1 so they show up on log scale.
            refline(0,   0.5, 0, 'k:');
            try, hh(2) = plot(j1zero, 0.5*ones(size(j1zero)),'+'); 
            catch, hh(2) = plot(0,NaN, '+');end
            try, hh(5) = plot(j2zero, 0.5*ones(size(j2zero)),'rx');
            catch, hh(5) = plot(0,NaN, 'rx');end
            hold off;
            set(gca,'yscale','log');
            try, refline(find(data(jj).status,1),1,Inf,'g-'); end
            xlabel('iteration #');
            ylabel('steps commanded');
            title('steps vs. iteration #');
            curylim = ylim;
            ylim([0.1 curylim(2)]);
            legend(hh,...
                   '\Delta\theta > 0',...
                   '\Delta\theta = 0',...
                   '\Delta\theta < 0',...
                   '\Delta\phi > 0',...
                   '\Delta\phi = 0',...
                   '\Delta\phi < 0');

            set(0,'DefaultLineMarkerSize','factory');
            drawnow;

% $$$             figure(pid+1e4)
% $$$             clf;
% $$$             % error to move ratio for moves outside of 10 um
% $$$             subplot(211);
% $$$             plot(rdiff(data(jj).J1err).*(~in10(1:end-1)'),'bo-'); hold on;
% $$$             plot(rdiff(data(jj).J2err).*(~in10(1:end-1)'),'ro-'); hold off;
% $$$             try, refline(find(in10,1),1,Inf,'k--'); end
% $$$             subplot(212);
% $$$             drawnow;
            
            if kk < length(targets)
                input('next..');
            end
        end
    end
    output = nsteps;
    
end

function output=CobraTPerror(maxerr, J2_t, arm1, arm2)
% calculate a locus of points in dtheta-dphi space corresponding to
% an error circle of radius maxerr in XY space
%
% maxerr: radius of error circle in XY space in MICRONS
% J2_t: target joint 2 angle in RADIANS
% arm1: J1 arm length in MICRONS, defaults to 2375
% arm2: J2 arm length in MICRONS, defaults to 2375

    if ~exist('arm1','var'), arm1 = 2375; end
    if ~exist('arm2','var'), arm2 = 2375; end
    
    if J2_t > pi
        warning('J2_t should be in radians, not degrees');
    end
    
    RR = sqrt(cos(J2_t) * 2 * arm1 * arm2 + arm1^2 + arm2^2);
    
    % for simplicity, assume target lies on the Y axis (X = 0)
    xyerror = i*RR + maxerr * exp(i*2*pi*(0:.01:1));

    J1_t = pi/2 - asin(arm2 * sin(J2_t) / RR);
    
    % force error bounds to be within the patrol region
    xyR = min(abs(xyerror), arm1+arm2);
    xyT = angle(xyerror);
    xyerror = xyR .* exp(i*xyT);

    PHI = acos((xyR.^2 - arm1^2 -arm2^2) / (2 * arm1 * arm2));
    THT = xyT - asin(arm2 * sin(PHI) ./ xyR);
    
    output = RR*(THT - J1_t) + i*arm2*(PHI - J2_t);
    output = output(:);
end
