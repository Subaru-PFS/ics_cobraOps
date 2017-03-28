clear all;
close all;

%% above code commented out because that functionality has been
%% pushed into loadCfgXml.
baseCfg = loadCfgXml;

% Start new config
newConfig = baseCfg;
cw = 1;
data = processMsimCenters(baseCfg);
makeMaps = 0;
cli_answer = input(sprintf('Do you want to create and store motor maps too? [y|N] '),'s');
if strcmp(cli_answer,'y')
    makeMaps = 0;
end
%keyboard;
for kk = 1:size(data, 2)
    
    pid = data(kk).pid;
    center = data(kk).j1center;
    thetaCentroidsCM = data(kk).j1centroids(data(kk).j1ok) - data(kk).j1center;
    thetaRadii = abs(thetaCentroidsCM);
    thetaRadius = mean(thetaRadii);
    data(kk).cold = getARMval(baseCfg,pid,1,'Global_base_pos_x') + 1i * getARMval(baseCfg,pid,1,'Global_base_pos_y');
    data(kk).l1old = getARMval(baseCfg,pid,1,'Link1_Link_Length');
    data(kk).l2old = getARMval(baseCfg,pid,1,'Link2_Link_Length');
    data(kk).ccwold = getARMval(baseCfg,pid,1,'CCW_Global_base_ori_z');
    data(kk).cwold = getARMval(baseCfg,pid,1,'CW_Global_base_ori_z');
    data(kk).amaxOld = getARMval(baseCfg,pid,1,'phiMax');
    data(kk).aminOld = getARMval(baseCfg,pid,1,'phiMin');
    data(kk).ot1fw = getARMval(baseCfg,pid,1,'Link1_fwd_Duration');
    data(kk).ot1rv = getARMval(baseCfg,pid,1,'Link1_rev_Duration');
    data(kk).ot2fw = getARMval(baseCfg,pid,1,'Link2_fwd_Duration');
    data(kk).ot2rv = getARMval(baseCfg,pid,1,'Link2_rev_Duration');
    
    
    data(kk).link1 = abs(data(kk).j1center - data(kk).j2center);
    j2centroidsCM = data(kk).j2centroids(data(kk).j2ok) - data(kk).j2center;
    phiRadii = abs(j2centroidsCM);
    data(kk).link2 = mean(phiRadii);
    stdRad = std(phiRadii);
    
    
    %% HERE CALCULATE PHI MIN and MAX angle
    a1 = angle(data(kk).j2center - data(kk).j1center);
    phiAngles = unwrap(angle(j2centroidsCM));
    %keyboard;
    a2 = min(phiAngles);
    a3 = max(phiAngles);
    
    % Something like this? yes, but not in degrees. you don't want
    % crazy numbers, so cut the angle at pi/2 (where phi should never
    % go) and then mod the result by 2*pi to make it contiguous.
    % With this setup, the unwrap in the defnition of phiAngles is
    % unnecessary.
    %
    % If you want the zero angle direction to be the folded
    % configuration, then change the sign of a1.  you may also have to
    % change the sign of the cut and the resulting angle.  I don't
    % recommend defining angles that way.  It's VERY UNNATURAL.  Euler
    % would NOT APPROVE!
    cut = pi/2;
    data(kk).phiIn = mod((a2 - a1) - cut, 2*pi) + cut; % - 2*pi;
    data(kk).phiOut = mod((a3 - a1) - cut, 2*pi) + cut; % - 2*pi;
    data(kk).amin = data(kk).phiIn * 180/ pi -180;
    data(kk).amax = data(kk).phiOut *180/pi -180 ;
    
    % Getting the theta phi angles for the data.
    qphi = XY2TP(data(kk).j2centroids - data(kk).j1center, data(kk).link1, data(kk).link2);
    qphi.tht = mod(qphi.tht, 2*pi);
    data(kk).j2phi = mod(real(qphi.phi),pi);
    data(kk).j2tht = qphi.tht;
    
    %Calculate neg. hardstop CCW
    thtm = nanmean(qphi.tht(5:30));
    thtstd = nanstd(qphi.tht(5:30));
    thtall = qphi.tht(qphi.tht>thtm-thtstd*2.5 & qphi.tht<thtm+thtstd*2.5);
    data(kk).thtm = nanmean(thtall);
    
    % Getting the theta phi angles for the pos hardstop.
    qphi2 = XY2TP(data(kk).j2centroids2 - data(kk).j1center, data(kk).link1, data(kk).link2);
    data(kk).j2phi2 = mod(real(qphi2.phi),pi);
    data(kk).j2tht2 = qphi2.tht;
    
    try
    %Calculate pos. hardstop CW
    thtm2 = nanmean(qphi2.tht(5:30));
    thtstd2 = nanstd(qphi2.tht(5:30));
    thtall2 = qphi2.tht(qphi2.tht>thtm2-thtstd2*2.5 & qphi2.tht<thtm2+thtstd2*2.5);
    data(kk).thtm2 = nanmean(thtall2);
    
    qtht = XY2TP(data(kk).j1centroids - data(kk).j1center, data(kk).link1, data(kk).link2);
    data(kk).j1phi = real(qtht.phi);
    data(kk).j1tht = unwrap(qtht.tht - data(kk).thtm);
    catch
        disp 'No CW hardstop in test data';
        
    end
end


centersfig = figure(1)
set (gca,'Ydir','reverse')
for kk = 1:size(data, 2)
    figure(1)
    %%
    pid =data(kk).pid;
    title(sprintf('PID %d',pid));
    
    %   title(sprintf('PID %d\n\\Theta angle range = %f\\tau',pid,...
    %                 range(unwrap(angle(thetaCentroidsCM))/(2*pi))));
    axis([real(data(kk).cold)-80 real(data(kk).cold)+80 imag(data(kk).cold)-80 imag(data(kk).cold)+80]);
    daspect([1,1,1]);
    hold on;
    
    h = plot(data(kk).j2centroids, 'ko');
    % set(h, 'color', [0.5 0.5 0.5])
    plot(data(kk).j2centroids(data(kk).j2ok),'gx');
    
    h = plot(data(kk).j1centroids, 'ko');
    % set(h, 'color', [0.5 0.5 0.5])
    plot(data(kk).j1centroids(data(kk).j1ok),'gx');
    
    
    xlabel('X [pix]');
    ylabel('Y [pix]');
    
    plot(data(kk).cold, 'kx');
    cmplx(@circle, data(kk).cold, data(kk).l1old, 'k');
    cmplx(@circle, data(kk).cold, data(kk).l1old + data(kk).l2old, 'k');
    
    plot([data(kk).j2center, data(kk).j2center + data(kk).link2 * exp(1i * (data(kk).phiIn + data(kk).thtm))],'k');
    plot([data(kk).j2center, data(kk).j2center + data(kk).link2 * exp(1i * (data(kk).phiOut + data(kk).thtm))],'b');
    
    plot([data(kk).j1center, data(kk).j2center + data(kk).link2 * exp(1i * (data(kk).thtm ))],'b--');
    
    plot(data(kk).j1center, 'rx');
    plot(data(kk).j2center, 'rx');
    cmplx(@circle, data(kk).j1center, data(kk).link1, 'b');
    cmplx(@circle, data(kk).j2center, data(kk).link2, 'b');
    cmplx(@circle, data(kk).j1center, data(kk).link2+data(kk).link1, 'b');
    cmplx(@text, data(kk).j1center + 10 + 10i, sprintf('PID %d', pid));
    % centerConf = std(thetaRadii)/(thetaRadius
    legend(sprintf('Center new: (%.4f,%.4f), (old: (%.4f,%.4f))',real(data(kk).j1center),imag(data(kk).j1center), real(data(kk).cold),imag(data(kk).cold)),...
        sprintf('Center Confidence %.4f', std(thetaRadii)/(thetaRadius * sqrt(length(thetaRadii)))), ...
        sprintf('Link1: %.2f (old: %.2f)', data(kk).link1, data(kk).l1old), ...
        sprintf('Link2: %.2f (old: %.2f)', data(kk).link2, data(kk).l2old), ...
        sprintf('Phi Angle Start:   %.2f (old: %.2f)', data(kk).amin, data(kk).aminOld), ...
        sprintf('Phi Angle End: %.2f (old: %.2f)', data(kk).amax, data(kk).amaxOld), ...
        sprintf('CCW Hardstop: %.2f (old: %.2f)', mod(data(kk).thtm,2*pi)*180/pi, data(kk).ccwold));%,... % Center confidence= how far was phi out for measuring the center.
       % sprintf('CW Hardstop: %.2f (old: %.2f)', mod(data(kk).thtm2,2*pi)*180/pi, data(kk).cwold)); % Center confidence= how far was phi out for measuring the center.
    
    
    
    hold off;
    
    
    if 0
        %%suggestions
        % use unwrap
        % set fwd limit to 6.3 rad, rev limit to .5 rad
        % XX = reshape(unrwap(.j1tht - .thtm), 200, 10)
        % [r c] = find(XX > 6.5)
        % etc.
        % may want XFwd = XX(1:100,:), XRev = XX(101:end,:)
        
        
        figure(pid)
        hold on;
        %   plot(data(kk).j2tht, 'r.');
        plot(data(kk).j1tht, 'rx');
        plot((data(kk).j2phi), 'b.');
        %   plot(data(kk).j1phi, 'bx');
        legend(sprintf('thttest THT'), ...
            sprintf('phitest PHI'));
        %legend(sprintf('phitest THT'),...
        %       sprintf('thttest THT'), ...
        %       sprintf('phitest PHI'), ...
        %       sprintf('thttest PHI'));
    end
    
    
    cli_answer = input(sprintf('Use these data for PID %d? [Y|n] ',pid),'s');
    if strcmp(cli_answer,'n')
        disp('Old centers and link lengths used.');
    else
        newConfig = setARMval(newConfig,pid, 1, 'Global_base_pos_x',real(data(kk).j1center));
        newConfig = setARMval(newConfig,pid, 1, 'Global_base_pos_y',imag(data(kk).j1center));
        newConfig = setARMval(newConfig,pid, 1, 'Link1_Link_Length',data(kk).link1);
        newConfig = setARMval(newConfig,pid, 1, 'Link2_Link_Length',data(kk).link2);
        ccwangle = mod(data(kk).thtm,2*pi)*180/pi;
        newConfig = setARMval(newConfig,pid,1,'CCW_Global_base_ori_z', ccwangle);
%        cwangle = mod(data(kk).thtm2,2*pi)*180/pi;
%        newConfig = setARMval(newConfig,pid,1,'CW_Global_base_ori_z', cwangle);

        %if(amin < 0) amin = 0.1; end % No neg angles for MSIM please.
        newConfig = setARMval(newConfig,pid, 1, 'Joint2_CCW_limit_angle',data(kk).amin);
        newConfig = setARMval(newConfig,pid, 1, 'Joint2_CW_limit_angle',data(kk).amax);
        disp(sprintf('Data for PID %d successfully updated!',pid));
        %   else
        %      disp('Please type y or n! Result not used.');
    end
end

%% Create MotorMaps.
if makeMaps
    for kk = 1:size(data, 2)
        
        s1.startAngle = data(kk).j1tht(1:end-1).*180/pi;
        s1.finishAngle = data(kk).j1tht(2:end).*180/pi;
        s1.moveSizes = diff(data(kk).j1tht).*180/pi;
        s1.mmap = s1.moveSizes ./60;
        
        reverseVals = s1.mmap < 0;
        s1f.startAngle = s1.startAngle(reverseVals ==0);
        s1f.finishAngle = s1.finishAngle(reverseVals ==0);
        s1f.moveSizes = s1.moveSizes(reverseVals ==0);
        s1f.mmap = s1.mmap(reverseVals ==0);
        
        s1r.startAngle = s1.startAngle(reverseVals);
        s1r.finishAngle = s1.finishAngle(reverseVals);
        s1r.moveSizes = s1.moveSizes(reverseVals);
        s1r.mmap = s1.mmap(reverseVals) .*-1;
        
        s2.startAngle = data(kk).j2phi(1:end-1).*180/pi;
        s2.finishAngle = data(kk).j2phi(2:end).*180/pi;
        s2.moveSizes = diff(data(kk).j2phi).*180/pi;
        s2.mmap = s2.moveSizes ./60;
        
        reverseVals = s2.mmap < 0;
        s2f.startAngle = s2.startAngle(reverseVals ==0);
        s2f.finishAngle = s2.finishAngle(reverseVals ==0);
        s2f.moveSizes = s2.moveSizes(reverseVals ==0);
        s2f.mmap = s2.mmap(reverseVals ==0);
        
        s2r.startAngle = s2.startAngle(reverseVals);
        s2r.finishAngle = s2.finishAngle(reverseVals);
        s2r.moveSizes = s2.moveSizes(reverseVals);
        s2r.mmap = s2.mmap(reverseVals).*-1;
        
        
        
        
        % Load previous motor maps from XML
        fwdMapS1 = getOldMap(1, 'fwd', data(kk).pid, baseCfg);
        rvsMapS1 = getOldMap(1, 'rev', data(kk).pid, baseCfg);
        fwdMapS2 = getOldMap(2, 'fwd', data(kk).pid, baseCfg);
        rvsMapS2 = getOldMap(2, 'rev', data(kk).pid, baseCfg);
        
        
        % Create Motor Maps from data
        fwdMapS1new = createMotorMap(s1f, fwdMapS1);
        rvsMapS1new = createMotorMap(s1r, rvsMapS1);
        fwdMapS2new = createMotorMap(s2f, fwdMapS2);
        rvsMapS2new = createMotorMap(s2r, rvsMapS2);
        
        % Plot results
        fh = figure('name', ['Motor Maps for Positioner ', num2str(data(kk).pid)] );
        clf;
        plotmmap('s1f', s1f, fwdMapS1new, fwdMapS1, 1, 1);
        plotmmap('s1r', s1r, rvsMapS1new, rvsMapS1, 1, 2);
        plotmmap('s2f', s2f, fwdMapS2new, fwdMapS2, 1, 3);
        plotmmap('s2r', s2r, rvsMapS2new, rvsMapS2, 1, 4);
        
        newSlow_J1_fwd_vmap = [length(fwdMapS1new) 200 fwdMapS1new(1,:);  length(fwdMapS1new) 200 fwdMapS1new(2,:)];
        newSlow_J1_rev_vmap = [length(rvsMapS1new) 200 rvsMapS1new(1,:);  length(rvsMapS1new) 200 rvsMapS1new(2,:)];
        newSlow_J2_fwd_vmap = [length(fwdMapS2new) 200 fwdMapS2new(1,:);  length(fwdMapS2new) 200 fwdMapS2new(2,:)];
        newSlow_J2_rev_vmap = [length(rvsMapS2new) 200 rvsMapS2new(1,:);  length(rvsMapS2new) 200 rvsMapS2new(2,:)];
        
        
        % Write new Map to XML
        newConfig = setMMap(newConfig,data(kk).pid,'SLOW',1,'fwd',newSlow_J1_fwd_vmap);
        newConfig = setMMap(newConfig,data(kk).pid,'SLOW',1,'rev',newSlow_J1_rev_vmap);
        newConfig = setMMap(newConfig,data(kk).pid,'SLOW',2,'fwd',newSlow_J2_fwd_vmap);
        newConfig = setMMap(newConfig,data(kk).pid,'SLOW',2,'rev',newSlow_J2_rev_vmap);
        
        newConfig = setMMap(newConfig,data(kk).pid,'FAST',1,'fwd',newSlow_J1_fwd_vmap);
        newConfig = setMMap(newConfig,data(kk).pid,'FAST',1,'rev',newSlow_J1_rev_vmap);
        newConfig = setMMap(newConfig,data(kk).pid,'FAST',2,'fwd',newSlow_J2_fwd_vmap);
        newConfig = setMMap(newConfig,data(kk).pid,'FAST',2,'rev',newSlow_J2_rev_vmap);
        
        
    end
end

figure(1)
axis([0 2000 0 2000]);


% Save new XML with adjusted maps
[xmlfile, xmlfilepath] = uiputfile('*.xml','Save CobraConfig XML file with new centers and link lenghts.');
cobraCfg2xml(newConfig,fullfile(xmlfilepath,xmlfile));
