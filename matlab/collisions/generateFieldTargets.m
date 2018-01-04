function output = generateFieldTargets(nfields,density)
% output = generateFieldTargets(nfields,density)
% Generates random targets of given density for cobra testing
% uses simFun with a predefined bench.
    
    bench = defineBenchGeometry([],1,1);

    ncobras = length(bench.pids);
    targets = zeros(ncobras,nfields);
    phiDT = zeros(ncobras,nfields);
    thtDT = zeros(ncobras,nfields);
    useP = logical(zeros(ncobras,nfields));

    for ff = 1:nfields
        % get the targets
        res = simFun(density,'',1,1,'UseThisBench',bench,'SkipTargetReplan',false);
            
        % extract info for the target list: x,y, tht delay, phi delay, direction
        targets(:,ff)  = res.targets;
        phiDT(:,ff) = res.Traj.phiDT;
        thtDT(:,ff) = (res.Traj.thtDTL .* res.Traj.useL +...
                       res.Traj.thtDT .* ~res.Traj.useL);
        useP(:,ff) = res.Traj.useP;
    end
    %% Create the files... 
    %% The zeros have to be waitsteps. 
    pids = bench.pids;
    
    for jj = 1:length(bench.pids)
        fname = sprintf('TargetList_mId_%d_pId_%d.txt',bench.mids(jj),bench.pids(jj));
        tfile = fopen(fname,'a');

        for qq = 1:length(targets(jj,:))
            fprintf(tfile,'%7.2f,%7.2f,%5d,%5d,%d\n',real(targets(jj,qq)),imag(targets(jj,qq)),...
                    round(thtDT(jj,qq)), round(phiDT(jj,qq)), useP(jj,qq)+1);
        end
        %        fprintf(1,'wrote %s\n',fname)
        fclose(tfile);
    end
    
    output = packstruct(targets, thtDT, phiDT, pids, bench);