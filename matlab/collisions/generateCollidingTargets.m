function output = generateCollidingTargets(nfields,density)

    ncobras = 13;
    toggle.include_nearest_neighbors = false;
    
    targets = zeros(ncobras,nfields);
    phiDT = zeros(ncobras,nfields);
    thtDT = zeros(ncobras,nfields);
    useP = logical(zeros(ncobras,nfields));
    bench = defineBenchGeometry([],1,1);

    ctr.loops = 0;
    
    for ff = 1:nfields
        tic;
        assigned = false(ncobras,1);
        while true
            ctr.loops = ctr.loops + 1;
            
            res = simFun(density,'',1,1,'UseThisBench',bench);
            
            nn = find(res.minDist < 2.5 * res(1).bench.minDist/2);
            
            assign_these = false(ncobras,1); % track assignments and update
                                             % after the for loop over pairs.
            for pair = nn'
                cobras = [res.bench.NN.row(pair) res.bench.NN.col(pair)];
                if toggle.include_nearest_neighbors
                    [rr cc] = find(res.bench.nnMap(cobras,:));
                    cobras = unique(cc);
                end
                if ~(sum(assigned(cobras))) % none of the cobras of interest are assigned
                    assign_these(cobras) = true;
                end
            end
            targets(assign_these,ff)  = res.targets(assign_these);
            phiDT(assign_these,ff) = res.Traj.phiDT(assign_these);
            thtDT(assign_these,ff) = (res.Traj.thtDTL(assign_these) .* res.Traj.useL(assign_these) +...
                                      res.Traj.thtDT(assign_these) .* ~res.Traj.useL(assign_these));
            useP(assign_these,ff) = res.Traj.useP(assign_these);
            assigned = assigned | assign_these;
            
            %% test for exiting the loop
            maxrun = 0;
            run = 0; % set to 1 if nn's are included
            runarr = [];
            for jj = 1:ncobras
                % increment run when there are consecutive unassigned cobras
                % reset to zero is there is an assigned one
                % to pad the end condition, add (jj==ncobras) to the first term
                run = (run + ~assigned(jj)) * ~assigned(jj);
                runarr(jj) = run;
                maxrun = max(runarr);
                if (maxrun >= 2), break, end;
            end
            if (maxrun < 2), break, end;
        end 
    
        figure(99);
        imagesc(assigned');
        title(sprintf('field %d',ff));
        axis equal;
        axis tight;
        drawnow;

        % assign the rest of the targets
        targets(~assigned,ff)  = res.targets(~assigned);
        phiDT(~assigned,ff) = res.Traj.phiDT(~assigned);
        % leave thtDT at zero for unassigned cobras.
        useP(~assigned,ff) = res.Traj.useP(~assigned);
        fprintf(1,'Field %d done in %f seconds\n', ff, toc);
    end
    %% Create the files... 
    %% The zeros have to be waitsteps. 
    pids = bench.pids;
    
    direction = {'N','P'};
    
    for jj = 1:length(bench.pids)
        fname = sprintf('TargetList_mId_%d_pId_%d.txt',bench.mids(jj),bench.pids(jj));
        tfile = fopen(fname,'a');

        %    targetoutput = [targets, zeros(length(targets)), zeros(length(targets))]; 
% $$$         plot(targets(jj,:),'b.');
% $$$         plot(bench.center(jj),'rx');
        for qq = 1:length(targets(jj,:))
            fprintf(tfile,'%7.2f,%7.2f,%5d,%5d,%d\n',real(targets(jj,qq)),imag(targets(jj,qq)),...
                    round(thtDT(jj,qq)), round(phiDT(jj,qq)), useP(jj,qq)+1);
        end
        %        fprintf(1,'wrote %s\n',fname)
        fclose(tfile);
    end
    
    output = packstruct(targets, thtDT, phiDT, pids, ctr);