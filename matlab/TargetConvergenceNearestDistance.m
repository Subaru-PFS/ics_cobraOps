function output=TargetConvergenceNearestDistance()
% calculate the nearest fiber to arm distance for each set of
% nearest neighbor pairs

    % load the as-is bench from the config file
    bench = defineBenchGeometry([],1,1);
    
    % load all test data
    mats = loadmats; 
    flds = fields(mats);
    
    for jj = 1:length(flds)
        if ~isempty(strfind(flds{jj},'_str'))
            data = mats.(flds{jj});
            pid = data(1).pid;
            cobraIndx = find(bench.pids == pid);
            try
                pos(cobraIndx,:,:) = [data.curPos];
            catch
                keyboard
            end
            status(cobraIndx,:,:) = [data.status];
            errRat1(cobraIndx,:,:) = rdiff([data.J1err]);
            errRat2(cobraIndx,:,:) = rdiff([data.J2err]);
            target(cobraIndx,:) = [data.targCmplx];
            havedata(cobraIndx) = true;
        end
    end
    clear pid cobraIndx data;
    
    % turn have data from logical array into cobra index array
    havedata = find(havedata); 
    
    % cut down the arrays
    pids   = bench.pids(havedata);
    pos    = pos(havedata,:,:);
    status   = status(havedata,:,:);
    status(isnan(status)) = 0; 
    center = bench.center(havedata);
    L1     = bench.L1(havedata);
    L2     = bench.L2(havedata);
    errRat1 = errRat1(havedata,:,:);
    errRat2 = errRat2(havedata,:,:);
    target  = target(havedata,:);
 
    % find the pos first-index corresponding arrays
    [rr cc] = find(bench.nnMap(havedata,havedata));
    UpperTriangle = cc > rr;
    row = rr(UpperTriangle);
    col = cc(UpperTriangle);
    idx = sub2ind([length(havedata), length(havedata)], row, col);
    
    clear UpperTriangle havedata;
    
    Niterations = size(pos,2);
    Ntargets    = size(pos,3);
    Nneighbors  = length(rr);
    
    mm = 1./bench.rf; %conversion to mm from this coordinate system.
    
    for tt = 1:Ntargets
        for ii = 1:Niterations
            fibers = pos(:,ii,tt);
            angles = XY2TP(fibers - center, L1, L2);
            %% this is a fix for NaN thetas maybe should be moved into XY2TP, not sure 
            BADTHT = find(isnan(angles.tht));
            if ~isempty(BADTHT)
                replacement_angle = ( angle(fibers-center) + ...
                                      pi * ( L2 > L1 & abs(fibers-center) < L1 ) );
                angles.tht(BADTHT) = replacement_angle(BADTHT);
            end
            %%
            elbows = center + L1 .* exp(i*angles.tht);

            % colltype:: collision type (1 = tip to elbow, 2 = tip to arm, 3 = tip to tip)
            [distances, colltype] = pt2linesegment(fibers(rr), elbows(cc), fibers(cc)); 
            distances = distances * mm;
    
            dist_matrix = sparse(rr,cc,distances);
          % coll_sign   = sparse(rr,cc,(colltype < 3) * 2 - 1);
            mindist = full(min(dist_matrix, dist_matrix.'));
            
            dist(:,ii,tt) = mindist(idx);
       end
    end
    clear fibers angles BADTHT replacement_angle elbows distances colltype mindist dist_matrix;
    clear ii jj tt idx;
    
    labels = subarray(strsplit(sprintf('(%d, %d)_', [pids(row) pids(col)]'),'_')',...
                      '1:end-1');

    output = packstruct(dist,labels,row,col);
    output.pid1 = pids(row);
    output.pid2 = pids(col);



%     for ii = 1:150
%     imagesc(output.dist(:,:,ii)<2.1)
%     pause(0.5)
%     end
  
    cmapfignum = 5000;
    figure(cmapfignum)
    
    threshold = 2.2;
    
    A = squeeze(sum(output.dist < threshold, 2));
    B = abs(squeeze(pos(:,end-1,:)) - target) < (0.01/mm);
    %B = squeeze(sum(status, 2)<1);
    AB = zeros(size([A;B]));
    
    AB(1:2:end,:) = B;
    AB(2:2:end,:) = A;
    imagesc(AB)

    xlabel('Field #');
    ylabel('pIDs');
    set(gca,'YTick',2:2:Nneighbors,'YTickLabel',output.labels);
    clear A B;
    
    %% 3-D arrays are PID|PAIR, ITERATION, TARGET
    isjammed = logical(zeros(Nneighbors/2,Ntargets));
    
    for tgt = 1:Ntargets
        nearies = find(min(dist(:,:,tgt),[],2) < threshold)';
        for jj = nearies;
            close_shave = (dist(jj,1:end-1,tgt) < threshold);
            % row(jj) and col(jj) are the PID-index of the positioner
            % pids(row(jj)) and pids(col(jj)) are the PIDs of the positioner
            pidpair = pids([row(jj) col(jj)])';
            cobraR = row(jj);
            cobraC = col(jj);
           
            %keyboard;
            jammalammaR = AmIJammed(errRat1(cobraR,:,tgt), errRat2(cobraR,:,tgt), ...
                                    pos(cobraR,:,tgt) - center(cobraR),...
                                    target(cobraR,tgt) - center(cobraR),...
                                    L1(cobraR), L2(cobraR), mm);

            jammalammaC = AmIJammed(errRat1(cobraC,:,tgt), errRat2(cobraC,:,tgt), ...
                                    pos(cobraC,:,tgt) - center(cobraC),...
                                    target(cobraC,tgt) - center(cobraC),...
                                    L1(cobraC), L2(cobraC), mm);
            
            jammed = (jammalammaR | jammalammaC) & close_shave;
            if sum(jammed)
                 jammed
                isjammed(jj, tgt) = true;

                %% mark the collision map figure
                figure(cmapfignum); hold on;
                plot(tgt,2*jj,'rx','MarkerSize',20,'Linewidth',3);
                xlim([-4.5 4.5] + tgt);
                hold off;
                drawnow;
                
% $$$             end
% $$$             if false


% $$$                 %% inspectTC on jammed pair
% $$$                 for pp = pidpair
% $$$                     inspectTC_CIT(pp,1,tgt);
% $$$                 end

                %% the 2-cobra figure:
                figure(jj+8000)
                plot(pos(cobraR,:,tgt),'bo--'); hold on;
                plot(pos(cobraC,:,tgt),'go--'); 
                plot(target([cobraR cobraC],tgt), 'kx','MarkerSize',20,'Linewidth',3);
                plot(cobraArms(pos(cobraR,jammed,tgt), center(cobraR), L1(cobraR), L2(cobraR)),...
                     'ro-','MarkerSize',10);
                plot(cobraArms(pos(cobraC,jammed,tgt), center(cobraC), L1(cobraC), L2(cobraC)),...
                     'ro-','MarkerSize',10);
                cmplx(@plotcircle, center(cobraR), L1(cobraR)+L2(cobraR),'k:');
                cmplx(@plotcircle, center(cobraC), L1(cobraC)+L2(cobraC),'k:');
                hold off;
                axis equal

                %% vs time
                figure(tgt+9000)
                
                C = abs(bsxfun(@minus, squeeze(pos(:,:,tgt)), target(:,tgt)));
                D = output.dist(:,:,tgt);
                CD = zeros(size([C;D]));
                CD(1:2:end,:) = C;
                CD(2:2:end,:) = D;
                imagesc(CD);
                  xlabel('Iteration #');
                 ylabel('pIDs');
                set(gca,'YTick',2:2:Nneighbors,'YTickLabel',output.labels);
              caxis([0, 10]);
                keyboard;
                drawnow;

% $$$                 close(pidpair)
% $$$                 close(jj + 8000)
            end
        end
    end

    output.isjammed = isjammed;

end

function output = cobraArms(XY, center, L1, L2)
% returns a 3XN matrix with center, elbow and fiber positions
    if iscolumn(XY)
        XY = XY.';
    end

    ang = XY2TP(XY-center, L1, L2);
    
    output = [ center * ones(size(XY)); ...
               center + L1 * exp(i*ang.tht);...
               XY];
end
    
    

function output = AmIJammed(ER1,ER2,XY,tgt,L1,L2,unitconv)
% determines if a cobra is jammed by it arms 1 and 2 error ratios, X,Y coordinates and target
% position
% 
% ER1 :  (N-1)x1 vectors of error ratios for arm 1
% ER2 :  (N-1)x1 vectors of error ratios for arm 2
% XY  :  Nx1 vector of fiber positions with cobra center (CC) as origin
% tgt :  target position in CC coordinates
% L1  :  length of arm 1
% L2  :  length of arm 2
% unitconv: unit conversion factor into mm
    
    MAXDIST = 0.01; % in mm
    ETHRESH = 0.8; % max good error ratio
    
    tgtang = XY2TP(tgt, L1, L2);
    fibang = XY2TP(XY(1:end-1), L1, L2);
    
    % note: abs(tgt) is the relevant lever for tht
    %       L2 is the lever for phi.
    
    rdtht = unitconv * (fibang.tht - tgtang.tht) * abs(tgt);
    rdphi = unitconv * (fibang.phi - tgtang.phi) * L2;
    dxy   = unitconv * abs(XY(1:end-1) - tgt);

    thtjam = (abs(rdtht) > MAXDIST) & (abs(ER1) > ETHRESH);
    phijam = (abs(rdphi) > MAXDIST) & (abs(ER2) > ETHRESH);

    output = (thtjam | phijam) & (dxy > MAXDIST);
    
end    
