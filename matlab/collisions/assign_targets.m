function output=assign_targets(tgt,bench,makefigs)
% usage output=assign_targets(target_list, bench_def, makefigs)
% first arg is target list (complex) or density (real scalar)
% specify bench if you have one, otherwise a full bench will be generated
%
% generate a uniform density field with ntgt target and try to
% assign them to cobras.

if exist('makefigs','var')
    toggle.makefigs = true;
else
    toggle.makefigs = false;
end

% bench definition
if ~exist('bench','var')
    centers = getCentersRails(14);
    centers = bsxfun(@times, [1 exp(i*2*pi/3) exp(i*4*pi/3)], centers);
    centers = reshape(centers,[],1);
    bench = defineBenchGeometry(centers,1,1);
    bench.alpha = 0;
    clear centers
end
ncobras = length(bench.center);

if isreal(tgt) & length(tgt) == 1
    tgtdensity = tgt;
    tgt = generateTargets(tgtdensity, bench);
end

ntgt = length(tgt);

% $$$ %% check min/max radius compliance.
% $$$ LL = sparse(abs(bsxfun(@minus, tgt, bench.center)) < max(bench.rMax)); ...
% $$$ %Who is in patrol area
% $$$ XY = LL*0; % initialize XY with zeros
% $$$ [cc tt] = find(LL); % ccobras and targets that are nonzero
% $$$ for jj = 1:length(cc) % for every cobra check assigment 
% $$$     xy   = tgt(tt(jj)) - bench.center(cc(jj)); % get the local coordinate
% $$$     dst  = abs(xy);
% $$$     if (dst > bench.rMin(cc(jj)) & dst < bench.rMax(cc(jj)))
% $$$         XY(cc(jj),tt(jj)) = xy;
% $$$     end
% $$$ end
% $$$ clear xy dst cc tt

% full calc of distances
XY   = tgt - bench.center;
dst  = abs(XY);
% filter for reachability
IN_PATROL_REGION = sparse( dst > bench.rMin & dst < bench.rMax );

% sparsify XY (local coordinate position) and dst (local radial position)
XY  = IN_PATROL_REGION .* XY;
dst = IN_PATROL_REGION .* dst;

nvalid = length(find(sum(IN_PATROL_REGION,1)));
% $$$ fprintf(1, 'Found %d valid targets (%5.1f %%)\n', nvalid, 100*nvalid/ntgt);

[srtDIST srtINDX] = sort(dst,2);
n_cobra_4_tgt = sum(IN_PATROL_REGION,1); % # of cobras for each target
n_tgt_4_cobra = sum(IN_PATROL_REGION,2);
n_tgt_4_cobra_initial = n_tgt_4_cobra; % bookkeeping variable
max_tgts_per_cobra = max(n_tgt_4_cobra);

%% reformat srtINDX so that we only have the relevant target indices for each cobra in distance
%% order.  To start, each cobra is assigned the nearest target. the targets are NOT uniquely
%% assigned.
for jj=1:ncobras
    if n_tgt_4_cobra(jj) > 0
        %srtDIST = 0 for all targets out of range, so the relevant targets are at the end of
        %every row of the matrix.  The next line here cycles them to the front of the list, so
        %the nearest neighbors are all in in column 1
        srtDIST(jj,:) = circshift(srtDIST(jj,:), [0 n_tgt_4_cobra(jj)]);
        srtINDX(jj,:) = circshift(srtINDX(jj,:), [0 n_tgt_4_cobra(jj)]);
    end
end
srtINDX = srtINDX .* sign(srtDIST); % zero out the indices of the unreachable targets
srtTGT_ID = srtINDX(:,1:(max_tgts_per_cobra + 1));
clear srtINDX

%%% reassign targets to cobras until assignments are unique.  The goal here is to make the first
%%% column of srtTGT_ID the target ID for each cobra.
Multiply_Assigned_Targets = 1;
ctr.uniquification_iterations = 0;
ctr.number_singly_assigned = [0 0 0 0];
while Multiply_Assigned_Targets > 0
    [uu it iu] = unique(srtTGT_ID(:,1)); % index of iu is a cobra, uu(iu) is a target ID.
                                         % srtTGT_ID = uu(iu), uu = strTGT(it,1)
    Multiply_Assigned_Targets = 0;
    ctr.uniquification_iterations = ctr.uniquification_iterations + 1;
% $$$     disp('---------------------------');
    for jj = 1:length(uu) % loop over assigned targets
        if uu(jj) == 0, continue; end; % 0 is the tgtID of the null target. skip it.
        cobras_4_this_tgt = find(iu == jj);
        if length(cobras_4_this_tgt) > 1 % this tgt [uu(jj)] is assigned to multiple cobras
                                         
            % of the cobras assigned to this target, find the number of targets available to
            % each cobra ...
            n_tgts_per_cobra = n_tgt_4_cobra(cobras_4_this_tgt);
            % ... and the cobras for which this is the only target
            single_tgt_cobras = cobras_4_this_tgt(find(n_tgts_per_cobra == 1));
            switch length(single_tgt_cobras)
                % decide action based on the # cobras for which there is only one choice
                % note that each target can be assigned to no more than 3 cobras (by geometry)
              case 0 % all cobras have multiple choices
                % HACK SOLUTION: give it to the first one in the list
                use_this_cobra = cobras_4_this_tgt(1);
                dont_use_these = cobras_4_this_tgt(2:end);
                ctr.number_singly_assigned(1) = ctr.number_singly_assigned(1) + 1;
              case 1 
                % give it to the singleton
                use_this_cobra = single_tgt_cobras;
                dont_use_these = cobras_4_this_tgt(find(n_tgts_per_cobra > 1));
                ctr.number_singly_assigned(2) = ctr.number_singly_assigned(2) + 1;
              case {2,3} % multiple cobras have only one choice
                % give it to the (radially) closest cobra
                [dst tmpindx] = min(srtDIST(single_tgt_cobras,1));
                use_this_cobra = single_tgt_cobras(tmpindx);
                dont_use_these = cobras_4_this_tgt(cobras_4_this_tgt ~= use_this_cobra);
                ctr.number_singly_assigned(3) = ctr.number_singly_assigned(3) + 1;
            end
% $$$             fprintf(1,'tgt%04d: use %04d(%d); don''t use %04d(%d)', uu(jj),use_this_cobra, ...
% $$$                     full(n_tgt_4_cobra(use_this_cobra)), dont_use_these(1), ...
% $$$                     full(n_tgt_4_cobra(dont_use_these(1))));
            for kk = 1:length(dont_use_these)
                not_me = dont_use_these(kk);
                % for cobras that will not be assigned to this target, shift the assignment to
                % the next nearest target.
                srtTGT_ID(not_me,:) = circshift(srtTGT_ID(not_me,:), [0 -1]);
                srtDIST(not_me,:)   = circshift(  srtDIST(not_me,:), [0 -1]);
                % also decrement the number of targets reachable by this cobra
                n_tgt_4_cobra(not_me) = n_tgt_4_cobra(not_me) - 1;
                % increment "Multiply_Assigned_Targets" so the while loop executes again.
                Multiply_Assigned_Targets = Multiply_Assigned_Targets + 1;
% $$$                 if kk > 1
% $$$                     fprintf(1,' %04d(%d)',dont_use_these(kk), full(n_tgt_4_cobra(dont_use_these(kk))));
% $$$                 end
            end
% $$$             disp(' ');
        end
        clear cobras_4_this_tgt dst tmpindx dont_use_these use_this_cobra;
    end
% $$$     NEVER HAPPENS.
% $$$     %% sometimes, 2 cobras will reassign to the same target, so
% $$$     %% check for uniqueness again (multiplicity in iu is bad)
% $$$     if Multiply_Assigned_Targets == 0
% $$$         [uu it iu] = unique(srtTGT_ID(:,1));
% $$$         hist_iu = hist(iu(iu>1),0:ncobras);
% $$$         if max(hist_iu) > 1
% $$$             disp('***CHECK FOR COLLISION ERROR***')
% $$$             %            Multiply_Assigned_Targets = 1;
% $$$         end
% $$$     end
end 
%% at this point, the targets are uniquely assigned to cobras.
clear uu it iu jj kk iuvalue n_tgts_per_cobra single_tgt_cobras Multiply_Assigned_Targets...
    not_me dont_use_these use_this_cobra

% $$$ fprintf(1,'targets assigned in %d iterations\n',ctr.uniquification_iterations);

collisions = ncobras;

tgt_home = bench.home0; % initially, home positions in SS move out. tgt_home can be changed by collisions
                        % on unassigned cobras.
ctr.reassigns=0;
while (collisions > 0)
    %% Find end-point collisions and reassign targets as needed.
    tgtID = srtTGT_ID(:,1); % ncobras X 1 list of target IDs
    isassigned = (tgtID>0); 
    
    % assigned target is either the moveable home or a real target
    tgt_assigned = tgt_home;
    tgt_assigned(isassigned) = tgt(tgtID(isassigned));
    
    % find targets that violate  cobra physical extents at target location
    coll_dist = CalcDistanceMatrix(tgt_assigned, bench);
    coll_dist = min(coll_dist,coll_dist');
    % "Rule 2" rows and columns
    [R2.row R2.col] = find(coll_dist > 0 & coll_dist < bench.minDist);
    collisions = length(R2.row);

    %% note, if  collisions == 0, then nothing else happens and this loop exits.
    
    %% first check and report on triple collisions
    occurences = hist(R2.row,1:ncobras);
    if length(find(occurences > 1)) > 0
        disp('Warning -- three-way collision(s) detected');
        disp(find(occurences > 1));
        %% currently, the for loop below handles only binary collisions.  a proper target assignment
        %% algorithm would handle ternary collisions cleanly.
    end
    %    for jj = 1:length(R2.row)
    %% solve the binary collisions
    for jj = find(R2.row < R2.col)'
% $$$         if R2.row(jj) < R2.col(jj) % only need to check half of each collision
            cobras = [R2.row(jj) R2.col(jj)];
            
            %% if both cobras are assigned, then weed out alternate targets that are already
            %% assigned, then move the cobra with more flexibility.
            for kk = cobras
                % tmpindx is a reversed list of alternate targets for cobra kk that are already
                % assigned to other cobras
                tmpindx = fliplr(find([0 ismember(srtTGT_ID(kk,2:n_tgt_4_cobra(kk)), srtTGT_ID(:,1))]));
                for tt = tmpindx
                    % for each entry in tmpindx, cycle that entry to the back of the list.
                    srtTGT_ID(kk,tt:end) = circshift(srtTGT_ID(kk,tt:end), [0 -1]);
                    srtDIST(kk,tt:end)   = circshift(  srtDIST(kk,tt:end), [0 -1]);
                end
                n_tgt_4_cobra(kk) = n_tgt_4_cobra(kk) - length(tmpindx);
                clear tmpindx
            end
            %% at this point, the alternate targets for each cobra are unassigned.

            ntgts = full(n_tgt_4_cobra(cobras));
% $$$             fprintf(1,'%4d: %d; %4d: %d.  ',cobras(1), ntgts(1), cobras(2), ntgts(2));
            
            [nmin jmin] = min(ntgts);
            %% resolve spatial conflicts
            if nmin > 0 % both cobras are assigned
                if ntgts(1) > ntgts(2)
                    not_me = cobras(1);
                elseif ntgts(2) > ntgts(1)
                    not_me = cobras(2);
                else % same # of degrees of freedom
                    if ntgts(1) == 1 % only one choice each - jettison the further-out one
                        [tmp tmpindx] = max(srtDIST(cobras,1));
                        not_me = cobras(tmpindx);
                    else % otherwise, choose randomly
                        not_me = cobras(ceil(rand(1) + 0.5));
                    end
                end
                srtTGT_ID(not_me,:) = circshift(srtTGT_ID(not_me,:), [0 -1]);
                srtDIST(not_me,:)   = circshift(  srtDIST(not_me,:), [0 -1]);
                n_tgt_4_cobra(not_me) = n_tgt_4_cobra(not_me) - 1;
            else % one cobra is unassigned.  unassigned cobras can't collide, so jmin specifies the
                 % unassigned cobra.
                not_me = cobras(jmin);
                my_angles = XY2TP(tgt_home(not_me) - bench.center(not_me), ...
                                  bench.L1(not_me), bench.L2(not_me));
                my_angles.tht = my_angles.tht + pi/3;
                tgt_home(not_me) = (bench.center(not_me) + ...
                                     bench.L1(not_me) * exp(i*my_angles.tht) + ...
                                     bench.L2(not_me) * exp(i* ...
                                                            (my_angles.tht + my_angles.phi)) );
            end
% $$$             fprintf(1, 'moved %4d\n',not_me);
            ctr.reassigns = ctr.reassigns + 1;
            clear jmin nmin ntgts tmp tmpindx tt
% $$$         end
    end
end

clear R2 cobras curr_collisions last_collisions not_me occurences

TP = XY2TP(tgt_assigned - bench.center, bench.L1, bench.L2);

n_assigned = sum(srtTGT_ID(:,1) > 0);
% $$$ fprintf(1,'Assigned %d targets (%5.1f %% of valid)\n', n_assigned, n_assigned/nvalid*100);
% $$$ fprintf(1,'# unassigned cobras: %d\n', sum(srtTGT_ID(:,1) <= 0));
% $$$ fprintf(1,'# reassigned cobras: %d\n', ctr.reassigns);
fprintf(1,'mean target offset from cobra axis: %.2f mm\n',full(mean(srtDIST(isassigned,1))));

cobra_arms = [bench.center ...
              bench.center + bench.L1.*exp(i*TP.tht) ...
              tgt_assigned].'; % cobra [center; elbow; fiber] array (3 x ncobras)

if toggle.makefigs
    figure(1001)
    plot(tgt,'bx');
    hold on;
    plot(cobra_arms, 'ko-');
    plot(tgt_assigned(tgtID > 0),'bo','markerface','g');
    plot(tgt_assigned(tgtID <= 0),'bo','markerface','k');
    hold off;
    cmplx(@plotcircle, bench.center, bench.rMax, 'k:');
% $$$     cmplx(@plotcircle, cobra_arms(2:3,R2.col), bench.minDist/2, 'r');
% $$$     cmplx(@text, bench.center(R2.col)+1, cellstr(num2str(R2.col)));
    axis equal;

    figure(1002)
    subplot(211);
    hh1 = histogram(n_tgt_4_cobra_initial);
    x_h = mean([hh1.BinEdges(1:end-1) ; hh1.BinEdges(2:end)]);
    est_density = sum(hh1.Values.*x_h) / sum(hh1.Values);
    if ~exist('tgtdensity','var')
        tgtdensity = est_density;
    end
    hold on;
    plot(x_h, est_density.^x_h * exp(-est_density) ./ factorial(x_h) * ncobras, 'r');
    hold off;
    ylabel('# cobras');
    xlabel('# targets');
    title(sprintf(['Distribution of targets per cobra with target ' ...
                   'density of %.1f per patrol area'],tgtdensity));
    legend('"observed" population distribution','Poisson distribution');
    
    subplot(212);
    hh2 = histogram(srtDIST(tgtID>0,1));
    hold on;
    xmax = mean(bench.rMax);
    slope = 2*n_assigned*hh2.BinWidth/xmax^2; 
    plot([0 xmax xmax], [0 xmax 0]*slope,'r');
    hold off;
    xlabel('target radial position in patrol area [mm]');
    ylabel('# targets');
    title(['target spatial distribution for n_t_g_t = ' num2str(ntgt)]);
    legend('"observed" population distribution','Uniform distribution','location','N');
end

output.tgt = tgt_assigned;
output.rem = tgt(~ismember(1:length(tgt), tgtID));
output.all = tgt;
output.isassigned = srtTGT_ID(:,1) > 0;
output.at_home    = abs(tgt_assigned - bench.home0) < 1e-9;
output.srtDIST = srtDIST;
output.srtTGT_ID = srtTGT_ID;

% $$$ histogram(n_cobra_4_tgt);
% $$$ xlabel('# cobras');
% $$$ ylabel('# targets');

% $$$ figure(213);
% $$$ histogram(bench.rMax);
% $$$ junk = allstats(bench.rMax);
% $$$ refline(junk.median,0,Inf,'r');
% $$$ refline(junk.mean,0,Inf,'r:');
% $$$ xlabel('rMax [mm]');
% $$$ title('patrol radius distribution');
% $$$ legend('distribution','median','mean');