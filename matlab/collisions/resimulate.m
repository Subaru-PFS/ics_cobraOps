function output=resimulate(nRepeats, simout, alpha)
% usage output=resimulate(nRepeats, sim_output, alpha)
% reruns the simulation on a set of targets for a given bench.
%
% assumes that the only strategies in play are early-late and
% late-late.  uses both hardstops to minimize fiber travel time.
% 
% run simFun like this:
% simout = simFun(1,'',1,1,'alpha',0)
%
% default alpha for resim is 0.07

if ~exist('alpha','var'), alpha = 0.07; end

if ~exist('simout','var')
    simout = simFun(1.5,'full',1,1,'alpha',alpha);
end

bench   = simout.bench;
targets = simout.targets;

numPos = length(bench.center);
%%%% UNDER CONSTUCTION

ctr = 0;
collTracker = [];
n_collisions = [];

tic
parfor kk=1:nRepeats
    proto = generateTrajectory2(targets, bench);
    Traj  = realizeTrajectory2(proto,bench,simout.Traj.useP,simout.Traj.useL);
    traj  = Traj.traj;
    Coll  = detectCollisionsSparse(traj, bench);

    if ~isempty(find(nonzeros(Coll.detected)))
        indx_collpair = find(sum(Coll.detected,2));
        n_collisions(kk) = length(indx_collpair);
        collTracker = [collTracker indx_collpair.'];
    end
    if ~(mod(kk,100))
        fprintf(1,'loop %d done\n',kk);
    end
end
tt = toc;
fprintf(1,'loop time = %f s\n',tt/nRepeats);

[col_pr ia indx_cp] = unique(collTracker);

occurences = hist(indx_cp, 1:length(col_pr));
ideal_minDist = simout.minDist(col_pr);


output.log = collTracker;
output.nc  = n_collisions;
output.reps = nRepeats;
output.sfout = simout;
output.x = occurences;
output.y = ideal_minDist;

figure(10591)
plot((occurences+randn(size(occurences))/10)/nRepeats, ideal_minDist,'o');
xlabel('Probability of collision');
ylabel('nearest approach on ideal trajectory [mm]');
title(['Nearest approach distributions vs. probability of ' ...
       'collision']);
