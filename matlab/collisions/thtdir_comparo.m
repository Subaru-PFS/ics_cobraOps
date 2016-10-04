function thtdir_comparo(nsim,psim)


% 1. generate a set of non-colliding targets assuming perfect
% behavior (alpha = 0)
%
% 2. set alpha t0 0.07, resimulate first move 60X
%
%% n60_07 = resimulate(60,-1); p60_07 = resimulate(60,1);

% The fiber index (1-2394) of the lower index fiber in each collision
% is stored in the output variable.  These are saved in
% resim_results.mat.

%load resim_results
if false
    nsim.log = n60_07;
    psim.log = p60_07;
    ntrials = 60;
end
%% full bench center positions
centers = nsim.sfout.bench.center;
ntrials = nsim.reps;

%% fiber frequency failures -- how often does a particular
%% positioner collide with any neighbor?

indx = 1:length(centers);

nfff = hist(nsim.log,indx);
pfff = hist(psim.log,indx);

%% failure distributions -- how are the frequencies of failure
%% distributed?

%nn = -ntrials/20:ntrials/20:ntrials;
nn=0:ntrials;

nffd = hist(nfff,nn);
pffd = hist(pfff,nn);

%% map of initial colliders
figure(1); clf;
plot(centers,'.'); hold on;
ps=plot(centers(psim.sfout.IR1_colliders),'ro','MarkerFace','r');
po=plot(centers(nsim.sfout.IR1_colliders),'bo','MarkerFace','b');
hold off;
title('Same-sense configuration requires more target replans')
xlabel('focal plane X position [mm]');
ylabel('focal plane Y position [mm]');
legend([po,ps],...
       {sprintf('opposite sense: %.1f replans',length(nsim.sfout.IR1_colliders)/2),...
        sprintf(    'same sense: %.1f replans',length(psim.sfout.IR1_colliders)/2)},...
       'location','NE');

%% collision frequency distribution
figure(2); clf;
bs=bar(fliplr(sort(pfff(pfff>0)))/ntrials,'r','EdgeColor','r'); hold on;
bo=bar(fliplr(sort(nfff(nfff>0)))/ntrials,'b','EdgeColor','b'); 
plot(fliplr(sort(pfff(pfff>0)))/ntrials,'r','LineWidth',2)
hold off;
title('Same sense out has more high-probability colliders')
xlabel('positioners sorted by collision probability');
ylabel('collision probability [P]');
legend([bo,bs],...
       {sprintf('opposite sense: %.1f collisions/field',length(nsim.log)/ntrials),...
        sprintf(    'same sense: %.1f collisions/field',length(psim.log)/ntrials)},...
       'location','SE');

%% number of collisions per trial distribution
figure(3); clf;
histogram(nsim.nc,'FaceColor','b'); hold on;
histogram(psim.nc,'FaceColor','r'); hold off;
xlabel('# collisions');
ylabel('# trials');
legend(sprintf('opposite sense: %s unplanned collisions',format_data(mean(nsim.nc),std(nsim.nc))),...
       sprintf('same sense: %s unplanned collisions',format_data(mean(psim.nc),std(psim.nc))))


%% collision frequency scatter plot
figure(4); clf;
plot(centers,'k.','MarkerSize',1);
hold on;
cmplx(@scatter,centers,max(nfff,1e-6)*0.5,'b','fill');
cmplx(@scatter,centers,max(pfff,1e-6)*0.5,'r','fill'); 
hold off;
title('Collision frequency map on a reachable field');
xlabel('focal plane X position [mm]');
ylabel('focal plane Y position [mm]');
drawnow;
% fix transparency
trans = [180 150];
markers = get(gca,'Children');
for jj=1:length(markers)
    if strcmp(markers(jj).Type,'scatter')
        mh = markers(jj).MarkerHandle;
        color = mh(1).get('FaceColorData');
        color(4) = uint8(trans(jj));
        for kk=1:length(mh)
            mh(kk).set('FaceColorData',color)
        end
    end
end
drawnow;
