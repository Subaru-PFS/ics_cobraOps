% demonstration script using config file "as is"

%% simulate the geometry in the XML configuation file
%% use real maps and arm lengths
%% show collisions
%% output is a structure that includes a set of targets that did
%% not collide in ONE simulation run.
idealSame = simFun(1,'full',0,0);
idealSameRes  = resimulate(1000,idealSame.targets, idealSame.bench);
idealOpp = simFun(1,'full',0,0,'thtDIR',1);
idealOppRes  = resimulate(1000,idealOpp.targets, idealOpp.bench);

realSame = simFun(1,'full',1,1,'alpha',0.07);
realSameRes  = resimulate(1000,realSame.targets, realSame.bench);
realOpp = simFun(1,'full',1,1,'alpha',0.07,'thtDIR',1);
realOppRes  = resimulate(1000,realOpp.targets, realOpp.bench);
% simfun args:
% simFun(triggers, 'type', useRealMaps, useRealLinks)
% option1: 'showMoves', 1== true
% option2: 'alpha', any number (zero for perfect maps)


disp('simFun done')
close all;

figure
cmplx(@plotcircle,EMsim.bench.center, EMsim.bench.L1+EMsim.bench.L2,'k:');
hold on;
plot(EMsim.targets,'o')
hold off;
title('Candidate field selection')
keyboard;


%% re-simulate the target selection a number of times in order to
%% gather statistics on the probability of failure
%% resimulate shows every collision
tic
EMre  = resimulate(1000, EMsim.targets, EMsim.bench);

figure
histogram(EMre);
xlabel('positioner ID');
ylabel('frequency of failure');

figure
pid = 1;%mode(EMre);
showMovementNN(EMsim.trajectories,EMsim.bench,EMsim.coll,pid,EMsim.targets)
plot(EMsim.bench.S1Rm(:,:)'); hold on;
plot(EMsim.bench.S1Rm(pid,:),'k','linewidth',5); hold off;
ylabel('steps to traverse angle?')
