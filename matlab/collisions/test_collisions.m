% script to run simFun at some density until a collision occurs.
% generate figures for collisions at that point
close all;
clear;

cat = [];
loops = 0;
dens = 2.0; 
ncollisions = 30;
tic
while length(cat) < ncollisions
    loops = loops + 1
    q = simFun(dens,'full',1,0,'alpha',0.0,'showFigures',0);
    for kk=q.caats(1:2:end)'
        % showMovementAfterSim(q,q.Coll.row(kk));
        figure; cat(end+1) = showCollision(q,kk);
        drawnow;
    end
% $$$     if ~isempty(kk)
% $$$         keyboard;
% $$$     end
end
toc
loops

% alpha = 0
% dens 1.0: 50/72
% dens 1.5: 52/56
% dens 1.8: 50/56
% dens 1.9: 50/49
% dens 2.0: 50/45, 30/43
% dens 2.1: 52/54
% dens 2.2: 50/67
% dens 3.0: 51/103

% alpha = 0.7
% dens 1.0: 239/234
% dens 1.5: 151/130, 100/82
% dens 2.0: 103/88
% dens 2.5: 15/20

