function animateTraj(pfisim)
% animate trajectories in theta-phi space

cc = pfisim.bench.center;
L1 = pfisim.bench.L1;
L2 = pfisim.bench.L2;
tgt = XY2TP(pfisim.targets - cc, L1, L2);

cut = 0.33;

nsteps = size(pfisim.trajectories,2);
for jj = 1:nsteps
   
    trj = XY2TP(pfisim.trajectories(:,jj) - cc, L1, L2);
    X = mod(trj.tht - tgt.tht - cut,2*pi) + cut - 2*pi;
    Y = trj.phi - tgt.phi;
    plot(X/pi, Y/pi, '.'); 
    axis([-2.1 .2 -.1 1]);
    xlabel('\Delta\theta/\pi')
    ylabel('\Delta\phi/\pi');
    title(sprintf('t = %03d/%d',jj,nsteps));
    drawnow;
    
end

pause(.25)
refline(0,0,0,'r');
refline(0,0,Inf,'r');