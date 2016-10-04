function miniresim(input,pid)
% inputs: simFun output, pid
% resimulates and plots trajectories around pid

proto = generateTrajectory2(input.targets, input.bench);
Traj  = realizeTrajectory2(proto, input.bench, input.Traj.useP, input.Traj.useL);
Coll  = detectCollisionsSparse(Traj.traj, input.bench);

showMovementNN(Traj.traj, input.bench, Coll, pid, input.targets);

size(Traj.traj)