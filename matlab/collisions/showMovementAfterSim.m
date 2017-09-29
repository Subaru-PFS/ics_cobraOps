function showMovementAfterSim(simFunOutput, pid)
% usage: showMovementAfterSim(simFun_output, positioner_ID)
    if ~ishold(gca)
        figure(1000+pid)
    end
    showMovementNN(simFunOutput.Traj.traj, simFunOutput.bench, simFunOutput.Coll, pid, simFunOutput.targets);
end