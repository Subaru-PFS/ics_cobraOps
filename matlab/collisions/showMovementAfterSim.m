function showMovementAfterSim(simFunOutput, pid)
    if ~ishold(gcf)
        figure(1e6+pid)
    end
    showMovementNN(simFunOutput.Traj.traj, simFunOutput.bench, simFunOutput.Coll, pid, simFunOutput.targets);
end