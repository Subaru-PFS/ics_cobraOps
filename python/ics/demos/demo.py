"""

This example demonstrates how to use the collisions simulation code.

"""

import numpy as np
import time as time

from ics.cobraOps import plotUtils
from ics.cobraOps import targetUtils
from ics.cobraOps.Bench import Bench
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector

# Define the target density to use
targetDensity = 1.5

# Load the cobras calibration product
calibrationProduct = CobrasCalibrationProduct("updatedMaps6.xml")

# Create the bench instance
bench = Bench(layout="full", calibrationProduct=calibrationProduct)
print("Number of cobras:", bench.cobras.nCobras)

# Generate the targets
targets = targetUtils.generateRandomTargets(targetDensity, bench)
print("Number of simulated targets:", targets.nTargets)

# Select the targets
selector = DistanceTargetSelector(bench, targets)
selector.run()
selectedTargets = selector.getSelectedTargets()

# Simulate an observation
start = time.time()
simulator = CollisionSimulator(bench, selectedTargets)
simulator.run()
print("Number of cobras involved in collisions:", simulator.nCollisions)
print("Number of cobras unaffected by end collisions: ", simulator.nCollisions - simulator.nEndPointCollisions)
print("Total simulation time (s):", time.time() - start)

# Plot the simulation results
simulator.plotResults(extraTargets=targets, paintFootprints=False)

# Animate one of the trajectory collisions
(problematicCobras,) = np.where(np.logical_and(simulator.collisions, simulator.endPointCollisions == False))

if len(problematicCobras) > 0:
    simulator.animateCobraTrajectory(problematicCobras[0], extraTargets=targets)

# Pause the execution to have time to inspect the figures
plotUtils.pauseExecution()
