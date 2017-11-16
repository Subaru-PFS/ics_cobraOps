"""

This file demonstrates how to use the collisions simulation code.

"""

import numpy as np
import time as time

import ics.cobraOps.plotUtils as plotUtils
import ics.cobraOps.targetUtils as targetUtils

from ics.cobraOps.Bench import Bench
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector


# Define the target density to use
targetDensity = 1.5

# Load the cobras calibration product
start = time.time()
calibrationProduct = CobrasCalibrationProduct("updatedMotorMapsFromThisRun2.xml")

# Create the bench instance
bench = Bench(layout="full", calibrationProduct=calibrationProduct)
print("Number of cobras:", bench.cobras.nCobras)

# Generate the targets
targets = targetUtils.generateTargets(targetDensity, bench)
print("Number of simulated targets:", targets.nTargets)

# Select the targets
selector = DistanceTargetSelector(bench, targets)
selector.run()
selectedTargets = selector.getSelectedTargets()

# Simulate an observation
simulator = CollisionSimulator(bench, selectedTargets)
simulator.run()
print("Number of cobras involved in collisions:", np.sum(simulator.collisionDetected))
print("Number of cobras unaffected by end collisions: ", np.sum(simulator.collisionDetected) - np.sum(simulator.endPointCollisionDetected))
print("Total computation time (s):", time.time() - start)

# Plot the simulation results
start = time.time()
plotUtils.createNewFigure("Collision simulation results", "x position (mm)", "y position (mm)")

# Set the axes limits
limRange = 1.05 * bench.radius * np.array([-1, 1])
xLim = bench.center.real + limRange
yLim = bench.center.imag + limRange
plotUtils.setAxesLimits(xLim, yLim)

# Draw the cobra patrol areas
patrolAreaColors = np.full((bench.cobras.nCobras, 4), [0.0, 0.0, 1.0, 0.15])
patrolAreaColors[simulator.collisionDetected] = [1.0, 0.0, 0.0, 0.3]
patrolAreaColors[simulator.endPointCollisionDetected] = [0.0, 1.0, 0.0, 0.5]
bench.cobras.addPatrolAreasToFigure(colors=patrolAreaColors)

# Draw the cobra links at the final fiber positions
linkColors = np.full((bench.cobras.nCobras, 4), [0.0, 0.0, 1.0, 0.5])
linkColors[simulator.assignedCobras == False] = [1.0, 0.0, 0.0, 0.25]
bench.cobras.addLinksToFigure(simulator.finalFiberPositions, colors=linkColors)

# Draw the cobra trajectories and the trajectory footprints of those that have
# a collision
footprintColors = np.zeros((bench.cobras.nCobras, simulator.trajectories.nSteps, 4))
footprintColors[simulator.collisionDetected, :] = [0.0, 0.0, 1.0, 0.05]
simulator.trajectories.addToFigure(paintFootprints=True, footprintColors=footprintColors)

# Draw all the targets
targets.addToFigure(colors=np.array([0.4, 0.4, 0.4, 1.0]))
selectedTargets.addToFigure(colors=np.array([1.0, 0.0, 0.0, 1.0]))

print("Plotting time (s):", time.time() - start)
plotUtils.pauseExecution()
