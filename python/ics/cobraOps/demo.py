"""

This file demonstrates how to use the collisions simulation code.
  
"""

import numpy as np

import ics.cobraOps.benchUtils as benchUtils
import ics.cobraOps.cobraUtils as cobraUtils
import ics.cobraOps.plotUtils as plotUtils
import ics.cobraOps.targetUtils as targetUtils
import ics.cobraOps.trajectoryUtils as trajectoryUtils


# Define the target density to use
targetDensity = 1.5

# Get the cobras central positions for the full PFI
cobraCenters = cobraUtils.getCobrasCenters("full")
print("Number of cobras:", len(cobraCenters))

# Define the bench geometry
bench = benchUtils.defineBenchGeometry(cobraCenters, True, True)

# Create a random sample of targets
targetPositions = targetUtils.generateTargets(targetDensity, bench)
print("Number of simulated targets:", len(targetPositions))

# Assign the targets to the cobras and get the fiber positions
(assignedTargets, fiberPositions) = targetUtils.assignTargets(targetPositions, bench)

targetPositions = targetUtils.generateOneTargetPerCobra(bench, 4.)
(assignedTargets, fiberPositions) = targetUtils.assignTargets(targetPositions, bench)
#fiberPositions = targetPositions.copy()
#assignedTargets = np.arange(0,len(cobraCenters))

# Calculate the cobra trajectories
(positiveMovement, nStepsP, nStepsN) = trajectoryUtils.defineThetaMovementDirection(fiberPositions, bench)
(elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)

# Calculate the cobra trajectory collisions
(problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)
print("Number of cobras affected by a collision (  I):", len(np.unique(problematicCobras)))

# Swap the theta movement direction for some of the problematic cobras
cobrasToSwap = trajectoryUtils.getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, True)
positiveMovement = trajectoryUtils.swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
(elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)
(problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)
print("Number of cobras affected by a collision ( II):", len(np.unique(problematicCobras)))

# Swap the theta movement direction for some of the problematic cobras
cobrasToSwap = trajectoryUtils.getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, False)
positiveMovement = trajectoryUtils.swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
(elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)
(problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)
print("Number of cobras affected by a collision (III):", len(np.unique(problematicCobras)))

# Swap the theta movement direction for some of the problematic cobras
cobrasToSwap = trajectoryUtils.getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, True)
positiveMovement = trajectoryUtils.swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
(elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)
(problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)
print("Number of cobras affected by a collision ( IV):", len(np.unique(problematicCobras)))

# Check how many of the collisions are not end point collisions
endCollisions = np.logical_and(trajectoryCollisions[problematicCobras, -1], trajectoryCollisions[nearbyProblematicCobras, -1])
print("Number of cobras unaffected by end collisions: ", len(np.unique(problematicCobras)) - len(np.unique(problematicCobras[endCollisions])))

# Plot the cobra trajectories
plotUtils.createNewFigure("Cobra trajectories", "x position", "y position")

patrolAreaColors = np.full((len(cobraCenters), 4), [0.0, 0.0, 1.0, 0.15])
patrolAreaColors[problematicCobras] = [1.0, 0.0, 0.0, 0.3]
patrolAreaColors[problematicCobras[endCollisions]] = [0.0, 1.0, 0.0, 0.5]
benchUtils.plotBenchGeometry(bench, patrolAreaColors)

footprintColors = np.zeros((elbowTrajectories.shape[0], fiberTrajectories.shape[1], 4))
footprintColors[problematicCobras, :] = [0.0, 0.0, 1.0, 0.05]
trajectoryUtils.plotTrajectories(elbowTrajectories, fiberTrajectories, bench, paintFootprints=True, footprintColors=footprintColors)

cobraColors = np.full((len(cobraCenters), 4), [0.0, 0.0, 1.0, 0.5])
cobraColors[assignedTargets == targetUtils.NULL_TARGET_INDEX] = [1.0, 0.0, 0.0, 0.25]
benchUtils.plotCobras(bench, fiberPositions, cobraColors)

targetColors = np.full((len(targetPositions), 4), [0.4, 0.4, 0.4, 1.0])
targetColors[assignedTargets[assignedTargets != targetUtils.NULL_TARGET_INDEX]] = [1.0, 0.0, 0.0, 1.0]
targetUtils.plotTargets(targetPositions, targetColors)

plotUtils.pauseExecution()

