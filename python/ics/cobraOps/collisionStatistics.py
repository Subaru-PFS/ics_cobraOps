"""

Calculates some cobra collision statistics.

"""

import os
import numpy as np

import ics.cobraOps.plotUtils as plotUtils
import ics.cobraOps.targetUtils as targetUtils

from ics.cobraOps.Bench import Bench
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector


# Decide if the code should try to solve cobra collisions assigning new targets
# to problematic cobras
solveCollisions = False

# Set the directory where the figures should be saved
outputDir = "/home/jgracia/"

# Set the figures suffix to use
suffix = "-solved" if solveCollisions else ""

# Define the maximum target-to-cobra distance steps
maxDistanceSteps = np.arange(2.5, 5.01, 0.05)

# Define the target density steps
targetDensitySteps = np.arange(2, 5.01, 1)

# Set the total number of repetitions for a given combination
repetitions = 5

# Create all the necessary arrays
maxDistanceArray = np.repeat(np.tile(maxDistanceSteps, len(targetDensitySteps)), repetitions)
targetDensityArray = np.repeat(np.repeat(targetDensitySteps, len(maxDistanceSteps)), repetitions)
unassignedCobrasArray = np.zeros(len(maxDistanceArray))
collisionsArray = np.zeros(len(maxDistanceArray))
endCollisionsArray = np.zeros(len(maxDistanceArray))

# Load the cobras calibration product
calibrationProduct = CobrasCalibrationProduct("updatedMotorMapsFromThisRun2.xml")

# Calculate the collisions for each maxDist-density combination
for i in range(len(maxDistanceArray)):
    print("TargetDensity", targetDensityArray[i], "MaxDistance", maxDistanceArray[i])
    
    # Create the bench instance
    bench = Bench(layout="full", calibrationProduct=calibrationProduct)
    
    # Create a random sample of targets
    targets = targetUtils.generateRandomTargets(targetDensityArray[i], bench)
    
    # Select the targets
    selector = DistanceTargetSelector(bench, targets)
    selector.run(maximumDistance=maxDistanceArray[i], solveCollisions=solveCollisions)
    selectedTargets = selector.getSelectedTargets()
    
    # Simulate an observation
    simulator = CollisionSimulator(bench, selectedTargets)
    simulator.run()
    
    # Fill the statistics arrays
    unassignedCobrasArray[i] = np.sum(simulator.assignedCobras == False)
    collisionsArray[i] = simulator.nCollisions
    endCollisionsArray[i] = simulator.nEndPointCollisions


# Plot the collision probabilities
colorMap = plotUtils.getColorMap()
colors = np.zeros((len(maxDistanceArray), 3))

for i, targetDensity in enumerate(targetDensitySteps):
    colors[targetDensityArray == targetDensity] = colorMap[i]

plotUtils.createNewFigure("Cobra unassignment probability", "Target distance", "Probability", size=(6, 6), aspectRatio="auto")
plotUtils.setAxesLimits([2.4, 5.1], [-0.01, 0.61])
plotUtils.plt.scatter(maxDistanceArray, unassignedCobrasArray / bench.cobras.nCobras, facecolor=colors, s=8)
plotUtils.saveFigure(os.path.join(outputDir, "unassignmentProbability" + suffix + ".png"))

plotUtils.createNewFigure("Cobra collision probability", "Target distance", "Probability", size=(6, 6), aspectRatio="auto")
plotUtils.setAxesLimits([2.4, 5.1], [-0.01, 0.16])
plotUtils.plt.scatter(maxDistanceArray, collisionsArray / bench.cobras.nCobras, facecolor=colors, s=8)
plotUtils.saveFigure(os.path.join(outputDir, "collisionProbability" + suffix + ".png"))

plotUtils.createNewFigure("Cobra trajectory collision probability", "Target distance", "Probability", size=(6, 6), aspectRatio="auto")
plotUtils.setAxesLimits([2.4, 5.1], [-0.01, 0.16])
plotUtils.plt.scatter(maxDistanceArray, (collisionsArray - endCollisionsArray) / bench.cobras.nCobras, facecolor=colors, s=8)
plotUtils.saveFigure(os.path.join(outputDir, "trajCollisionProbability" + suffix + ".png"))

plotUtils.pauseExecution()

