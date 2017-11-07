"""

Calculates some cobra collision statistics.

"""

import os
import numpy as np

import ics.cobraOps.benchUtils as benchUtils
import ics.cobraOps.cobraUtils as cobraUtils
import ics.cobraOps.plotUtils as plotUtils
import ics.cobraOps.targetUtils as targetUtils
import ics.cobraOps.trajectoryUtils as trajectoryUtils


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

# Get the cobras central positions for the full PFI
cobraCenters = cobraUtils.getCobrasCenters("full")

# Calculate the collisions for each maxDist-density combination
for i in range(len(maxDistanceArray)):
    print("TargetDensity", targetDensityArray[i], "MaxDistance", maxDistanceArray[i])
    
    # Define the bench geometry
    bench = benchUtils.defineBenchGeometry(cobraCenters, True, True)

    # Create a random sample of targets
    targetPositions = targetUtils.generateTargets(targetDensityArray[i], bench)

    # Get the indices and distances of the targets that can be reached by each
    # cobra, making sure that their distance is smaller than maxDistance
    (targetIndices, targetDistances) = targetUtils.getAccesibleTargets(targetPositions, bench, maxDistanceArray[i])

    # Assign targets to cobras and set the fiber positions
    if solveCollisions:
        # Assign a single target to each cobra based on their separation
        assignedTargets = targetUtils.assignTargetsByDistance(targetIndices, targetDistances)
        
        # Solve possible cobra collision reassigning some of the cobras and
        # set the cobra fiber position
        fiberPositions = targetUtils.solveCobraCollisions(assignedTargets, targetIndices, targetPositions, bench)
    else:
        # Assign a single target to each cobra randomly       
        assignedTargets = targetUtils.assignTargetsRandomly(targetIndices, targetDistances)
        
        # Set the cobra fiber positions
        cobrasWithTarget = assignedTargets != targetUtils.NULL_TARGET_INDEX
        fiberPositions = bench["home0"].copy()
        fiberPositions[cobrasWithTarget] = targetPositions[assignedTargets[cobrasWithTarget]]  

    # Calculate the cobra trajectories
    (positiveMovement, nStepsP, nStepsN) = trajectoryUtils.defineThetaMovementDirection(fiberPositions, bench)
    (elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)

    # Calculate the cobra trajectory collisions
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)

    # Swap the theta movement direction for some of the problematic cobras
    cobrasToSwap = trajectoryUtils.getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, True)
    positiveMovement = trajectoryUtils.swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
    (elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)

    # Swap the theta movement direction for some of the problematic cobras
    cobrasToSwap = trajectoryUtils.getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, False)
    positiveMovement = trajectoryUtils.swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
    (elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)

    # Swap the theta movement direction for some of the problematic cobras
    cobrasToSwap = trajectoryUtils.getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, True)
    positiveMovement = trajectoryUtils.swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
    (elbowTrajectories, fiberTrajectories) = trajectoryUtils.calculateTrajectories(fiberPositions, positiveMovement, bench)
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = trajectoryUtils.detectCollisions(elbowTrajectories, fiberTrajectories, bench)
    
    # Fill the statistics arrays
    unassignedCobrasArray[i] = np.sum(assignedTargets == targetUtils.NULL_TARGET_INDEX)
    collisionsArray[i] = len(np.unique(problematicCobras))
    endCollisions = np.logical_and(trajectoryCollisions[problematicCobras, -1], trajectoryCollisions[nearbyProblematicCobras, -1])
    endCollisionsArray[i] = len(np.unique(problematicCobras[endCollisions]))


# Plot the collision probabilities
colorMap = plotUtils.getColorMap()
colors = np.zeros((len(maxDistanceArray), 3))

for i, targetDensity in enumerate(targetDensitySteps):
    colors[targetDensityArray == targetDensity] = colorMap[i]

plotUtils.createNewFigure("Cobra unassignment probability", "Target distance", "Probability", size=(6, 6), aspectRatio="auto")
plotUtils.setAxesLimits([2.4, 5.1], [-0.01, 0.61])
plotUtils.plt.scatter(maxDistanceArray, unassignedCobrasArray / len(cobraCenters), facecolor=colors, s=8)
plotUtils.saveFigure(os.path.join(outputDir, "unassignmentProbability" + suffix + ".png"))

plotUtils.createNewFigure("Cobra collision probability", "Target distance", "Probability", size=(6, 6), aspectRatio="auto")
plotUtils.setAxesLimits([2.4, 5.1], [-0.01, 0.16])
plotUtils.plt.scatter(maxDistanceArray, collisionsArray / len(cobraCenters), facecolor=colors, s=8)
plotUtils.saveFigure(os.path.join(outputDir, "collisionProbability" + suffix + ".png"))

plotUtils.createNewFigure("Cobra trajectory collision probability", "Target distance", "Probability", size=(6, 6), aspectRatio="auto")
plotUtils.setAxesLimits([2.4, 5.1], [-0.01, 0.16])
plotUtils.plt.scatter(maxDistanceArray, (collisionsArray - endCollisionsArray) / len(cobraCenters), facecolor=colors, s=8)
plotUtils.saveFigure(os.path.join(outputDir, "trajCollisionProbability" + suffix + ".png"))

plotUtils.pauseExecution()

