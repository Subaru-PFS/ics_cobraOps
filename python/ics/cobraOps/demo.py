"""

This file demonstrates how to use the collisions simulation code.
  
"""

import ics.cobraOps.cobraUtils as cobraUtils
import ics.cobraOps.targetUtils as targetUtils
import ics.cobraOps.benchUtils as benchUtils
import ics.cobraOps.plotUtils as plotUtils


# Define the target density to use
targetDensity = 1.5

# Get the cobras central positions for the full PFI
centers = cobraUtils.getPFICenters()
print("Number of cobras: " + str(len(centers)))

# Define the bench geometry
bench = benchUtils.defineBenchGeometry(centers, True, True)

# Create a random sample of targets
targetPositions = targetUtils.generateTargets(targetDensity, bench)
print("Number of simulated targets: " + str(len(targetPositions)))

# Assign the targets to the cobras and get the cobra positions
(assignedTargets, fiberPositions) = targetUtils.assignTargets(targetPositions, bench)

# Get the cobras for which the collision could not solved
(problematicCobras, nearbyProblematicCobras) = targetUtils.getProblematicCobras(fiberPositions, bench)
print("Number of unsolved collisions: " + str(len(problematicCobras)/2))

# Plot the cobra-target associations
targetUtils.plotCobraTargetAssociations(fiberPositions, problematicCobras, assignedTargets, targetPositions, bench)
plotUtils.pauseExecution()
