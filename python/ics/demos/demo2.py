"""

This example demonstrates how to use the collisions simulation code with the
most recent calibration data.

"""

import pathlib
import time
import numpy as np

from ics.cobraCharmer.pfiDesign import PFIDesign
from ics.cobraOps import plotUtils
from ics.cobraOps import targetUtils
from ics.cobraOps.Bench import Bench
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector

# Define the target density to use
targetDensity = 1.5

# Get the calibration product from the XML from Chi-Hung
calibrationProduct = PFIDesign(
    pathlib.Path("/home/jgracia/2020-10-15-theta-slow.xml"))

# Transform the calibration product cobra centers and link lengths units from
# pixels to millimeters
calibrationProduct.centers -= 5048.0 + 3597.0j
calibrationProduct.centers *= np.exp(1j * np.deg2rad(1.0)) / 13.02
calibrationProduct.L1 /= 13.02
calibrationProduct.L2 /= 13.02

# Update the status information for the cobras that are known to be invisible
invisibleCobras = [(15, 1), (15, 55), (17, 37), (29, 57), (31, 14)]

for moduleId, cobraId in invisibleCobras:
    calibrationProduct.setCobraStatus(cobraId, moduleId, invisible=True)

# Update the status information for the cobras that are known to be bad
badCobras = [
    (1, 47), (4, 22), (7, 5), (7, 19), (15, 23), (21, 10), (22, 13), (28, 41),
    (34, 22), (37, 1), (42, 15)]

for moduleId, cobraId in badCobras:
    calibrationProduct.setCobraStatus(
        cobraId, moduleId, brokenTheta=True, brokenPhi=True)

# Update the status information for the cobras that had problems measuring the
# link lengths because they hit fiducial fibers
interferencingCobras = [
    (14, 13), (27, 38), (29, 41), (33, 12), (34, 1), (42, 43)]

for moduleId, cobraId in interferencingCobras:
    calibrationProduct.setCobraStatus(
        cobraId, moduleId, brokenTheta=True, brokenPhi=True)

# Use the median value link lengths in those cobras with zero link lengths
zeroLinkLengths = np.logical_or(
    calibrationProduct.L1 == 0, calibrationProduct.L2 == 0)
calibrationProduct.L1[zeroLinkLengths] = np.median(
    calibrationProduct.L1[~zeroLinkLengths])
calibrationProduct.L2[zeroLinkLengths] = np.median(
    calibrationProduct.L2[~zeroLinkLengths])

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
print("Number of cobras unaffected by end collisions: ",
      simulator.nCollisions - simulator.nEndPointCollisions)
print("Total simulation time (s):", time.time() - start)

# Plot the simulation results
simulator.plotResults(extraTargets=targets, paintFootprints=False)

# Animate one of the trajectory collisions
(problematicCobras,) = np.where(
    np.logical_and(simulator.collisions, simulator.endPointCollisions == False))

if len(problematicCobras) > 0:
    simulator.animateCobraTrajectory(problematicCobras[0], extraTargets=targets)

# Pause the execution to have time to inspect the figures
plotUtils.pauseExecution()
