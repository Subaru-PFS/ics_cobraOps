"""

This example demonstrates how to use the collisions simulation code with the
most recent calibration data.

"""

import os
import time
import numpy as np

from ics.cobraOps import plotUtils
from ics.cobraOps import targetUtils
from ics.cobraOps.Bench import Bench
from ics.cobraOps.BlackDotsCalibrationProduct import BlackDotsCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector as TargetSelector
from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach

# Initialize the cobra coach instance
os.environ["PFS_INSTDATA_DIR"] = "/home/jgracia/github/pfs_instdata"
cobraCoach = CobraCoach(
    loadModel=True, trajectoryMode=True, rootDir="/home/jgracia/testPFI/")

# Get the calibration product
calibrationProduct = cobraCoach.calibModel

# Print the number of cobras and bad cobras
badCobras = calibrationProduct.status != calibrationProduct.COBRA_OK_MASK
print(f"Number of cobras: {calibrationProduct.nCobras}")
print(f"Number of bad cobras: {np.sum(badCobras)}")

# Fix the phi and tht angles for some of the cobras
wrongAngles = calibrationProduct.phiIn == 0
calibrationProduct.phiIn[wrongAngles] = -np.pi
calibrationProduct.phiOut[wrongAngles] = 0
calibrationProduct.tht0[wrongAngles] = 0
calibrationProduct.tht1[wrongAngles] = (2.1 * np.pi) % (2 * np.pi)
print(f"Number of cobras with wrong phi and tht angles: {np.sum(wrongAngles)}")

# Check if there is any cobra with too short or too long link lengths
tooShortLinks = np.logical_or(
    calibrationProduct.L1 < 1, calibrationProduct.L2 < 1)
tooLongLinks = np.logical_or(
    calibrationProduct.L1 > 5, calibrationProduct.L2 > 5)
print(f"Number of cobras with too short link lenghts: {np.sum(tooShortLinks)}")
print(f"Number of cobras with too long link lenghts: {np.sum(tooLongLinks)}")

# Load the black dots calibration file
calibrationFileName = os.path.join(
    os.environ["PFS_INSTDATA_DIR"],"data/pfi/dot", "black_dots_mm.csv")
blackDotsCalibrationProduct = BlackDotsCalibrationProduct(calibrationFileName)

# Create the bench instance
bench = Bench(cobraCoach, blackDotsCalibrationProduct)

# Generate the targets
targetDensity = 1.5
targets = targetUtils.generateRandomTargets(targetDensity, bench)
print(f"Number of simulated targets: {targets.nTargets}")

# Select the targets
safetyMargin = 0.25
selector = TargetSelector(bench, targets)
selector.run(safetyMargin=safetyMargin)
selectedTargets = selector.getSelectedTargets()

# Simulate an observation
start = time.time()
simulator = CollisionSimulator(bench, selectedTargets)
simulator.run()
nTrajectoryCollisions = simulator.nCollisions - simulator.nEndPointCollisions
print(f"Number of cobras involved in collisions: {simulator.nCollisions}")
print(f"Number of cobras unaffected by end collisions: {nTrajectoryCollisions}")
print(f"Total simulation time (s): {time.time() - start}")

# Plot the simulation results
simulator.plotResults(extraTargets=targets, paintFootprints=False)

# Pause the execution to have time to inspect the figures
plotUtils.pauseExecution()
