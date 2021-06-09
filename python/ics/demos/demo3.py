"""

This example demonstrates how to use the collisions simulation code with the
most recent calibration data.

"""

import os
import time
import logging
import numpy as np

from ics.cobraOps import plotUtils
from ics.cobraOps import targetUtils
from ics.cobraOps.Bench import Bench
from ics.cobraOps.CollisionSimulator2 import CollisionSimulator2
from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector
from procedures.moduleTest.cobraCoach import CobraCoach

# Disable the matplotlit warnings
logging.getLogger("matplotlib.font_manager").disabled = True

# Initialize the cobra coach instance
os.environ["PFS_INSTDATA_DIR"] = "/home/jgracia/github/pfs_instdata"
cobraCoach = CobraCoach(
    "fpga", loadModel=False, trajectoryMode=True,
    rootDir="/home/jgracia/testPFI/")
cobraCoach.loadModel(version="ALL", moduleVersion="final_20210512")

# Get the calibration product
calibrationProduct = cobraCoach.calibModel

# Set some dummy center positions and phi angles for those cobras that have
# zero centers
zeroCenters = calibrationProduct.centers == 0
calibrationProduct.centers[zeroCenters] = np.arange(np.sum(zeroCenters)) * 300j
calibrationProduct.phiIn[zeroCenters] = -np.pi
calibrationProduct.phiOut[zeroCenters] = 0
print("Cobras with zero centers: %i" % np.sum(zeroCenters))

# Transform the calibration product cobra centers and link lengths units from
# pixels to millimeters
calibrationProduct.centers -= 5048.0 + 3597.0j
calibrationProduct.centers *= np.exp(1j * np.deg2rad(1.0)) / 13.02
calibrationProduct.L1 /= 13.02
calibrationProduct.L2 /= 13.02

# Use the median value link lengths in those cobras with zero link lengths
zeroLinkLengths = np.logical_or(
    calibrationProduct.L1 == 0, calibrationProduct.L2 == 0)
calibrationProduct.L1[zeroLinkLengths] = np.median(
    calibrationProduct.L1[~zeroLinkLengths])
calibrationProduct.L2[zeroLinkLengths] = np.median(
    calibrationProduct.L2[~zeroLinkLengths])
print("Cobras with zero link lenghts: %i" % np.sum(zeroLinkLengths))

# Use the median value link lengths in those cobras with too long link lengths
tooLongLinkLengths = np.logical_or(
    calibrationProduct.L1 > 100, calibrationProduct.L2 > 100)
calibrationProduct.L1[tooLongLinkLengths] = np.median(
    calibrationProduct.L1[~tooLongLinkLengths])
calibrationProduct.L2[tooLongLinkLengths] = np.median(
    calibrationProduct.L2[~tooLongLinkLengths])
print("Cobras with too long link lenghts: %i" % np.sum(tooLongLinkLengths))

# Create the bench instance
bench = Bench(layout="calibration", calibrationProduct=calibrationProduct)
print("Number of cobras:", bench.cobras.nCobras)

# Generate the targets
targetDensity = 1.5
targets = targetUtils.generateRandomTargets(targetDensity, bench)
print("Number of simulated targets:", targets.nTargets)

# Select the targets
selector = DistanceTargetSelector(bench, targets)
selector.run()
selectedTargets = selector.getSelectedTargets()

# Simulate an observation
start = time.time()
simulator = CollisionSimulator2(bench, cobraCoach, selectedTargets)
simulator.run()
print("Number of cobras involved in collisions:", simulator.nCollisions)
print("Number of cobras unaffected by end collisions: ",
      simulator.nCollisions - simulator.nEndPointCollisions)
print("Total simulation time (s):", time.time() - start)

# Plot the simulation results
simulator.plotResults(extraTargets=targets, paintFootprints=False)

# Pause the execution to have time to inspect the figures
plotUtils.pauseExecution()
