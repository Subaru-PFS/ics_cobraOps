"""

This example demonstrates how to load the most recent calibration data files.

"""

import pathlib
import numpy as np

from ics.cobraCharmer.pfiDesign import PFIDesign
from ics.cobraOps import plotUtils
from ics.cobraOps.Bench import Bench

"""
# Define the location of the XML calibration file for each cobra module
os.environ["PFS_INSTDATA_DIR"] = "/home/jgracia/github/pfs_instdata"
modules_directory = pathlib.Path(
    os.environ["PFS_INSTDATA_DIR"]) / "data" / "pfi" / "modules"
modules_names = [("SC%2i" % i).replace(" ", "0") for i in range(1, 43)]
modules_file_names = [
    modules_directory / module_name / (
        module_name + "_final.xml") for module_name in modules_names]

# Get the calibration product
calibrationProduct = PFIDesign(None)
calibrationProduct.loadModelFiles(modules_file_names)
"""

# Get the calibration product from the pfs_instdata repository
calibrationProduct = PFIDesign(pathlib.Path(
    "/home/jgracia/github/pfs_instdata/data/pfi/modules/ALL/ALL_final.xml"))

# Set some dummy center positions and phi angles for those cobras that have
# zero centers
zeroCenters = calibrationProduct.centers == 0
calibrationProduct.centers[zeroCenters] = np.arange(np.sum(zeroCenters)) * 300j
calibrationProduct.phiIn[zeroCenters] = -np.pi
calibrationProduct.phiOut[zeroCenters] = 0

# Use the median value link lengths in those cobras with zero link lengths
zeroLinkLengths = np.logical_or(
    calibrationProduct.L1 == 0, calibrationProduct.L2 == 0)
calibrationProduct.L1[zeroLinkLengths] = np.median(
    calibrationProduct.L1[~zeroLinkLengths])
calibrationProduct.L2[zeroLinkLengths] = np.median(
    calibrationProduct.L2[~zeroLinkLengths])

# Transform the calibration product cobra centers and link lengths units from
# pixels to millimeters
calibrationProduct.centers -= 5048.0 + 3597.0j
calibrationProduct.centers *= np.exp(1j * np.deg2rad(1.0)) / 13.02
calibrationProduct.L1 /= 13.02
calibrationProduct.L2 /= 13.02

# Create the bench instance
bench = Bench(layout="calibration", calibrationProduct=calibrationProduct)
print("Number of cobras:", bench.cobras.nCobras)

print(bench.getCobrasNeighbors(np.where(bench.cobras.hasProblem)[0]))
# Plot the bench
plotUtils.createNewFigure(
    "Calibration data", "x position (mm)", "y position (mm)")
limRange = 1.05 * bench.radius * np.array([-1, 1])
xLim = bench.center.real + limRange
yLim = bench.center.imag + limRange
plotUtils.setAxesLimits(xLim, yLim)
patrolAreaColors = np.full((bench.cobras.nCobras, 4), [0.0, 0.0, 1.0, 0.15])
patrolAreaColors[
    (bench.cobras.status & PFIDesign.COBRA_INVISIBLE_MASK) != 0] = [
        1.0, 0.0, 0.0, 0.3]
patrolAreaColors[
    (bench.cobras.status & PFIDesign.COBRA_BROKEN_MOTOR_MASK) != 0] = [
        0.0, 1.0, 0.0, 0.5]
bench.cobras.addPatrolAreasToFigure(
    colors=patrolAreaColors, paintHardStops=True, paintBlackDots=False)
linkColors = np.full(patrolAreaColors.shape, [0.0, 0.0, 1.0, 0.5])
bench.cobras.addLinksToFigure(bench.cobras.home0, colors=linkColors)

# Pause the execution to have time to inspect the figures
plotUtils.pauseExecution()
