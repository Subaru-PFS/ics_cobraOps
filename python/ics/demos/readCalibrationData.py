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
bench = Bench(layout="calibration", calibrationProduct=calibrationProduct)
print("Number of cobras:", bench.cobras.nCobras)

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
