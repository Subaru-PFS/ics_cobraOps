"""

This file demonstrates how to load the cobras calibrated data XML files.

"""

import os
import pathlib
import numpy as np

from ics.cobraCharmer.pfiDesign import PFIDesign
from ics.cobraOps import plotUtils
from ics.cobraOps.Bench import Bench

# Define the location of the XML calibration file for each cobra module
os.environ["PFS_INSTDATA_DIR"] = "/home/jgracia/github/pfs_instdata"
modules_directory = pathlib.Path(os.environ["PFS_INSTDATA_DIR"]) / "data" / "pfi" / "modules"
modules_names = [("SC%2i" % i).replace(" ", "0") for i in range(1, 43)]
modules_file_names = [modules_directory / module_name / (module_name + "_final.xml") for module_name in modules_names]

# Load the cobras calibration product from the pfs_instdata project
# calibrationProduct = PFIDesign.loadPfi(moduleVersion="final")
calibrationProduct = PFIDesign(None)
calibrationProduct.loadModelFiles(modules_file_names)

# Save the link length information
L1 = calibrationProduct.L1
L2 = calibrationProduct.L2

# Load the latest XML calibration product form Chi-Hung
calibrationProduct = PFIDesign(pathlib.Path("/home/jgracia/2020-10-15-theta-slow.xml"))

# Replace the link lengths in those cobras with zero values
wrong_link_lengths = np.logical_or(calibrationProduct.L1 == 0, calibrationProduct.L2 == 0)
calibrationProduct.L1[wrong_link_lengths] = L1[wrong_link_lengths]
calibrationProduct.L2[wrong_link_lengths] = L2[wrong_link_lengths]

# Transform the calibration product cobra centers and link lengths units from
# pixels to millimeters
calibrationProduct.centers -= 5048.0 + 3597.0j
calibrationProduct.centers *= np.exp(1j * np.deg2rad(1.0)) / 13.02
calibrationProduct.L1 /= 13.02
calibrationProduct.L2 /= 13.02

# Update the status information for the cobras that are known to be invisible
invisibleCobras = [("SC15", "01"), ("SC15", "55"), ("SC17", "37"),
                   ("SC29", "57"), ("SC31", "14")]

for moduleId, cobraId in invisibleCobras:
    calibrationProduct.setCobraStatus(int(cobraId), int(moduleId[2:]), invisible=True)

# Update the status information for the cobras that are known to be bad
badCobras = [("SC01", "47"), ("SC04", "22"), ("SC07", "05"), ("SC07", "19"),
             ("SC15", "23"), ("SC21", "10"), ("SC22", "13"), ("SC28", "41"),
             ("SC34", "22"), ("SC37", "01"), ("SC42", "15")]

for moduleId, cobraId in badCobras:
    calibrationProduct.setCobraStatus(int(cobraId), int(moduleId[2:]), brokenTheta=True, brokenPhi=True)

# Update the status information for the cobras that had problems measuring the
# link lengths because they hit fiducial fibers
interferencingCobras = [("SC14", "13"), ("SC27", "38"), ("SC29", "41"),
                        ("SC33", "12"), ("SC34", "01"), ("SC42", "43")]

for moduleId, cobraId in interferencingCobras:
    calibrationProduct.setCobraStatus(int(cobraId), int(moduleId[2:]), brokenTheta=True, brokenPhi=True)

# Create the bench instance
bench = Bench(layout="calibration", calibrationProduct=calibrationProduct)
print("Number of cobras:", bench.cobras.nCobras)

# Update the rMin and rMax values for the cobras that have problems so they are
# visible in the figure
rMinMedian = np.median(bench.cobras.rMin)
rMaxMedian = np.median(bench.cobras.rMax)
bench.cobras.rMin[(calibrationProduct.status & PFIDesign.COBRA_INVISIBLE_MASK) != 0] = rMinMedian
bench.cobras.rMax[(calibrationProduct.status & PFIDesign.COBRA_INVISIBLE_MASK) != 0] = rMaxMedian
bench.cobras.rMin[(calibrationProduct.status & PFIDesign.COBRA_BROKEN_MOTOR_MASK) != 0] = rMinMedian
bench.cobras.rMax[(calibrationProduct.status & PFIDesign.COBRA_BROKEN_MOTOR_MASK) != 0] = rMaxMedian

# Plot the bench
plotUtils.createNewFigure("Calibration data", "x position (mm)", "y position (mm)")
limRange = 1.05 * bench.radius * np.array([-1, 1])
xLim = bench.center.real + limRange
yLim = bench.center.imag + limRange
plotUtils.setAxesLimits(xLim, yLim)
patrolAreaColors = np.full((bench.cobras.nCobras, 4), [0.0, 0.0, 1.0, 0.15])
patrolAreaColors[(calibrationProduct.status & PFIDesign.COBRA_INVISIBLE_MASK) != 0] = [1.0, 0.0, 0.0, 0.3]
patrolAreaColors[(calibrationProduct.status & PFIDesign.COBRA_BROKEN_MOTOR_MASK) != 0] = [0.0, 1.0, 0.0, 0.5]
bench.cobras.addPatrolAreasToFigure(colors=patrolAreaColors, paintHardStops=True, paintBlackDots=False)

# Pause the execution to have time to inspect the figures
plotUtils.pauseExecution()
