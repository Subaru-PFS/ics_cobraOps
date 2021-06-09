import os
import logging
import numpy as np

from ics.cobraCharmer.pfiDesign import PFIDesign
from ics.cobraOps import plotUtils
from ics.cobraOps import targetUtils
from ics.cobraOps.Bench import Bench
from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector
from procedures.moduleTest.cobraCoach import CobraCoach
from procedures.moduleTest import engineer

logging.getLogger('matplotlib.font_manager').disabled = True

os.environ['PFS_INSTDATA_DIR'] = "/home/jgracia/github/pfs_instdata"
mod = 'ALL'
cc = CobraCoach('fpga', loadModel=False, trajectoryMode=True, rootDir="/home/jgracia/testPFI/")
cc.loadModel(version=mod, moduleVersion="final")

# Get the calibration product from cobra coach
calibrationProduct = cc.calibModel

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

# Generate the targets
targetDensity = 5.5
targets = targetUtils.generateRandomTargets(targetDensity, bench)

# Select the targets
selector = DistanceTargetSelector(bench, targets)
selector.run()
selectedTargets = selector.getSelectedTargets()

# Get the cobra angles at the target positions
cond = np.full(cc.nCobras, False)
cond[cc.goodIdx] = True
cond = np.logical_and(cond, selectedTargets.notNull)
finalPositions = cc.pfi.anglesToPositions(cc.allCobras, np.zeros(cc.nCobras), np.deg2rad(-150) - calibrationProduct.phiIn )
finalPositions[cond] = selectedTargets.positions[cond]
thetas, phis, flags = cc.pfi.positionsToAngles(cc.allCobras[cond], selectedTargets.positions[cond])
thetas = thetas[:, 0]
phis = phis[:, 0]
cobraIdx = np.where(cond)[0]

engineer.setCobraCoach(cc)

#mmTheta = np.load('/home/jgracia/testPFI/SP01_mmThetaFast5.npy')
#mmPhi = np.load('/home/jgracia/testPFI/SP01_mmPhiFast5.npy')
#mmThetaSlow = np.load('/home/jgracia/testPFI/SP01_mmThetaSlow.npy')
#mmPhiSlow = np.load('/home/jgracia/testPFI/SP01_mmPhiSlow.npy')
#engineer.setConstantSpeedMaps(mmTheta, mmPhi, mmThetaSlow, mmPhiSlow)

engineer.setConstantOntimeMode(maxSteps=3000)
#engineer.setConstantSpeedMode(maxSegments=10, maxSteps=100)

traj, moves = engineer.createTrajectory(cobraIdx, thetas, phis, tries=8, twoSteps=True, threshold=20.0, timeStep=20)

fiberPositions = traj.calculateFiberPositions(cc)
elbowPositions = traj.calculateElbowPositions(cc)

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

bench.cobras.addLinksToFigure(finalPositions, colors=linkColors)
plotUtils.addTrajectories(fiberPositions)
plotUtils.addTrajectories(elbowPositions, colors=np.array([1.0, 0.0, 0.0, 1.0]))

# Draw the targets assigned to the cobras
#selectedTargets.addToFigure(colors=np.array([0.0, 1.0, 0.0, 1.0]))
plotUtils.addPoints(targets.positions, s=3, facecolor=np.array([0.7, 0.7, 0.7, 1.0]))
plotUtils.addPoints(selectedTargets.positions[cond], s=3, facecolor=np.array([0.0, 1.0, 0.0, 1.0]))

# Pause the execution to have time to inspect the figures
plotUtils.pauseExecution()


