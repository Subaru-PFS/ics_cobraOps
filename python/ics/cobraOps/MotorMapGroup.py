"""

MotorMapGroup class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from . import plotUtils
from .cobraConstants import HOMES_THETA_DISTANCE
from .cobraConstants import MOTOR_MAP_ANGULAR_STEP
from .cobraConstants import MOTOR1_STEP_SIZE
from .cobraConstants import MOTOR2_STEP_SIZE


class MotorMapGroup():
    """Class describing the properties of a group of cobra motor maps.

    """

    def __init__(self, nMaps):
        """Constructs a new MotorMapGroup instance with default properties.

        Parameters
        ----------
        nMaps: int
            The number of motor maps in the group. Should be the same as the
            number or cobras.

        Returns
        -------
        object
            The MotorMapGroup instance.

        """
        # Save the number of motor maps
        self.nMaps = nMaps

        # Set the default angular step sizes between measurements
        self.angularSteps = np.full(self.nMaps, MOTOR_MAP_ANGULAR_STEP)

        # Set the default step maps
        nTht = int(np.ceil(HOMES_THETA_DISTANCE / MOTOR_MAP_ANGULAR_STEP))
        nPhi = int(np.ceil(np.pi / MOTOR_MAP_ANGULAR_STEP))
        self.S1Pm = np.full((self.nMaps, nTht), np.rad2deg(MOTOR_MAP_ANGULAR_STEP) / MOTOR1_STEP_SIZE)
        self.S1Nm = np.full((self.nMaps, nTht), np.rad2deg(MOTOR_MAP_ANGULAR_STEP) / MOTOR1_STEP_SIZE)
        self.S2Pm = np.full((self.nMaps, nPhi), np.rad2deg(MOTOR_MAP_ANGULAR_STEP) / MOTOR2_STEP_SIZE)
        self.S2Nm = np.full((self.nMaps, nPhi), np.rad2deg(MOTOR_MAP_ANGULAR_STEP) / MOTOR2_STEP_SIZE)
        self.F1Pm = self.S1Pm.copy()
        self.F1Nm = self.S1Nm.copy()
        self.F2Pm = self.S2Pm.copy()
        self.F2Nm = self.S2Nm.copy()

        # Set the theta and phi offset arrays
        self.thtOffsets = np.arange(nTht + 1) * self.angularSteps[:, np.newaxis]
        self.phiOffsets = np.arange(nPhi + 1) * self.angularSteps[:, np.newaxis]

        # Calculate the integrated step maps
        self.calculateIntegratedStepMaps()

    def calculateIntegratedStepMaps(self, useSlowMaps=True):
        """Calculates the integrated theta and phi step maps.

        Parameters
        ----------
        useSlowMaps: bool, optional
            If True (False), the slow (fast) motor maps will be used in the
            calculation. Default is True.

        """
        # Calculate the cumulative step sums
        zeros = np.zeros((self.nMaps, 1))

        if useSlowMaps:
            self.posThtSteps = np.hstack((zeros, np.cumsum(self.S1Pm, axis=1)))
            self.negThtSteps = np.hstack((zeros, np.cumsum(self.S1Nm, axis=1)))
            self.posPhiSteps = np.hstack((zeros, np.cumsum(self.S2Pm, axis=1)))
            self.negPhiSteps = np.hstack((zeros, np.cumsum(self.S2Nm, axis=1)))
        else:
            self.posThtSteps = np.hstack((zeros, np.cumsum(self.F1Pm, axis=1)))
            self.negThtSteps = np.hstack((zeros, np.cumsum(self.F1Nm, axis=1)))
            self.posPhiSteps = np.hstack((zeros, np.cumsum(self.F2Pm, axis=1)))
            self.negPhiSteps = np.hstack((zeros, np.cumsum(self.F2Nm, axis=1)))

    def useCalibrationProduct(self, calibrationProduct):
        """Updates the motor map properties with the calibration product ones.

        Parameters
        ----------
        calibrationProduct: object
            The cobras calibration product containing the motor map properties.

        """
        if calibrationProduct.nCobras == self.nMaps:
            # Use directly the calibration product arrays
            self.angularSteps = calibrationProduct.angularSteps.copy()
            self.S1Pm = calibrationProduct.S1Pm.copy()
            self.S2Pm = calibrationProduct.S2Pm.copy()
            self.S1Nm = calibrationProduct.S1Nm.copy()
            self.S2Nm = calibrationProduct.S2Nm.copy()
            self.F1Pm = calibrationProduct.F1Pm.copy()
            self.F2Pm = calibrationProduct.F2Pm.copy()
            self.F1Nm = calibrationProduct.F1Nm.copy()
            self.F2Nm = calibrationProduct.F2Nm.copy()
        else:
            # Assign the motor maps properties randomly
            indices = np.random.randint(calibrationProduct.nCobras, size=self.nMaps)
            self.angularSteps = calibrationProduct.angularSteps[indices]
            self.S1Pm = calibrationProduct.S1Pm[indices]
            self.S2Pm = calibrationProduct.S2Pm[indices]
            self.S1Nm = calibrationProduct.S1Nm[indices]
            self.S2Nm = calibrationProduct.S2Nm[indices]
            self.F1Pm = calibrationProduct.F1Pm[indices]
            self.F2Pm = calibrationProduct.F2Pm[indices]
            self.F1Nm = calibrationProduct.F1Nm[indices]
            self.F2Nm = calibrationProduct.F2Nm[indices]

        # Set the theta and phi offset arrays
        self.thtOffsets = np.arange(self.S1Pm.shape[1] + 1) * self.angularSteps[:, np.newaxis]
        self.phiOffsets = np.arange(self.S2Pm.shape[1] + 1) * self.angularSteps[:, np.newaxis]

        # Update the integrated step maps
        self.calculateIntegratedStepMaps()

    def calculateSteps(self, deltaTht, startPhi, deltaPhi):
        """Calculates the total number of motor steps required to move the
        cobra fibers the given theta and phi delta angles.

        Parameters
        ----------
        deltaTht: object
            A numpy array with the theta delta offsets relative to the home
            positions.
        startPhi: object
            A numpy array with the starting phi angle positions.
        deltaPhi: object
            A numpy array with the phi delta offsets relative to the starting
            phi positions.

        Returns
        -------
        tuple
            A python tuple with the total number of motor steps for the theta
            and phi angles.

        """
        # Get the integrated step maps for the theta angle
        thtSteps = self.negThtSteps.copy()
        thtSteps[deltaTht >= 0] = self.posThtSteps[deltaTht >= 0]

        # Get the integrated step maps for the phi angle
        phiSteps = self.negPhiSteps.copy()
        phiSteps[deltaPhi >= 0] = self.posPhiSteps[deltaPhi >= 0]

        # Calculate the theta and phi offsets relative to the home positions
        thtOffset = np.zeros(self.nMaps)
        phiOffset = np.abs(startPhi)
        phiOffset[deltaPhi >= 0] = np.pi + startPhi[deltaPhi >= 0]

        # Calculate the total number of motor steps for each angle
        nThtSteps = np.empty(self.nMaps)
        nPhiSteps = np.empty(self.nMaps)

        for c in range(self.nMaps):
            # Calculate the total number of motor steps for the theta movement
            stepsRange = np.interp([thtOffset[c], thtOffset[c] + np.abs(deltaTht[c])], self.thtOffsets[c], thtSteps[c])
            nThtSteps[c] = stepsRange[1] - stepsRange[0]

            # Calculate the total number of motor steps for the phi movement
            stepsRange = np.interp([phiOffset[c], phiOffset[c] + np.abs(deltaPhi[c])], self.phiOffsets[c], phiSteps[c])
            nPhiSteps[c] = stepsRange[1] - stepsRange[0]

        return (nThtSteps, nPhiSteps)

    def plot(self, useSlowMaps=True, indices=None):
        """Plots the cobras motor maps on a new figure.

        Parameters
        ----------
        useSlowMaps: bool, optional
            If True (False), the slow (fast) motor maps will be plotted.
            Default is True.
        indices: object, optional
            A numpy array with the map indices to use. If it is set to None,
            all the maps will be used. Default is None.

        """
        # Use all the maps if the indices parameter was not provided
        if indices is None:
            indices = np.arange(self.nMaps)

        # Create a new figure
        plotUtils.createNewFigure("Cobra motor maps", "angle offsets (degrees)", " 1/speed (angle bin/ motor steps)", aspectRatio="auto")

        # Loop over all the map indices
        for m in indices:
            # Draw the theta motor maps
            angleOffsets = np.rad2deg(self.thtOffsets[m, :-1])

            if useSlowMaps:
                plotUtils.addLine(angleOffsets, self.S1Pm[m], color=[0, 0, 0, 0.4])
                plotUtils.addLine(angleOffsets, self.S1Nm[m], color=[1, 0, 0, 0.4])
            else:
                plotUtils.addLine(angleOffsets, self.F1Pm[m], color=[0, 0, 0, 0.4])
                plotUtils.addLine(angleOffsets, self.F1Nm[m], color=[1, 0, 0, 0.4])

            # Draw the phi motor maps
            angleOffsets = np.rad2deg(self.phiOffsets[m, :-1])

            if useSlowMaps:
                plotUtils.addLine(angleOffsets, self.S2Pm[m], color=[0, 1, 0, 0.4])
                plotUtils.addLine(angleOffsets, self.S2Nm[m], color=[0, 0, 1, 0.4])
            else:
                plotUtils.addLine(angleOffsets, self.F2Pm[m], color=[0, 1, 0, 0.4])
                plotUtils.addLine(angleOffsets, self.F2Nm[m], color=[0, 0, 1, 0.4])
