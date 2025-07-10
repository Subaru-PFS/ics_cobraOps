"""

CobraGroup class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from ics.cobraCharmer.pfiDesign import PFIDesign

from . import plotUtils
from .AttributePrinter import AttributePrinter
from .cobraConstants import COBRA_LINK_RADIUS


class CobraGroup(AttributePrinter):
    """Class describing the properties of a group of cobras.

    """

    def __init__(self, calibrationProduct):
        """Constructs a new CobraGroup instance from the calibration product.

        Parameters
        ----------
        calibrationProduct: object
            The cobras calibration product containing the cobra properties.

        Returns
        -------
        object
            The CobraGroup instance.

        """
        # Get the number of cobras and their central positions
        self.nCobras = len(calibrationProduct.centers)
        self.centers = calibrationProduct.centers.copy()

        # Get their status information
        self.status = calibrationProduct.status.copy()
        self.hasProblem = self.status != PFIDesign.COBRA_OK_MASK

        # Get the theta home angles
        self.tht0 = calibrationProduct.tht0.copy()
        self.tht1 = calibrationProduct.tht1.copy()

        # Get the phi home angles
        self.phiIn = calibrationProduct.phiIn.copy()
        self.phiOut = calibrationProduct.phiOut.copy()

        # Get the link lengths and radius
        self.L1 = calibrationProduct.L1.copy()
        self.L2 = calibrationProduct.L2.copy()
        self.linkRadius = np.full(self.nCobras, COBRA_LINK_RADIUS)

        # Calculate the patrol areas minimum and maximum radii
        self.calculatePatrolAreaRadii()

        # Calculate the home positions
        self.calculateHomePositions()

    def calculatePatrolAreaRadii(self):
        """Calculates the minimum and maximum radius that the cobras can reach.

        """
        self.rMin = np.abs(
            self.L1 + self.L2 * np.exp(1j * np.maximum(-np.pi, self.phiIn)))
        self.rMax = np.abs(
            self.L1 + self.L2 * np.exp(1j * np.minimum(self.phiOut, 0)))

    def calculateHomePositions(self, phiOffset=0.00001):
        """Calculates the cobras home positions.

        Parameters
        ----------
        phiOffset: float, optional
            Cobra homes phi angle offset in radians, relative to their phiIn
            values. Default is 0.00001 radians.

        """
        # Calculate the phi home angles
        self.phiHome = np.maximum(-np.pi, self.phiIn) + phiOffset

        # Calculate the home positions (0 = same sense, 1 = opposite sense)
        self.home0 = self.calculateFiberPositions(self.tht0, self.phiHome)
        self.home1 = self.calculateFiberPositions(self.tht1, self.phiHome)

    def calculateFiberPositions(self, tht, phi, indices=None):
        """Calculates the cobras fiber positions for the given rotation angles.

        Parameters
        ----------
        tht: object
            A numpy array with the cobra theta angles in radians.
        phi: object
            A numpy array with the cobra phi angles in radians.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all cobras will be used. Default is None.

        Returns
        -------
        object
            A complex numpy array with the cobra fiber positions.

        """
        # Set the angles to the home position angles for cobras with problems
        tht = tht.copy()
        phi = phi.copy()
        tht[self.hasProblem] = self.tht0[self.hasProblem]
        phi[self.hasProblem] = self.phiHome[self.hasProblem]

        # Extract some useful information
        centers = self.centers
        L1 = self.L1
        L2 = self.L2

        # Select a subset of the cobras if necessary
        if indices is not None:
            centers = centers[indices]
            L1 = L1[indices]
            L2 = L2[indices]
            tht = tht[indices]
            phi = phi[indices]

        # Return the fiber positions
        return centers + L1 * np.exp(1j * tht) + L2 * np.exp(1j * (tht + phi))

    def calculateElbowPositions(self, fiberPositions, indices=None,
                                useNegativePhi=True):
        """Calculates the cobra elbow positions for the given fiber positions.

        The code assumes that the cobras can reach the given positions.

        Parameters
        ----------
        fiberPositions: object
            A complex numpy array with the cobra fiber positions.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all cobras will be used. Default is None.
        useNegativePhi: bool, optional
            If True the phi angle values will be negative. If False, the phi
            angles will be positive. Default is True.

        Returns
        -------
        object
            A complex numpy array with the cobra elbow positions.

        """
        # Set the fiber positions to the home position for cobras with problems
        fiberPositions = fiberPositions.copy()
        fiberPositions[self.hasProblem] = self.home0[self.hasProblem]

        # Extract some useful information
        centers = self.centers
        L1 = self.L1
        L2 = self.L2

        # Select a subset of the cobras if necessary
        if indices is not None:
            centers = centers[indices]
            L1 = L1[indices]
            L2 = L2[indices]
            fiberPositions = fiberPositions[indices]

        # Calculate the cobras theta angles applying the law of cosines
        relativePositions = fiberPositions - centers
        distance = np.abs(relativePositions)
        distanceSq = distance ** 2
        L1Sq = L1 ** 2
        L2Sq = L2 ** 2
        phiSign = -1 if useNegativePhi else +1
        tht = np.angle(relativePositions) - phiSign * np.arccos(
            -(L2Sq - L1Sq - distanceSq) / (2 * L1 * distance))

        # Return the elbow positions
        return centers + L1 * np.exp(1j * tht)

    def calculateMultipleElbowPositions(self, finalPositions, cobraIndices,
                                        targetIndices, useNegativePhi=True):
        """Calculates the cobra elbow positions for a list of cobra-target
        associations.

        The code assumes that the cobras can reach the given target positions.

        Parameters
        ----------
        finalPositions: object
            A complex numpy array with the target positions.
        cobraIndices: object
            A numpy array with the cobra indices to use.
        targetIndices: object
            A numpy array with the target indices to use.
        useNegativePhi: bool, optional
            If True the phi angle values will be negative. If False, the phi
            angles will be positive. Default is True.

        Returns
        -------
        object
            A complex numpy array with the cobra elbow positions.

        """
        # Extract some useful information
        centers = self.centers[cobraIndices]
        hasProblem = self.hasProblem[cobraIndices]
        L1 = self.L1[cobraIndices]
        L2 = self.L2[cobraIndices]
        home0 = self.home0[cobraIndices]
        finalPositions = finalPositions[targetIndices]

        # Set the target positions to the home position for cobras with problems
        finalPositions[hasProblem] = home0[hasProblem]

        # Calculate the cobras theta angles applying the law of cosines
        relativePositions = finalPositions - centers
        distance = np.abs(relativePositions)
        distanceSq = distance ** 2
        L1Sq = L1 ** 2
        L2Sq = L2 ** 2
        phiSign = -1 if useNegativePhi else +1
        tht = np.angle(relativePositions) - phiSign * np.arccos(
            -(L2Sq - L1Sq - distanceSq) / (2 * L1 * distance))

        # Return the elbow positions
        return centers + L1 * np.exp(1j * tht)

    def calculateCobraElbowPositions(self, cobraIndex, fiberPositions,
                                     useNegativePhi=True):
        """Calculates the elbow positions for a given cobra and a list of fiber
        positions.

        The code assumes that the cobra can reach the given positions.

        Parameters
        ----------
        cobraIndex: int
            The cobra index
        fiberPositions: object
            A complex numpy array with the fiber positions for the given cobra.
        useNegativePhi: bool, optional
            If True the phi angle values will be negative. If False, the phi
            angles will be positive. Default is True.

        Returns
        -------
        object
            A complex numpy array with the cobra elbow positions.

        """
        # Set the fiber positions to the home position if the cobra has problems
        if self.hasProblem[cobraIndex]:
            fiberPositions = fiberPositions.copy()
            fiberPositions[:] = self.home0[cobraIndex]

        # Extract the cobra information
        center = self.centers[cobraIndex]
        L1 = self.L1[cobraIndex]
        L2 = self.L2[cobraIndex]

        # Calculate the cobra theta angles applying the law of cosines
        relativePositions = fiberPositions - center
        distance = np.abs(relativePositions)
        distanceSq = distance ** 2
        L1Sq = L1 ** 2
        L2Sq = L2 ** 2
        phiSign = -1 if useNegativePhi else +1
        tht = np.angle(relativePositions) - phiSign * np.arccos(
            -(L2Sq - L1Sq - distanceSq) / (2 * L1 * distance))

        # Return the elbow positions
        return center + L1 * np.exp(1j * tht)

    def calculateRotationAngles(self, fiberPositions, indices=None,
                                useNegativePhi=True):
        """Calculates the cobra rotation angles for the given fiber positions.

        The code assumes that the cobras can reach the given positions.

        Parameters
        ----------
        fiberPositions: object
            A complex numpy array with the fiber positions.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all the cobras will be used. Default is None.
        useNegativePhi: bool, optional
            If True the phi angle values will be negative. If False, the phi
            angles will be positive. Default is True.

        Returns
        -------
        tuple
            A python tuple with the cobra rotation angles (theta, phi).

        """
        # Set the fiber positions to the home position for cobras with problems
        fiberPositions = fiberPositions.copy()
        fiberPositions[self.hasProblem] = self.home0[self.hasProblem]

        # Extract some useful information
        centers = self.centers
        L1 = self.L1
        L2 = self.L2

        # Select a subset of the cobras if necessary
        if indices is not None:
            centers = centers[indices]
            L1 = L1[indices]
            L2 = L2[indices]
            fiberPositions = fiberPositions[indices]

        # Calculate the cobras rotation angles applying the law of cosines
        relativePositions = fiberPositions - centers
        distance = np.abs(relativePositions)
        distanceSq = distance ** 2
        L1Sq = L1 ** 2
        L2Sq = L2 ** 2
        phiSign = -1 if useNegativePhi else +1
        phi = phiSign * np.arccos((distanceSq - L1Sq - L2Sq) / (2 * L1 * L2))
        tht = np.angle(relativePositions) - phiSign * np.arccos(
            -(L2Sq - L1Sq - distanceSq) / (2 * L1 * distance))

        # Force tht to go from -pi to pi
        tht = (tht - np.pi) % (2 * np.pi) - np.pi

        return (tht, phi)

    def addPatrolAreasToFigure(self, colors=np.array([0.0, 0.0, 1.0, 0.15]),
                               indices=None, paintHardStops=True,
                               paintBlackDots=True):
        """Draws the cobra patrol areas on top of an existing figure.

        Parameters
        ----------
        colors: object, optional
            The patrol area colors. Default is very light blue.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all the cobras will be used. Default is None.
        paintHardStops: bool, optional
            True if the cobra hard stops should be painted. Default is True.

        """
        # Extract some useful information
        centers = self.centers
        rMin = self.rMin
        rMax = self.rMax
        tht0 = self.tht0
        tht1 = self.tht1

        # Select a subset of the cobras if necessary
        if indices is not None:
            centers = centers[indices]
            rMin = rMin[indices]
            rMax = rMax[indices]
            tht0 = tht0[indices]
            tht1 = tht1[indices]

            if colors.ndim == 2:
                colors = colors[indices]

        # Draw the cobra patrol areas using ring shapes
        plotUtils.addRings(centers, rMin, rMax, facecolors=colors)

        # Add the theta hard stops if necessary
        if paintHardStops:
            plotUtils.addLines(
                centers, centers + rMax * np.exp(1j * tht0),
                linewidths=1, linestyles="dashed", color="0.3")
            plotUtils.addLines(
                centers, centers + rMax * np.exp(1j * tht1),
                linewidths=1, linestyles="dashdot", color="0.3")

    def addLinksToFigure(self, fiberPositions,
                         colors=np.array([0.0, 0.0, 1.0, 0.5]), indices=None):
        """Draws the cobras links on top of an existing figure.

        Parameters
        ----------
        fiberPositions: object
            A complex numpy array with the fiber positions.
        colors: object, optional
            The link colors. Default is light blue.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all the cobras will be used. Default is None.

        """
        # Set the fiber positions to the home position for cobras with problems
        fiberPositions = fiberPositions.copy()
        fiberPositions[self.hasProblem] = self.home0[self.hasProblem]

        # Calculate the elbow positions
        elbowPositions = self.calculateElbowPositions(
            fiberPositions, indices=indices)

        # Extract some useful information
        centers = self.centers
        linkRadius = self.linkRadius

        # Select a subset of the cobras if necessary
        if indices is not None:
            centers = centers[indices]
            linkRadius = linkRadius[indices]
            fiberPositions = fiberPositions[indices]

            if colors.ndim == 2:
                colors = colors[indices]

        # Draw the cobras using a combination of thin and thick lines
        plotUtils.addLines(
            centers, elbowPositions, edgecolor=colors, linewidths=2)
        plotUtils.addThickLines(
            elbowPositions, fiberPositions, linkRadius, facecolors=colors)
