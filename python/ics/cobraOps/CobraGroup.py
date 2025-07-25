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

from . import plotUtils
from .AttributePrinter import AttributePrinter


class CobraGroup(AttributePrinter):
    """Class describing the properties of a group of cobras.

    """

    COBRA_LINK_RADIUS = 1.0
    """The default cobra link radius in mm."""


    def __init__(self, cobraCoach):
        """Constructs a new CobraGroup instance from the cobra coach instance.

        Parameters
        ----------
        cobraCoach: object
            The cobra coach instance containing the cobra properties.

        Returns
        -------
        object
            The CobraGroup instance.

        """
        # Save the cobra coach instance
        self.cobraCoach = cobraCoach

        # Get the calibration product from the cobra coach
        calibrationProduct = self.cobraCoach.calibModel

        # Get the number of cobras and their central positions
        self.nCobras = len(calibrationProduct.centers)
        self.centers = calibrationProduct.centers.copy()

        # Get their status information
        self.status = calibrationProduct.status.copy()
        self.isGood = np.full(self.nCobras, False)
        self.isGood[self.cobraCoach.goodIdx] = True

        # Get the theta home angles
        self.tht0 = calibrationProduct.tht0.copy()
        self.tht1 = calibrationProduct.tht1.copy()

        # Get the phi limit angles
        self.phiIn = calibrationProduct.phiIn.copy()
        self.phiOut = calibrationProduct.phiOut.copy()

        # Get the link lengths and radius
        self.L1 = calibrationProduct.L1.copy()
        self.L2 = calibrationProduct.L2.copy()
        self.linkRadius = np.full(self.nCobras, CobraGroup.COBRA_LINK_RADIUS)

        # Calculate the patrol areas minimum and maximum radii
        self.rMin = np.abs(
            self.L1 + self.L2 * np.exp(1j * np.maximum(-np.pi, self.phiIn)))
        self.rMax = np.abs(
            self.L1 + self.L2 * np.exp(1j * np.minimum(self.phiOut, 0)))

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

    def calculateElbowPositions(self, fiberPositions, indices=None):
        """Calculates the cobra elbow positions for the given fiber positions.

        The code assumes that the cobras can reach the given positions.

        Parameters
        ----------
        fiberPositions: object
            A complex numpy array with the cobra fiber positions.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all cobras will be used. Default is None.

        Returns
        -------
        object
            A complex numpy array with the cobra elbow positions.

        """
        # Select a subset of the cobras if necessary
        cobras = self.cobraCoach.allCobras

        if indices is not None:
            cobras = cobras[indices]
            fiberPositions = fiberPositions[indices]

        # Calculate the cobras theta angles at the fiber positions
        thetaAngles, _, _ = self.cobraCoach.pfi.positionsToAngles(
            cobras, fiberPositions)

        # Select the first angle solution
        thetaAngles = thetaAngles[:, 0]

        # Return the cobras elbow positions
        return self.cobraCoach.pfi.anglesToElbowPositions(cobras, thetaAngles)

    def calculateCobraElbowPositions(self, cobraIndex, fiberPositions):
        """Calculates the elbow positions for a given cobra and a list of fiber
        positions.

        The code assumes that the cobra can reach the given positions.

        Parameters
        ----------
        cobraIndex: int
            The cobra index
        fiberPositions: object
            A complex numpy array with the fiber positions for the given cobra.

        Returns
        -------
        object
            A complex numpy array with the cobra elbow positions.

        """
        # Select the cobra
        cobra = self.cobraCoach.allCobras[[cobraIndex]]

        # Calculate the cobra elbow positions at each fiber position
        elbowPositions = np.empty(fiberPositions.shape, dtype=complex)

        for i, fiberPosition in enumerate(fiberPositions):
            # Calculate the cobra theta angles at the fiber position
            thetaAngles, _, _ = self.cobraCoach.pfi.positionsToAngles(
                cobra, np.array([fiberPosition], dtype=complex))

            # Select the first angle solution
            thetaAngles = thetaAngles[:, 0]

            # Calculate the elbow position
            elbowPositions[i] = self.cobraCoach.pfi.anglesToElbowPositions(
                cobra, thetaAngles)[0]

        return elbowPositions

    def addPatrolAreasToFigure(self, colors=np.array([0.0, 0.0, 1.0, 0.15]),
                               indices=None, paintHardStops=True):
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

    def addLinksToFigure(self, fiberPositions, elbowPositions=None,
                         colors=np.array([0.0, 0.0, 1.0, 0.5]), indices=None):
        """Draws the cobras links on top of an existing figure.

        Parameters
        ----------
        fiberPositions: object
            A complex numpy array with the fiber positions.
        elbowPositions: object, optional
            A complex numpy array with the elbow positions. If it is set to
            None, the elbow positions will be calculated from the fiber
            positions. Default is None.
        colors: object, optional
            The link colors. Default is light blue.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all the cobras will be used. Default is None.

        """
        # Calculate the elbow positions if necessary
        if elbowPositions is None:
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
