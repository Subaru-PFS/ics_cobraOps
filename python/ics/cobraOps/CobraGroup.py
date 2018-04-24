"""

Cobra group class.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

from . import plotUtils

from .AttributePrinter import AttributePrinter
from .MotorMapGroup import MotorMapGroup
from .cobraConstants import (COBRA_LINK_LENGTH,
                                         COBRA_LINK_RADIUS,
                                         PHI_SAFETY_ANGLE,
                                         HOMES_THETA_DISTANCE,
                                         BLACK_DOT_RELATIVE_POSITION,
                                         BLACK_DOT_RADIUS)


class CobraGroup(AttributePrinter):
    """

    Class describing the properties of a group of cobras.

    """

    def __init__(self, centers):
        """Constructs a new cobra group instance with default properties.

        Parameters
        ----------
        centers: object
            A complex numpy array with the cobras central positions.

        Returns
        -------
        object
            The cobra group instance.

        """
        # Set the number of cobras and their central positions
        self.nCobras = len(centers)
        self.centers = centers.copy()

        # Set the theta home angles randomly
        self.tht0 = 2 * np.pi * np.random.random(self.nCobras)
        self.tht1 = (self.tht0 + HOMES_THETA_DISTANCE) % (2 * np.pi)

        # Set the phi home angles randomly
        self.phiOut = -PHI_SAFETY_ANGLE * (1.0 + 0.2 * np.random.random(self.nCobras))
        self.phiIn = self.phiOut - (np.pi - 2 * PHI_SAFETY_ANGLE)

        # Set the default link lengths
        self.L1 = np.full(self.nCobras, COBRA_LINK_LENGTH)
        self.L2 = self.L1.copy()

        # Set the default link radius
        self.linkRadius = np.full(self.nCobras, COBRA_LINK_RADIUS)

        # Set the default black dot position and radius
        self.blackDotPosition = np.full(self.nCobras, BLACK_DOT_RELATIVE_POSITION)
        self.blackDotRadius = np.full(self.nCobras, BLACK_DOT_RADIUS)

        # Set the default motor maps
        self.motorMaps = MotorMapGroup(self.nCobras)

        # Calculate the patrol areas minimum and maximum radii
        self.calculatePatrolAreaRadii()

        # Calculate the cobras home positions
        self.calculateHomePositions()


    def calculatePatrolAreaRadii(self):
        """Calculates the minimum and maximum radius that the cobras can reach.

        """
        self.rMin = np.abs(self.L1 + self.L2 * np.exp(1j * self.phiIn))
        self.rMax = np.abs(self.L1 + self.L2 * np.exp(1j * self.phiOut))


    def calculateHomePositions(self, phiOffset=0.5):
        """Calculates the cobras home positions.

        Parameters
        ----------
        phiOffset: float, optional
            Cobra homes phi angle offset in radians, relative to their phiIn
            values. Default is 0.5 radians (28.6 degrees).

        """
        # Calculate the phi home angles
        phiHome = self.phiIn + phiOffset

        # Calculate the home positions (0 = same sense, 1 = opposite sense)
        self.home0 = self.calculateFiberPositions(self.tht0, phiHome)
        self.home1 = self.calculateFiberPositions(self.tht1, phiHome)


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
            all the cobras will be used. Default is None.

        Returns
        -------
        object
            A complex numpy array with the fiber positions.

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


    def calculateElbowPositions(self, fiberPositions, indices=None, useNegativePhi=True):
        """Calculates the cobra elbow positions for the given fiber positions.

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
        object
            A complex numpy array with the elbow positions.

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
            fiberPositions = fiberPositions[indices]

        # Calculate the cobras theta angles applying the law of cosines
        relativePositions = fiberPositions - centers
        distance = np.abs(relativePositions)
        distanceSq = distance ** 2
        L1Sq = L1 ** 2
        L2Sq = L2 ** 2
        phiSign = -1 if useNegativePhi else +1
        tht = np.angle(relativePositions) - phiSign * np.arccos(-(L2Sq - L1Sq - distanceSq) / (2 * L1 * distance))

        # Return the elbow positions
        return centers + L1 * np.exp(1j * tht)


    def calculateCobraElbowPositions(self, cobraIndex, fiberPositions, useNegativePhi=True):
        """Calculates the elbow positions for a given cobra and a list of fiber positions.

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
        tht = np.angle(relativePositions) - phiSign * np.arccos(-(L2Sq - L1Sq - distanceSq) / (2 * L1 * distance))

        # Return the elbow positions
        return center + L1 * np.exp(1j * tht)


    def calculateRotationAngles(self, fiberPositions, indices=None, useNegativePhi=True):
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
            A python tuple with the cobras rotation angles (theta, phi).

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
            fiberPositions = fiberPositions[indices]

        # Calculate the cobras rotation angles applying the law of cosines
        relativePositions = fiberPositions - centers
        distance = np.abs(relativePositions)
        distanceSq = distance ** 2
        L1Sq = L1 ** 2
        L2Sq = L2 ** 2
        phiSign = -1 if useNegativePhi else +1
        phi = phiSign * np.arccos((distanceSq - L1Sq - L2Sq) / (2 * L1 * L2))
        tht = np.angle(relativePositions) - phiSign * np.arccos(-(L2Sq - L1Sq - distanceSq) / (2 * L1 * distance))

        # Force tht to go from -pi to pi
        tht = (tht - np.pi) % (2 * np.pi) - np.pi

        return (tht, phi)


    def useCalibrationProduct(self, calibrationProduct, useRealLinks=True, useRealMaps=True):
        """Updates the cobra properties with the calibration product ones.

        Parameters
        ----------
        calibrationProduct: object
            The cobras calibration product containing the cobra properties.
        useRealLinks: bool, optional
            If True, the cobras link properties will be updated with the
            calibration product values. Default is True.
        useRealMaps: bool, optional
            If True, the cobras motor maps will be updated with the calibration
            product values. Default is True.

        """
        # Check if we should use the calibration product link properties
        if useRealLinks:
            if calibrationProduct.nCobras == self.nCobras:
                # Use directly the calibration product arrays
                self.tht0 = calibrationProduct.tht0.copy()
                self.tht1 = calibrationProduct.tht1.copy()
                self.phiIn = calibrationProduct.phiIn.copy()
                self.phiOut = calibrationProduct.phiOut.copy()
                self.L1 = calibrationProduct.L1.copy()
                self.L2 = calibrationProduct.L2.copy()
            else:
                # Randomize the calibration cobra indices
                indices = np.random.randint(calibrationProduct.nCobras, size=(3, self.nCobras))

                # Assign random link properties to each cobra
                self.tht0 = calibrationProduct.tht0[indices[0]]
                self.tht1 = calibrationProduct.tht1[indices[0]]
                self.phiIn = calibrationProduct.phiIn[indices[1]]
                self.phiOut = calibrationProduct.phiOut[indices[1]]
                self.L1 = calibrationProduct.L1[indices[2]]
                self.L2 = calibrationProduct.L2[indices[2]]

            # Add some safety range to the phi limits
            self.phiIn += PHI_SAFETY_ANGLE
            self.phiOut -= PHI_SAFETY_ANGLE

            # Update the patrol areas minimum and maximum radii
            self.calculatePatrolAreaRadii()

            # Update the cobras home positions
            self.calculateHomePositions()

        # Check if we should use the calibration product motor maps
        if useRealMaps:
            self.motorMaps.useCalibrationProduct(calibrationProduct)


    def addPatrolAreasToFigure(self, colors=np.array([0.0, 0.0, 1.0, 0.15]), indices=None, paintHardStops=True, paintBlackDots=True):
        """Draws the cobras patrol areas on top of an existing figure.

        Parameters
        ----------
        colors: object, optional
            The patrol area colors. Default is very light blue.
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all the cobras will be used. Default is None.
        paintHardStops: bool, optional
            True if the cobra hard stops should be painted. Default is True.
        paintBlackDots: bool, optional
            True if the cobra black dots should be painted. Default is True.

        """
        # Extract some useful information
        centers = self.centers
        rMin = self.rMin
        rMax = self.rMax
        tht0 = self.tht0
        tht1 = self.tht1
        blackDotPosition = self.blackDotPosition
        blackDotRadius = self.blackDotRadius
        
        # Select a subset of the cobras if necessary
        if indices is not None:
            centers = centers[indices]
            rMin = rMin[indices]
            rMax = rMax[indices]
            tht0 = tht0[indices]
            tht1 = tht1[indices]
            blackDotPosition = blackDotPosition[indices]
            blackDotRadius = blackDotRadius[indices]

            if colors.ndim == 2:
                colors = colors[indices]

        # Draw the cobra patrol areas using ring shapes
        plotUtils.addRings(centers, rMin, rMax, facecolors=colors)
        
        # Draw the cobra black dots if necessary
        if paintBlackDots:
            plotUtils.addCircles(centers + blackDotPosition, blackDotRadius, facecolors=np.array([0.0, 0.0, 0.0, 0.15]))

        # Add the theta hard stops if necessary
        if paintHardStops:
            plotUtils.addLines(centers, centers + rMax * np.exp(1j * tht0), linewidths=1, linestyles="dashed", color="0.3")
            plotUtils.addLines(centers, centers + rMax * np.exp(1j * tht1), linewidths=1, linestyles="dashdot", color="0.3")


    def addLinksToFigure(self, fiberPositions, colors=np.array([0.0, 0.0, 1.0, 0.5]), indices=None):
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
        # Calculate the elbow positions
        elbowPositions = self.calculateElbowPositions(fiberPositions, indices=indices)

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
        plotUtils.addLines(centers, elbowPositions, edgecolor=colors, linewidths=2)
        plotUtils.addThickLines(elbowPositions, fiberPositions, linkRadius, facecolors=colors)
