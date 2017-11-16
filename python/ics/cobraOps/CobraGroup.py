"""

Cobra group class.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

import ics.cobraOps.plotUtils as plotUtils

from ics.cobraOps.AttributePrinter import AttributePrinter
from ics.cobraOps.cobraConstants import (COBRA_LINK_LENGTH,
                                         COBRA_LINK_RADIUS,
                                         PHI_SAFETY_ANGLE,
                                         HOMES_THETA_DISTANCE)


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
        
        # Leave the motor maps empty
        self.angularStep = None
        self.S1Pm = None
        self.S1Nm = None
        self.S2Pm = None
        self.S2Nm = None
        
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
            A python tuple with the cobras rotation angles (tht, phi).
        
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
        
        # Force tht to go from -pi to pi, instead of from 0 to 2pi
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
                self.phiIn = calibrationProduct.phiIn.copy()
                self.phiOut = calibrationProduct.phiOut.copy()
                self.L1 = calibrationProduct.L1.copy()
                self.L2 = calibrationProduct.L2.copy()
            else:
                # Randomize the calibration cobra indices
                indices = np.random.randint(calibrationProduct.nCobras, size=(4, self.nCobras))
                
                # Assign random link properties to each cobra
                self.phiIn = calibrationProduct.phiIn[indices[0]] 
                self.phiOut = calibrationProduct.phiOut[indices[1]] 
                self.L1 = calibrationProduct.L1[indices[2]]
                self.L2 = calibrationProduct.L2[indices[3]]
            
            # Add some safety range to the phi limits
            self.phiIn += PHI_SAFETY_ANGLE
            self.phiOut -= PHI_SAFETY_ANGLE
            
            # Update the patrol areas minimum and maximum radii
            self.calculatePatrolAreaRadii()
            
            # Update the cobras home positions
            self.calculateHomePositions()
        
        # Check if we should use the calibration product motor maps
        if useRealMaps:
            if calibrationProduct.nCobras == self.nCobras:
                # Use directly the calibration product arrays
                self.angularStep = calibrationProduct.angularStep.copy()
                self.S1Pm = calibrationProduct.S1Pm.copy()
                self.S2Pm = calibrationProduct.S2Pm.copy()
                self.S1Nm = calibrationProduct.S1Nm.copy()
                self.S2Nm = calibrationProduct.S2Nm.copy()
            else:
                # Randomize the calibration cobra indices
                indices = np.random.randint(calibrationProduct.nCobras, size=(5, self.nCobras))
                
                # Assign a random motor map to each cobra
                self.angularStep = calibrationProduct.angularStep[indices[0]]
                self.S1Pm = calibrationProduct.S1Pm[indices[1]]
                self.S2Pm = calibrationProduct.S2Pm[indices[2]]
                self.S1Nm = calibrationProduct.S1Nm[indices[3]]
                self.S2Nm = calibrationProduct.S2Nm[indices[4]]
    
    
    def addPatrolAreasToFigure(self, colors=[0.0, 0.0, 1.0, 0.15], indices=None, paintHardStops=True):
        """Draws the cobras patrol areas on top off an existing figure.
        
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
        
        # Add the stage 1 theta hard stops if necessary
        if paintHardStops:
            plotUtils.addLines(centers, centers + rMax * np.exp(1j * tht0), linewidths=1, linestyles="dashed", color="0.3")
            plotUtils.addLines(centers, centers + rMax * np.exp(1j * tht1), linewidths=1, linestyles="dashdot", color="0.3")
    
    
    def addLinksToFigure(self, fiberPositions, colors=[0.0, 0.0, 1.0, 0.5], indices=None):
        """Draws the cobras links on top off an existing figure.
        
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
    
    
    def addMotorMapsToFigure(self, indices=None):
        """Draws the cobras motor maps on top off an existing figure.
        
        Parameters
        ----------
        indices: object, optional
            A numpy array with the cobra indices to use. If it is set to None,
            all the cobras will be used. Default is None.
        
        """
        # Return if the motors maps are not defined
        if self.S1Pm is None:
            return
        
        # Use all the cobras if the indices parameter was not provided
        if indices is None:
            indices = np.arange(self.nCobras)
        
        # Loop over all the cobra indices
        for c in indices:
            # Draw the theta motor maps
            angles = np.arange(0, 360, self.angularStep[c])
            plotUtils.addLine(angles, self.S1Pm[c, :len(angles)], color=[0, 0, 0, 0.4])
            plotUtils.addLine(angles, self.S1Nm[c, :len(angles)], color=[1, 0, 0, 0.4])
            
            # Draw the phi motor maps
            angles = np.arange(0, 180, self.angularStep[c])
            plotUtils.addLine(angles, self.S2Pm[c, :len(angles)], color=[0, 1, 0, 0.4])
            plotUtils.addLine(angles, self.S2Nm[c, :len(angles)], color=[0, 0, 1, 0.4])
