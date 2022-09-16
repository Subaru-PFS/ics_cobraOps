"""

Bench class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from .BlackDotGroup import BlackDotGroup
from .cobraConstants import COBRAS_SEPARATION
from .cobraConstants import MODULE_FIRST_LINE_LENGTH
from .cobraConstants import MODULE_SECOND_LINE_LENGTH
from .cobraConstants import MODULES_PER_SECTOR
from .CobraGroup import CobraGroup


class Bench:
    """Class describing the properties of a PFI bench.

    """

    def __init__(self, cobraCenters=None, layout="calibration",
                 calibrationProduct=None, blackDotsCalibrationProduct=None,
                 blackDotsMargin=1.0
                 ):
        """Constructs a new Bench instance.

        Parameters
        ----------
        cobraCenters: object, optional
            A complex numpy array with the cobras central positions. If it is
            set to None, the cobra centers will follow the specified bench
            layout. Default is None.
        layout: str, optional
            The bench layout to use: full, hex, line, rails or calibration.
            This parameter is only used if the cobraCenters parameter is set to
            None. Default is "calibration".
        calibrationProduct: object, optional
            The cobras calibration product with the cobra properties to use.
            Default is None.
        blackDotsCalibrationProduct: object, optional
            The black dots calibration product with the black dots properties
            to use. Default is None.
        blackDotsMargin: real, optional
            The margin factor in radius of the black dots to avoid in fiber allocation. Default is 1.0.

        Returns
        -------
        object
            The Bench instance.

        """
        # Check that all the calibration products are provided when we use the
        # calibration mode
        if layout == "calibration":
            if calibrationProduct is None:
                raise Exception(
                    "The cobras calibration product needs to be provided")

            if blackDotsCalibrationProduct is None:
                raise Exception(
                    "The black dots calibration product needs to be provided")

        # Calculate the cobra centers if they have not been provided
        if cobraCenters is None:
            # Check if the centers from the calibration product should be used
            if layout == "calibration":
                cobraCenters = calibrationProduct.centers.copy()
            else:
                cobraCenters = Bench.calculateCobraCenters(layout)

        # Create the cobra group instance
        self.cobras = CobraGroup(cobraCenters)

        # Update the cobra properties if necessary
        if calibrationProduct is not None:
            self.cobras.useCalibrationProduct(calibrationProduct)

        # Calculate the bench center
        self.center = np.mean(self.cobras.centers[~self.cobras.hasProblem])

        # Calculate the bench radius
        self.radius = np.max(np.abs(
            self.cobras.centers[~self.cobras.hasProblem] - self.center) + 
            self.cobras.rMax[~self.cobras.hasProblem])

        # Calculate the cobra nearest neighbors associations array
        self.calculateCobraAssociations()

        # Create the black dot group instance
        self.blackDots = BlackDotGroup(self.cobras.centers)

        # Update the black dots properties if necessary
        if blackDotsCalibrationProduct is not None:
            self.blackDots.useCalibrationProduct(blackDotsCalibrationProduct)

        # Apply margin factor in radius for the black dot avoidance
        self.blackDots.radius *= blackDotsMargin

    def calculateCobraAssociations(self):
        """Calculates the cobras nearest neighbors associations array.

        Each column in the array contains a different cobra association.

        """
        # Calculate the cobras distance matrix
        distanceMatrix = np.abs(
            self.cobras.centers[:, np.newaxis] - self.cobras.centers)

        # Mask the diagonal because it contains distances to the same cobras
        distanceMatrix[np.arange(self.cobras.nCobras),
                       np.arange(self.cobras.nCobras)] = np.Inf

        # Calculate the median minimum distance between cobras
        medianMinDistance = np.median(
            np.min(distanceMatrix, axis=1)[~self.cobras.hasProblem])

        # Obtain the nearest neighbors indices
        (cobrasIndices, nearbyCobrasIndices) = np.where(
            distanceMatrix < 1.5 * medianMinDistance)

        # Remove all the duplicated cobra associations
        uniqueAssociations = cobrasIndices < nearbyCobrasIndices
        cobrasIndices = cobrasIndices[uniqueAssociations]
        nearbyCobrasIndices = nearbyCobrasIndices[uniqueAssociations]

        # Save the cobra associations in a single array
        self.cobraAssociations = np.vstack((cobrasIndices, nearbyCobrasIndices))

    def getCobraNeighbors(self, cobraIndex):
        """Returns the indices of the cobras that are neighbors to a given
        cobra.

        Parameters
        ----------
        cobraIndex: int
            The cobra index.

        Returns
        -------
        object
            A numpy array with the cobra neighbor indices.

        """
        # Get the two cobra neighbor association possibilities
        firstNeighborGroup = self.cobraAssociations[
            1][self.cobraAssociations[0] == cobraIndex]
        secondNeighborGroup = self.cobraAssociations[
            0][self.cobraAssociations[1] == cobraIndex]

        return np.concatenate((firstNeighborGroup, secondNeighborGroup))

    def getCobrasNeighbors(self, cobraIndices):
        """Returns the indices of the cobras that are neighbors to the given
        cobras.

        Parameters
        ----------
        cobraIndices: object
            A numpy array with the cobra indices.

        Returns
        -------
        object
            A numpy array with the cobras neighbor indices.

        """
        # Get the two neighbor association possibilities
        firstNeighborGroup = self.cobraAssociations[
            1][np.in1d(self.cobraAssociations[0], cobraIndices)]
        secondNeighborGroup = self.cobraAssociations[
            0][np.in1d(self.cobraAssociations[1], cobraIndices)]

        return np.unique(
            np.concatenate((firstNeighborGroup, secondNeighborGroup)))

    def getCollisionsForCobra(self, cobraIndex, fiberPositions):
        """Calculates the total number of collisions for a given cobra.

        Parameters
        ----------
        cobraIndex: int
            The cobra index.
        fiberPositions: object
            A complex numpy array with the cobras fiber positions.

        Returns
        -------
        int
            The number of collisions between the cobra and its neighbors.

        """
        # Get the indices of the nearby cobras
        nearbyCobrasIndices = self.getCobraNeighbors(cobraIndex)
        nNearbyCobras = len(nearbyCobrasIndices)

        # Set the fiber positions to the home position for cobras with problems
        fiberPositions = fiberPositions.copy()
        fiberPositions[self.cobras.hasProblem] = self.cobras.home0[
            self.cobras.hasProblem]

        # Calculate the cobras elbow positions
        allIndices = np.append(cobraIndex, nearbyCobrasIndices)
        elbowPositions = self.cobras.calculateElbowPositions(
            fiberPositions, indices=allIndices)

        # Calculate the distances between the cobra link and the nearby cobras
        # links
        startPoints1 = np.repeat(fiberPositions[cobraIndex], nNearbyCobras)
        endPoints1 = np.repeat(elbowPositions[0], nNearbyCobras)
        startPoints2 = fiberPositions[nearbyCobrasIndices]
        endPoints2 = elbowPositions[1:]
        distances = Bench.distancesBetweenLineSegments(
            startPoints1, endPoints1, startPoints2, endPoints2)

        # Get the cobra collisions for the current configuration
        collisions = distances < (
            self.cobras.linkRadius[cobraIndex] + self.cobras.linkRadius[
                nearbyCobrasIndices])

        return np.sum(collisions)

    def calculateCobraAssociationCollisions(self, fiberPositions,
                                            associationIndices=None):
        """Calculates which cobra associations are involved in a collision.

        Parameters
        ----------
        fiberPositions: object
            A complex numpy array with the cobras fiber positions.
        associationIndices: object, optional
            A numpy array with the cobra association indices to use. If it is
            set to None, all the cobra associations will be used. Default is
            None.

        Returns
        -------
        object
            A boolean numpy array indicating which cobra associations are
            involved in a collision at the given fiber positions.

        """
        # Extract some useful information
        cobraAssociations = self.cobraAssociations
        linkRadius = self.cobras.linkRadius

        # Select a subset of the cobra associations if necessary
        if associationIndices is not None:
            cobraAssociations = cobraAssociations[:, associationIndices]

        # Set the fiber positions to the home position for cobras with problems
        fiberPositions = fiberPositions.copy()
        fiberPositions[self.cobras.hasProblem] = self.cobras.home0[
            self.cobras.hasProblem]

        # Calculate the cobras elbow positions
        elbowPositions = self.cobras.calculateElbowPositions(fiberPositions)

        # Calculate the distances between the cobras links
        startPoints1 = fiberPositions[cobraAssociations[0]]
        endPoints1 = elbowPositions[cobraAssociations[0]]
        startPoints2 = fiberPositions[cobraAssociations[1]]
        endPoints2 = elbowPositions[cobraAssociations[1]]
        distances = Bench.distancesBetweenLineSegments(
            startPoints1, endPoints1, startPoints2, endPoints2)

        # Return the cobra associations collisions
        return distances < (
            linkRadius[cobraAssociations[0]] + linkRadius[cobraAssociations[1]])

    def getProblematicCobraAssociations(self, fiberPositions):
        """Returns the indices of the cobra associations involved in a
        collision.

        Parameters
        ----------
        fiberPositions: object
            A complex numpy array with the cobras fiber positions.

        Returns
        -------
        object
            An integer numpy array with the indices of the cobras involved in a
            collision. Each column in the array contains a different cobra
            association.

        """
        # Calculate the cobra association collisions
        collisions = self.calculateCobraAssociationCollisions(fiberPositions)

        # Return the indices of the cobras involved in the collisions
        return self.cobraAssociations[:, collisions]

    @staticmethod
    def calculateCobraCenters(layout):
        """Calculates the cobras central positions for a given bench layout.

        Parameters
        ----------
        layout: str
            The bench layout to use: hex, line, rails or full.

        Returns
        -------
        object
            Complex numpy array with the cobras central positions.

        """
        if layout == "hex":
            # The centers should follow an hexagon pattern
            cobraCenters = np.empty(7, dtype="complex")
            cobraCenters[0] = 0.0
            cobraCenters[1:] = COBRAS_SEPARATION * np.exp(
                np.arange(6) * 1j * np.pi / 3)
        elif layout == "line":
            # Line of cobras
            cobraCenters = COBRAS_SEPARATION * np.arange(27) + 1j
        elif layout == "rails":
            # One PFI sector with 14 cobra modules
            cobraCenters = Bench.calculateFirstSectorCenters()
        elif layout == "full":
            # Full PFI bench (2394 cobras distributed in 3 rotated sectors)
            cobraCenters = Bench.calculatePFICenters()
        else:
            raise Exception(
                "%s is not a valid bench layout. Use full, hex, line, rails or "
                "calibration." % layout)

        return cobraCenters

    @staticmethod
    def calculateFirstSectorCenters():
        """Calculates the cobras central positions for the first PFI sector.

        Returns
        -------
        object
            Complex numpy array with the cobras central positions.

        """
        # Create the cobras centers array
        cobrasPerModule = MODULE_FIRST_LINE_LENGTH + MODULE_SECOND_LINE_LENGTH
        cobraCenters = np.empty(
            cobrasPerModule * MODULES_PER_SECTOR, dtype="complex")

        # Fill the first module
        firstModule = cobraCenters[:cobrasPerModule]
        firstModule[:MODULE_FIRST_LINE_LENGTH] = COBRAS_SEPARATION * np.arange(
            MODULE_FIRST_LINE_LENGTH)
        firstModule[MODULE_FIRST_LINE_LENGTH:] = COBRAS_SEPARATION * np.arange(
            MODULE_SECOND_LINE_LENGTH) + COBRAS_SEPARATION * np.exp(
                1j * np.pi / 3)
        firstModule += COBRAS_SEPARATION * np.exp(1j * 2 * np.pi / 3)

        # Order the first module centers by the x coordinate
        firstModule.sort()

        # Fill the rest of the modules
        modulesOffset = 2 * COBRAS_SEPARATION * np.exp(1j * 2 * np.pi / 3)

        for i in range(1, MODULES_PER_SECTOR):
            cobraCenters[i * cobrasPerModule:(i + 1) * cobrasPerModule] = (
                firstModule + i * modulesOffset)

        return cobraCenters

    @staticmethod
    def calculatePFICenters():
        """Calculates the cobras central positions for the full Prime Focus
        Instrument.

        Returns
        -------
        object
            Complex numpy array with the cobras central positions.

        """
        # Get the first sector cobra centers
        firstSector = Bench.calculateFirstSectorCenters()

        # Create the cobras centers array
        cobrasPerSector = len(firstSector)
        cobraCenters = np.empty(3 * cobrasPerSector, dtype="complex")

        # Add the first sector
        cobraCenters[:cobrasPerSector] = firstSector

        # Add the second sector rotating the first sector 120 degrees
        cobraCenters[cobrasPerSector:-cobrasPerSector] = firstSector * np.exp(
            1j * 2 * np.pi / 3)

        # Add the third sector rotating the first sector 240 degrees
        cobraCenters[-cobrasPerSector:] = firstSector * np.exp(
            1j * 4 * np.pi / 3)

        return cobraCenters

    @staticmethod
    def distancesBetweenLineSegments(startPoints1, endPoints1, startPoints2,
                                     endPoints2):
        """Calculates the minimum distances between line segments.

        Parameters
        ----------
        startPoints1: object
            A complex numpy array with the first line segments start
            coordinates.
        endPoints1: object
            A complex numpy array with the first line segments end coordinates.
        startPoints2: object
            A complex numpy array with the second line segments start
            coordinates.
        endPoints2: object
            A complex numpy array with the second line segments end
            coordinates.

        Returns
        -------
        object
            A numpy array with the minimum distance between the line segments.

        """
        # Calculate the minimum distances for each point to segment combination
        distances1 = Bench.distancesToLineSegments(
            startPoints1, startPoints2, endPoints2)
        distances2 = Bench.distancesToLineSegments(
            endPoints1, startPoints2, endPoints2)
        distances3 = Bench.distancesToLineSegments(
            startPoints2, startPoints1, endPoints1)
        distances4 = Bench.distancesToLineSegments(
            endPoints2, startPoints1, endPoints1)

        # Return the minimum distances
        return np.min((distances1, distances2, distances3, distances4), axis=0)

    @staticmethod
    def distancesToLineSegments(points, startPoints, endPoints):
        """Calculates the minimum distances between points and line segments.

        Parameters
        ----------
        points: object
            A complex numpy array with the point coordinates.
        startPoints: object
            A complex numpy array with the line segments start coordinates.
        endPoints: object
            A complex numpy array with the line segments end coordinates.

        Returns
        -------
        object
            A numpy array with the minimum distances between the points and the
            line segments.

        """
        # Translate the points and the line segment end points to the line
        # segment starting points
        translatedPoints = points - startPoints
        translatedEndPoints = endPoints - startPoints

        # Rotate the translated points to have the line segment on the x axis
        rotatedPoints = translatedPoints * np.exp(-1j * np.angle(
            translatedEndPoints))

        # Define 3 regions for the points: left of the origin, over the line
        # segments, and right of the line segments
        x = rotatedPoints.real
        lineLengths = np.abs(translatedEndPoints)
        (region1,) = np.where(x <= 0)
        (region2,) = np.where(np.logical_and(x > 0 , x < lineLengths))
        (region3,) = np.where(x >= lineLengths)

        # Calculate the minimum distances in each region
        distances = np.empty(len(points))
        distances[region1] = np.abs(rotatedPoints[region1])
        distances[region2] = np.abs(rotatedPoints[region2].imag)
        distances[region3] = np.abs(
            rotatedPoints[region3] - lineLengths[region3])

        return distances
