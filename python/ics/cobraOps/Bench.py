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
from .CobraGroup import CobraGroup


class Bench:
    """Class describing the properties of a PFI bench.

    """

    def __init__(self, cobraCoach, blackDotsCalibrationProduct,
                 blackDotsMargin=1.0):
        """Constructs a new Bench instance.

        Parameters
        ----------
        cobraCoach: object
            The cobra coach instance.
        blackDotsCalibrationProduct: object
            The black dots calibration product with the black dots properties
            to use.
        blackDotsMargin: real, optional
            The margin factor in radius of the black dots to avoid in fiber
            allocation. Default is 1.0.

        Returns
        -------
        object
            The Bench instance.

        """
        # Create the cobra group instance
        self.cobras = CobraGroup(cobraCoach)

        # Calculate the bench center
        self.center = np.mean(self.cobras.centers[self.cobras.isGood])

        # Calculate the bench radius
        self.radius = np.max(
            np.abs(self.cobras.centers[self.cobras.isGood] - self.center) +
            self.cobras.rMax[self.cobras.isGood])

        # Calculate the cobra nearest neighbors associations array
        self.calculateCobraAssociations()

        # Create the black dot group instance
        self.blackDots = BlackDotGroup(
            blackDotsCalibrationProduct, blackDotsMargin)

    def calculateCobraAssociations(self):
        """Calculates the cobras nearest neighbors associations array.

        Each column in the array contains a different cobra association.

        """
        # Calculate the cobras distance matrix
        distanceMatrix = np.abs(
            self.cobras.centers[:, np.newaxis] - self.cobras.centers)

        # Mask the diagonal because it contains distances to the same cobras
        distanceMatrix[np.arange(self.cobras.nCobras),
                       np.arange(self.cobras.nCobras)] = np.inf

        # Calculate the median minimum distance between cobras
        medianMinDistance = np.median(
            np.min(distanceMatrix, axis=1)[self.cobras.isGood])

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
