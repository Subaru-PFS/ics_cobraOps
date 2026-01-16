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
                 blackDotsMargin=1.0, fiducials=None):
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
        fiducials: object, optional
            A pandas DataFrame with the fiducials data. Default is None.

        Returns
        -------
        object
            The Bench instance.

        """
        # Create the cobra group instance
        self.cobras = CobraGroup(cobraCoach)

        # Create the black dot group instance
        self.blackDots = BlackDotGroup(
            blackDotsCalibrationProduct, blackDotsMargin)

        # Save the fiducials data
        self.fiducials = fiducials

        # Calculate the bench center
        self.center = np.mean(self.cobras.centers)

        # Calculate the bench radius
        self.radius = np.max(
            np.abs(self.cobras.centers - self.center) + self.cobras.rMax)

        # Calculate the cobra nearest neighbors associations array
        self.calculateCobraAssociations()

    def calculateCobraAssociations(self):
        """Calculates the cobras nearest neighbors associations array.

        Each column in the array contains a different cobra association.

        """
        # Calculate all the possible cobra associations
        indices = np.arange(self.cobras.nCobras)
        (cobraIndices, associatedCobraIndices) = np.where(
            indices[:, np.newaxis] < indices)

        # Calculate the cobra association distances
        cobraAssociationsDistances = np.abs(
            self.cobras.centers[cobraIndices] -
            self.cobras.centers[associatedCobraIndices])

        # Check which associations could result in possible collisions
        cobrasMaximumExtension = self.cobras.rMax + self.cobras.linkRadius
        possibleCollisions = cobraAssociationsDistances < 1.1 * (
            cobrasMaximumExtension[cobraIndices] +
            cobrasMaximumExtension[associatedCobraIndices])

        # Save the relevant cobra associations in a single array
        self.cobraAssociations = np.vstack((
            cobraIndices[possibleCollisions],
            associatedCobraIndices[possibleCollisions]))

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
        firstNeighborGroup = self.cobraAssociations[1][
            self.cobraAssociations[0] == cobraIndex]
        secondNeighborGroup = self.cobraAssociations[0][
            self.cobraAssociations[1] == cobraIndex]

        return np.concatenate((firstNeighborGroup, secondNeighborGroup))
