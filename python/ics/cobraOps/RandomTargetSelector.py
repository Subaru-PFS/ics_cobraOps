"""

RandomTargetSelector class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from .TargetGroup import TargetGroup
from .TargetSelector import TargetSelector


class RandomTargetSelector(TargetSelector):
    """Subclass of the TargetSelector class used to select optimal targets for
    a given PFI bench. The selection criteria is completely random.

    """

    def run(self, maximumDistance=np.inf, safetyMargin=0):
        """Runs the whole target selection process assigning a single target to
        each cobra in the bench.

        Parameters
        ----------
        maximumDistance: float, optional
            The maximum radial distance allowed between the targets and the
            cobra centers. Default is no limit (the maximum radius that the
            cobra can reach).
        safetyMargin: float, optional
            Safety margin in mm added to Rmin and subtracted from Rmax to take
            into account possible effects that could change the effective cobra
            patrol area. Default is 0.

        """
        # Construct a KD tree if the target density is large enough
        if self.targets.nTargets / self.bench.cobras.nCobras > 50:
            self.constructKDTree()

        # Obtain the accessible targets for each cobra ordered by distance
        self.calculateAccessibleTargets(maximumDistance, safetyMargin)

        # Order the accessible targets randomly
        self.orderAccessibleTargetsRandomly()

        # Select a single target for each cobra
        self.selectTargets()

    def orderAccessibleTargetsRandomly(self):
        """Orders the accessible targets arrays randomly.

        """
        # Loop over the cobras
        for i in range(self.bench.cobras.nCobras):
            # Get the accessible target indices, distances and elbow positions
            indices = self.accessibleTargetIndices[i]
            distances = self.accessibleTargetDistances[i]
            elbows = self.accessibleTargetElbows[i]
            nTargets = np.sum(indices != TargetGroup.NULL_TARGET_INDEX)

            # Randomize the targets order to remove the distance order
            randomOrder = np.random.permutation(nTargets)
            indices[:nTargets] = indices[randomOrder]
            distances[:nTargets] = distances[randomOrder]
            elbows[:nTargets] = elbows[randomOrder]

    def selectTargets(self):
        """Selects a single random target for each cobra.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        # Create the array that will contain the assigned target indices
        self.assignedTargetIndices = np.full(
            self.bench.cobras.nCobras, TargetGroup.NULL_TARGET_INDEX)

        # Calculate the number of accessible targets per cobra
        nTargetsPerCobra = np.sum(
            self.accessibleTargetIndices != TargetGroup.NULL_TARGET_INDEX,
            axis=1)

        # Assign random targets to cobras, starting with those cobras with fewer
        # accessible targets
        cobraIndices = np.argsort(nTargetsPerCobra)
        freeTargets = np.full(self.targets.nTargets, True)

        for i in cobraIndices:
            # Get the indices of the accessible targets to this cobra
            indices = self.accessibleTargetIndices[i]
            indices = indices[indices != TargetGroup.NULL_TARGET_INDEX]

            # Select only those targets that are free
            indices = indices[freeTargets[indices]]

            # Check that there is at least one target available
            if len(indices) > 0:
                # Select the first target
                targetIndex = indices[0]

                # Assign the target to the cobra
                self.assignedTargetIndices[i] = targetIndex
                freeTargets[targetIndex] = False
