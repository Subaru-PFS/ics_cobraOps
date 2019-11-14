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

from .cobraConstants import NULL_TARGET_INDEX
from .TargetSelector import TargetSelector


class RandomTargetSelector(TargetSelector):
    """

    Subclass of the TargetSelector class used to select optimal targets for a
    given PFI bench. The selection criteria is completely random.

    """

    def run(self, maximumDistance=np.Inf, solveCollisions=True):
        """Runs the whole target selection process assigning a single target to
        each cobra in the bench.

        Parameters
        ----------
        maximumDistance: float, optional
            The maximum radial distance allowed between the targets and the
            cobra centers. Default is no limit (the maximum radius that the
            cobra can reach).
        solveCollisions: bool, optional
            If True, the selector will try to solve cobra collisions assigning
            them alternative targets. Default is True.

        """
        # Construct a KD tree if the target density is large enough
        if self.targets.nTargets / self.bench.cobras.nCobras > 50:
            self.constructKDTree()

        # Obtain the accessible targets for each cobra ordered by distance
        self.calculateAccessibleTargets(maximumDistance)

        # Randomize the accessible targets order
        self.randomizeAccessibleTargetsOrder()

        # Select a single target for each cobra
        self.selectTargets()

        # Try to solve end point collisions
        if solveCollisions:
            self.solveEndPointCollisions()

    def randomizeAccessibleTargetsOrder(self):
        """Randomizes the accessible targets order.

        """
        # Loop over the cobras
        for i in range(self.bench.cobras.nCobras):
            # Get the accessible target indices, distances and elbow positions
            indices = self.accessibleTargetIndices[i]
            distances = self.accessibleTargetDistances[i]
            elbows = self.accessibleTargetElbows[i]

            # Randomize the targets order
            nTargets = np.sum(indices != NULL_TARGET_INDEX)
            newOrder = np.random.permutation(nTargets)
            indices[:nTargets] = indices[newOrder]
            distances[:nTargets] = distances[newOrder]
            elbows[:nTargets] = elbows[newOrder]

    def selectTargets(self):
        """Selects a single random target for each cobra.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        nTargets = self.targets.nTargets

        # Assign random targets to cobras
        self.assignedTargetIndices = np.full(nCobras, NULL_TARGET_INDEX)
        freeTargets = np.full(nTargets, True)

        for i in range(nCobras):
            # Get the accessible targets to this cobra
            targetIndices = self.accessibleTargetIndices[i]

            # Remove from the list the NULL_TARGET_INDEX value if it's present
            targetIndices = targetIndices[targetIndices != NULL_TARGET_INDEX]

            # Select free targets only
            targetIndices = targetIndices[freeTargets[targetIndices]]

            # Check that there is at least one target available
            if len(targetIndices) > 0:
                # Select the first target
                targetIndex = targetIndices[0]

                # Assign the target to the selected cobra
                self.assignedTargetIndices[i] = targetIndex
                freeTargets[targetIndex] = False
