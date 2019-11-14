"""

DistanceTargetSelector class.

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


class DistanceTargetSelector(TargetSelector):
    """

    Subclass of the TargetSelector class used to select optimal targets for a
    given PFI bench. The selection criteria is based on the target to cobra
    distance.

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

        # Obtain the accessible targets for each cobra
        self.calculateAccessibleTargets(maximumDistance)

        # Select a single target for each cobra
        self.selectTargets()

        # Try to solve end point collisions
        if solveCollisions:
            self.solveEndPointCollisions()

    def selectTargets(self):
        """Selects a single target for each cobra based on their distance to
        the cobra center.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        nTargets = self.targets.nTargets
        maxTargetsPerCobra = self.accessibleTargetIndices.shape[1]

        # Assign targets to cobras looping from the closest to the more far
        # away ones
        self.assignedTargetIndices = np.full(nCobras, NULL_TARGET_INDEX)
        freeCobras = np.full(nCobras, True)
        freeTargets = np.full(nTargets, True)

        for i in range(maxTargetsPerCobra):
            # Get a list with the unique targets in the given column
            columnTargetIndices = self.accessibleTargetIndices[:, i]
            uniqueTargetIndices = np.unique(columnTargetIndices[freeCobras])

            # Remove from the list the NULL_TARGET_INDEX value if it's present
            uniqueTargetIndices = uniqueTargetIndices[uniqueTargetIndices != NULL_TARGET_INDEX]

            # Select free targets only
            uniqueTargetIndices = uniqueTargetIndices[freeTargets[uniqueTargetIndices]]

            # Loop over the unique target indices
            for targetIndex in uniqueTargetIndices:
                # Get the free cobras for which this target is the closest in
                # the current column
                (associatedCobras,) = np.where(np.logical_and(columnTargetIndices == targetIndex, freeCobras))

                # Check how many associated cobras we have
                if len(associatedCobras) == 1:
                    # Use this single cobra for this target
                    cobraToUse = associatedCobras[0]
                else:
                    # Select the cobras for which this is the only target
                    accessibleTargets = self.accessibleTargetIndices[associatedCobras, i:]
                    targetIsAvailable = np.logical_and(accessibleTargets != NULL_TARGET_INDEX, freeTargets[accessibleTargets])
                    nAvailableTargets = np.sum(targetIsAvailable, axis=1)
                    singleTargetCobras = associatedCobras[nAvailableTargets == 1]

                    # Decide depending on how many of these cobras we have
                    if len(singleTargetCobras) == 0:
                        # All cobras have multiple targets: select the closest
                        # cobra
                        distances = self.accessibleTargetDistances[associatedCobras, i]
                        cobraToUse = associatedCobras[distances.argmin()]
                    elif len(singleTargetCobras) == 1:
                        # Assign the target to the cobra that can only reach this
                        # target
                        cobraToUse = singleTargetCobras[0]
                    else:
                        # Assign the target to the closest single target cobra
                        distances = self.accessibleTargetDistances[singleTargetCobras, i]
                        cobraToUse = singleTargetCobras[distances.argmin()]

                # Assign the target to the selected cobra
                self.assignedTargetIndices[cobraToUse] = targetIndex
                freeCobras[cobraToUse] = False
                freeTargets[targetIndex] = False
