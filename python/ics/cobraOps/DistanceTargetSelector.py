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
    """Subclass of the TargetSelector class used to select optimal targets for
    a given PFI bench. The selection criteria is based on the target-to-cobra
    distance.

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

        # Select a single target for each cobra
        self.selectTargets()

    def selectTargets(self):
        """Selects a single target for each cobra based on their distance to
        the cobra center.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        # Create the array that will contain the assigned target indices
        self.assignedTargetIndices = np.full(
            self.bench.cobras.nCobras, NULL_TARGET_INDEX)

        # Assign targets to cobras, starting from the closest to the more far
        # away ones
        freeCobras = np.full(self.bench.cobras.nCobras, True)
        freeTargets = np.full(self.targets.nTargets, True)
        maxTargetsPerCobra = self.accessibleTargetIndices.shape[1]

        for i in range(maxTargetsPerCobra):
            # Get a list with the unique targets in this column
            indices = self.accessibleTargetIndices[:, i]
            uniqueIndices = np.unique(indices[freeCobras])

            # Remove from the list the NULL target index value if it's present
            uniqueIndices = uniqueIndices[uniqueIndices != NULL_TARGET_INDEX]

            # Select only those targets that are free
            uniqueIndices = uniqueIndices[freeTargets[uniqueIndices]]

            # Loop over the unique target indices
            for targetIndex in uniqueIndices:
                # Get the free cobras for which this target is the closest
                (cobraIndices,) = np.where(
                    np.logical_and(indices == targetIndex, freeCobras))

                # Check how many cobras we have
                if len(cobraIndices) == 1:
                    # Use this single cobra for this target
                    cobraToUse = cobraIndices[0]
                else:
                    # Select the cobras for which this is the only target
                    accessibleTargets = self.accessibleTargetIndices[
                        cobraIndices, i:]
                    targetIsAvailable = accessibleTargets != NULL_TARGET_INDEX
                    targetIsAvailable[targetIsAvailable] = freeTargets[
                        accessibleTargets[targetIsAvailable]]
                    nAvailableTargets = np.sum(targetIsAvailable, axis=1)
                    singleTargetCobras = cobraIndices[nAvailableTargets == 1]

                    # Decide depending on how many of these cobras we have
                    if len(singleTargetCobras) == 0:
                        # All cobras have multiple targets: select the closest
                        distances = self.accessibleTargetDistances[
                            cobraIndices, i]
                        cobraToUse = cobraIndices[distances.argmin()]
                    elif len(singleTargetCobras) == 1:
                        # Assign the target to the only single target cobra
                        cobraToUse = singleTargetCobras[0]
                    else:
                        # Assign the target to the closest single target cobra
                        distances = self.accessibleTargetDistances[
                            singleTargetCobras, i]
                        cobraToUse = singleTargetCobras[distances.argmin()]

                # Assign the target to the selected cobra
                self.assignedTargetIndices[cobraToUse] = targetIndex
                freeCobras[cobraToUse] = False
                freeTargets[targetIndex] = False
