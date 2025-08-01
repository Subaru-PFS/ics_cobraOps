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

    def calculateAccessibleTargets(self, maximumDistance=np.inf,
                                   safetyMargin=0):
        """Calculates the targets that each cobra can reach.

        The accessible targets are ordered randomly.

        This method should always be run before the selecTargets method.

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
        self._calculateAccessibleTargets(
            maximumDistance, safetyMargin, orderRandomly=True)

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
