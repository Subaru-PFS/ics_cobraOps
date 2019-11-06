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
        # Obtain the accessible targets for each cobra
        self.calculateAccessibleTargets(maximumDistance)

        # Select a single target for each cobra
        self.selectTargets()

        # Try to solve end point collisions
        if solveCollisions:
            self.solveEndPointCollisions()


    def selectTargets(self):
        """Selects a single random target for each cobra.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        nTargets = self.targets.nTargets

        # Assign random targets to cobras
        self.assignedTargetIndices = np.full(nCobras, NULL_TARGET_INDEX, dtype="int")
        freeTargets = np.full(nTargets, True)

        for c in range(nCobras):
            # Get the accessible targets to this cobra
            targetIndices = self.accessibleTargetIndices[c]

            # Remove from the list the NULL_TARGET_INDEX value if it's present
            targetIndices = targetIndices[targetIndices != NULL_TARGET_INDEX]

            # Select free targets only
            targetIndices = targetIndices[freeTargets[targetIndices]]

            # Check that there is at least one target available
            if len(targetIndices) > 0:
                # Select one target randomly
                targetIndex = np.random.choice(targetIndices)

                # Assign the target to the selected cobra
                self.assignedTargetIndices[c] = targetIndex
                freeTargets[targetIndex] = False
