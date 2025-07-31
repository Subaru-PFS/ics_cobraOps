"""

PriorityTargetSelector class.

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


class PriorityTargetSelector(TargetSelector):
    """Subclass of the TargetSelector class used to select optimal targets for
    a given PFI bench. The selection criteria is based on the target priority.

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

        # Order the accessible targets by their priority
        self.orderAccessibleTargetsByPriority()

        # Select a single target for each cobra
        self.selectTargets()

    def orderAccessibleTargetsByPriority(self):
        """Orders the accessible targets arrays by decreasing target priority.

        """
        # Create the array that will contain the accessible target priorities
        arrayShape = self.accessibleTargetIndices.shape
        self.accessibleTargetPriorities = np.zeros(arrayShape)

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

            # Order the targets by their priority
            priorities = self.targets.priorities[indices[:nTargets]]
            priorityOrder = np.argsort(priorities)[::-1]
            indices[:nTargets] = indices[priorityOrder]
            distances[:nTargets] = distances[priorityOrder]
            elbows[:nTargets] = elbows[priorityOrder]

            # Save the targets priorities
            self.accessibleTargetPriorities[i, :nTargets] = priorities[
                priorityOrder]

    def selectTargets(self):
        """Selects a single target for each cobra based on their priorities.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        # Create the array that will contain the assigned target indices
        self.assignedTargetIndices = np.full(
            self.bench.cobras.nCobras, TargetGroup.NULL_TARGET_INDEX)

        # Assign targets to cobras, starting from the highest priorities
        freeCobras = np.full(self.bench.cobras.nCobras, True)
        freeTargets = np.full(self.targets.nTargets, True)
        maxTargetsPerCobra = self.accessibleTargetIndices.shape[1]

        for i in range(maxTargetsPerCobra):
            # Get a list with the unique targets in this column
            indices = self.accessibleTargetIndices[:, i]
            uniqueIndices = np.unique(indices[freeCobras])

            # Remove from the list the NULL target index value if it's present
            uniqueIndices = uniqueIndices[
                uniqueIndices != TargetGroup.NULL_TARGET_INDEX]

            # Select only those targets that are free
            uniqueIndices = uniqueIndices[freeTargets[uniqueIndices]]

            # Loop over the unique target indices
            for targetIndex in uniqueIndices:
                # Get the free cobras for which this target has the highest
                # priority
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
                    targetIsAvailable = (
                        accessibleTargets != TargetGroup.NULL_TARGET_INDEX)
                    targetIsAvailable[targetIsAvailable] = freeTargets[
                        accessibleTargets[targetIsAvailable]]
                    nAvailableTargets = np.sum(targetIsAvailable, axis=1)
                    singleTargetCobras = cobraIndices[nAvailableTargets == 1]

                    # Decide depending on how many single target cobras we have
                    if len(singleTargetCobras) == 0:
                        # All cobras have multiple targets: select the cobra
                        # with the lowest priority sum
                        prioritySum = np.sum(self.accessibleTargetPriorities[
                            cobraIndices, i:] * targetIsAvailable, axis=1)
                        cobraToUse = cobraIndices[prioritySum.argmin()]
                    elif len(singleTargetCobras) == 1:
                        # Assign the target to the only single target cobra
                        cobraToUse = singleTargetCobras[0]
                    else:
                        # Assign the target to the single target cobra with the
                        # lowest priority sum
                        prioritySum = np.sum(self.accessibleTargetPriorities[
                            cobraIndices, i:] * targetIsAvailable, axis=1)
                        prioritySum = prioritySum[nAvailableTargets == 1]
                        cobraToUse = singleTargetCobras[prioritySum.argmin()]

                # Assign the target to the selected cobra
                self.assignedTargetIndices[cobraToUse] = targetIndex
                freeCobras[cobraToUse] = False
                freeTargets[targetIndex] = False
