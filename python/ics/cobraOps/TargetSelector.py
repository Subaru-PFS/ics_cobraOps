"""

TargetSelector abstract class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

from abc import ABC, abstractmethod

import numpy as np
from scipy.spatial import KDTree

from .TargetGroup import TargetGroup


class TargetSelector(ABC):
    """Abstract class used to select optimal targets for a given PFI bench.

    """

    def __init__(self, bench, targets):
        """Constructs a new TargetSelector instance.

        Parameters
        ----------
        bench: object
            The PFI bench instance.
        targets: object
            The TargetGroup instance.

        Returns
        -------
        object
            The TargetSelector instance.

        """
        # Save the bench and the targets
        self.bench = bench
        self.targets = targets

        # Construct a KD tree if the target density is large enough
        if self.targets.nTargets / self.bench.cobras.nCobras > 50:
            self.kdTree = self.constructKDTree()
        else:
            self.kdTree = None

        # Define some internal variables that will be used by the
        # calculateAccessibleTargets and selectTargets methods
        self.accessibleTargetIndices = None
        self.accessibleTargetDistances = None
        self.accessibleTargetElbows = None
        self.accessibleTargetPriorities = None
        self.assignedTargetIndices = None

    @abstractmethod
    def calculateAccessibleTargets(self, maximumDistance=np.inf,
                                   safetyMargin=0):
        """Calculates the targets that each cobra can reach.

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
        pass

    @abstractmethod
    def selectTargets(self):
        """Selects a single target for each cobra.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        pass

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
        # Obtain the accessible targets for each cobra
        self.calculateAccessibleTargets(maximumDistance, safetyMargin)

        # Select a single target for each cobra
        self.selectTargets()

    def constructKDTree(self, leafSize=None):
        """Constructs a K-dimensional tree using the target positions.

        The KD tree will then be used for the target-to-cobra distance
        calculations. The KD tree method starts to be efficient compared to
        brute force for target-to-cobra densities above 50.

        Parameters
        ----------
        leafSize: float, optional
            The KD tree leaf size. If None, an optimal leaf size value will be
            used. Default is None.

        Returns
        -------
        object
            K-dimensional tree constructed using the target positions.

        """
        # Calculate the KD tree leaf size if it is not provided
        if leafSize is None:
            leafSize = 3 * self.targets.nTargets // self.bench.cobras.nCobras
            leafSize = max(2, leafSize)

        # Return the KD tree
        return KDTree(
            np.column_stack(
                (self.targets.positions.real, self.targets.positions.imag)),
            leafsize=leafSize)

    def _getTargetsInsidePatrolArea(self, cobraIndex, maximumDistance=np.inf,
                                    safetyMargin=0):
        """Calculates the targets that fall inside a given cobra patrol area.

        This method doesn't consider if the cobra is broken or not.

        Parameters
        ----------
        cobraIndex: int
            The cobra index.
        maximumDistance: float, optional
            The maximum radial distance allowed between the targets and the
            cobra center. Default is no limit (the maximum radius that the
            cobra can reach).
        safetyMargin: float, optional
            Safety margin in mm added to Rmin and subtracted from Rmax to take
            into account possible effects that could change the effective cobra
            patrol area. Default is 0.

        Returns
        -------
        tuple
            A python tuple with the indices, positions and distances of the
            targets that fall inside the cobra patrol area. The arrays are
            ordered by their distance to the cobra center (closer targets
            appear first).

        """
        # Get the cobra center and the patrol area limits
        cobraCenter = self.bench.cobras.centers[cobraIndex]
        rMin = self.bench.cobras.rMin[cobraIndex] + safetyMargin
        rMax = min(
            self.bench.cobras.rMax[cobraIndex] - safetyMargin, maximumDistance)

        # If available, use the KD tree for the distance calculations
        if self.kdTree is not None:
            # Get all the targets that the cobra can reach. Remember to
            # remove any possible NULL targets that might exist
            distances, indices = self.kdTree.query(
                [cobraCenter.real, cobraCenter.imag], k=None,
                distance_upper_bound=rMax)
            indices = np.array(indices, dtype=np.int)
            distances = np.array(distances)
            validTargets = np.logical_and(
                self.targets.notNull[indices], distances > rMin)
            indices = indices[validTargets]
            distances = distances[validTargets]
            positions = self.targets.positions[indices]
        else:
            # Get all the targets that the cobra can reach. Remember to
            # remove any possible NULL targets that might exist
            xDistances = np.abs(cobraCenter.real - self.targets.positions.real)
            (indices,) = np.where(
                np.logical_and(self.targets.notNull, xDistances < rMax))
            positions = self.targets.positions[indices]
            yDistances = np.abs(cobraCenter.imag - positions.imag)
            validTargets = yDistances < rMax
            indices = indices[validTargets]
            positions = positions[validTargets]
            distances = np.abs(cobraCenter - positions)
            validTargets = np.logical_and(distances > rMin, distances < rMax)
            indices = indices[validTargets]
            positions = positions[validTargets]
            distances = distances[validTargets]

            # Sort the targets by their distance to the cobra center
            sortedIndices = distances.argsort()
            indices = indices[sortedIndices]
            positions = positions[sortedIndices]
            distances = distances[sortedIndices]

        return indices, positions, distances

    def _calculateAccessibleTargets(self, maximumDistance, safetyMargin,
                                    orderRandomly=False, orderByPriority=False):
        """Calculates the targets that each cobra can reach.

        By default accessible targets are ordered by their distance to the cobra
        center, unless orderRandomly or orderByPriority are set to True.

        Parameters
        ----------
        maximumDistance: float
            The maximum radial distance allowed between the targets and the
            cobra centers.
        safetyMargin: float
            Safety margin in mm added to Rmin and subtracted from Rmax to take
            into account possible effects that could change the effective cobra
            patrol area.
        orderRandomly: bool, optional
            If True, accessible targets will be ordered randomly. Default is
            False, which means that the targets will be ordered by their
            distance to the cobra center.
        orderByPriority: bool, optional
            If True, accessible targets will be ordered by their priority.
            Default is False, which means that the targets will be ordered
            by their distance to the cobra center.

        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        blackDotsPositions = self.bench.blackDots.centers
        blackDotsRadius = self.bench.blackDots.radius

        # Obtain the cobra-target associations
        associations = []
        maxTargetsPerCobra = 0

        for i in range(nCobras):
            # Get the targets that fall inside the cobra patrol area
            indices, positions, distances = self._getTargetsInsidePatrolArea(
                i, maximumDistance, safetyMargin)

            # Invalidate all the targets if the cobra has a problem
            if not self.bench.cobras.isGood[i]:
                indices = indices[[]]
                positions = positions[[]]
                distances = distances[[]]

            # Invalidate targets falling inside the black dots
            blackDotdistances = np.abs(
                positions[:, np.newaxis] - blackDotsPositions)
            validTargets = np.all(blackDotdistances > blackDotsRadius, axis=1)
            indices = indices[validTargets]
            positions = positions[validTargets]
            distances = distances[validTargets]

            # Save the cobra-target association information
            associations.append((i, indices, positions, distances))
            maxTargetsPerCobra = max(maxTargetsPerCobra, len(indices))

        # Create the accessible target arrays
        arrayShape = (nCobras, maxTargetsPerCobra)
        self.accessibleTargetIndices = np.full(
            arrayShape, TargetGroup.NULL_TARGET_INDEX)
        self.accessibleTargetDistances = np.zeros(arrayShape)
        self.accessibleTargetElbows = np.zeros(arrayShape, dtype=complex)
        self.accessibleTargetPriorities = np.zeros(arrayShape)

        # Fill the arrays with the cobra-target association information
        for i, indices, positions, distances in associations:
            # Get the total number of accessible targets for this cobra
            nTargets = len(indices)

            # Calculate the elbow positions at the target positions
            elbows = self.bench.cobras.calculateCobraElbowPositions(
                i, positions)

            # Get the target priorities
            priorities = self.targets.priorities[indices]

            # Check if we need to order the accessible targets randomly
            if orderRandomly or orderByPriority:
                # Randomize the targets order to remove the distance order
                randomOrder = np.random.permutation(nTargets)
                indices = indices[randomOrder]
                distances = distances[randomOrder]
                elbows = elbows[randomOrder]
                priorities = priorities[randomOrder]

            # Check if we need to order the accessible targets by their priority
            if orderByPriority:
                # Order the targets by their priority
                priorityOrder = np.argsort(priorities)[::-1]
                indices = indices[priorityOrder]
                distances = distances[priorityOrder]
                elbows = elbows[priorityOrder]
                priorities = priorities[priorityOrder]

            # Fill the arrays
            self.accessibleTargetIndices[i, :nTargets] = indices
            self.accessibleTargetDistances[i, :nTargets] = distances
            self.accessibleTargetElbows[i, :nTargets] = elbows
            self.accessibleTargetPriorities[i, :nTargets] = priorities

    def getAccessibleTargetsInformation(self, cobraIndex):
        """Returns the indices, distances to the cobra center and cobra elbow
        positions for those targets that can be accessed by the given cobra.

        This method should not be run before the calculateAccessibleTargets
        method.

        Returns
        -------
        tuple
            A python tuple with the indices, distances to the cobra centers and
            cobra elbow positions for those targets that can be accessed by the
            cobra. The arrays are ordered by the target distance to the cobra
            center (closer targets appear first).

        """
        # Get the cobra accessible targets information
        indices = self.accessibleTargetIndices[cobraIndex]
        distances = self.accessibleTargetDistances[cobraIndex]
        elbows = self.accessibleTargetElbows[cobraIndex]

        # Return only the valid targets
        notNull = indices != TargetGroup.NULL_TARGET_INDEX

        return indices[notNull], distances[notNull], elbows[notNull]

    def getSelectedTargets(self):
        """Returns a new target group with the selected target for each cobra.

        This method should not be run before the selecTargets method.

        Returns
        -------
        object
            A new target group instance with the selected target for each cobra.

        """
        return self.targets.select(self.assignedTargetIndices)
