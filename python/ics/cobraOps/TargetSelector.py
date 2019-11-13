"""

TargetSelector abstract class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np
from abc import ABC, abstractmethod
from scipy.spatial import KDTree

from .cobraConstants import NULL_TARGET_INDEX


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

        # Define some internal variables that will be used by the
        # computeAccessibleTargets and selectTargets methods
        self.kdTree = None
        self.accessibleTargetIndices = None
        self.accessibleTargetDistances = None
        self.accessibleTargetElbows = None
        self.assignedTargetIndices = None

    @abstractmethod
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
        pass

    @abstractmethod
    def selectTargets(self):
        """Selects a single target for each cobra.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        pass

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

        """
        # Calculate the KD tree leaf size if it is not provided
        if leafSize is None:
            leafSize = 3 * self.targets.nTargets // self.bench.cobras.nCobras
            leafSize = max(2, leafSize)

        # Construct the KD tree
        self.kdTree = KDTree(np.column_stack((self.targets.positions.real,
                                              self.targets.positions.imag)),
                                              leafsize=leafSize)

    def getTargetsInsidePatrolArea(self, cobraIndex, maximumDistance=np.Inf):
        """Calculates the targets that fall inside a given cobra patrol area.

        Parameters
        ----------
        cobraIndex: int
            The cobra index.
        maximumDistance: float, optional
            The maximum radial distance allowed between the targets and the
            cobra center. Default is no limit (the maximum radius that the
            cobra can reach).

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
        rMin = self.bench.cobras.rMin[cobraIndex]
        rMax = min(self.bench.cobras.rMax[cobraIndex], maximumDistance)

        # If available, use the KD tree for the distance calculations
        if self.kdTree is not None:
            # Get all the targets that the cobra can reach. Remember to
            # invalidate any possible NULL targets that might exist
            distances, indices = self.kdTree.query(
                [cobraCenter.real, cobraCenter.imag], k=None,
                distance_upper_bound=rMax)
            indices = np.array(indices, dtype=np.int)
            distances = np.array(distances)
            validTargets = np.logical_and(self.targets.notNull[indices],
                                          distances > rMin)
            indices = indices[validTargets]
            distances = distances[validTargets]
            positions = self.targets.positions[indices]
        else:
            # Get all the targets that the cobra can reach. Remember to
            # invalidate any possible NULL targets that might exist
            xDistances = np.abs(cobraCenter.real - self.targets.positions.real)
            (indices,) = np.where(np.logical_and(self.targets.notNull,
                                                 xDistances < rMax))
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

    def calculateAccessibleTargets(self, maximumDistance=np.Inf):
        """Calculates the targets that each cobra can reach.

        The results are saved in the accessibleTargetIndices,
        accesssibleTargetDistances and accessibleTargetElbows internal arrays.

        This method should always be run before the selecTargets method.

        Parameters
        ----------
        maximumDistance: float, optional
            The maximum radial distance allowed between the targets and the
            cobra centers. Default is no limit (the maximum radius that the
            cobra can reach).

        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        cobraCenters = self.bench.cobras.centers
        cobraBlackDotPosition = self.bench.cobras.blackDotPosition
        cobraBlackDotRadius = self.bench.cobras.blackDotRadius

        # Obtain the cobra-target associations
        associations = []
        maxTargetsPerCobra = 0

        for i in range(nCobras):
            # Get the targets that fall inside the cobra patrol area
            indices, positions, distances = self.getTargetsInsidePatrolArea(
                i, maximumDistance)

            # Invalidate targets falling in the cobra black dots
            blackDotdistances = np.abs(
                cobraCenters[i] + cobraBlackDotPosition[i] - positions)
            validTargets = blackDotdistances > cobraBlackDotRadius[i]
            indices = indices[validTargets]
            positions = positions[validTargets]
            distances = distances[validTargets]

            # Save the cobra-target association information
            associations.append((i, indices, positions, distances))
            maxTargetsPerCobra = max(maxTargetsPerCobra, len(indices))

        # Create the accessible target arrays
        arrayShape = (nCobras, maxTargetsPerCobra)
        self.accessibleTargetIndices = np.full(arrayShape, NULL_TARGET_INDEX)
        self.accessibleTargetDistances = np.zeros(arrayShape)
        self.accessibleTargetElbows = np.zeros(arrayShape, dtype="complex")

        # Fill the arrays with the cobra-target association information
        for i, indices, positions, distances in associations:
            # Calculate the elbow positions at the target positions
            elbows = self.bench.cobras.calculateCobraElbowPositions(
                i, positions)

            # Fill the arrays
            nTargets = len(indices)
            self.accessibleTargetIndices[i, :nTargets] = indices
            self.accessibleTargetDistances[i, :nTargets] = distances
            self.accessibleTargetElbows[i, :nTargets] = elbows

    def solveEndPointCollisions(self):
        """Detects and solves cobra collisions assigning them alternative
        targets.

        This method should always be run after the selectTargets method.

        """
        # Extract some useful information
        nTargets = self.targets.nTargets
        targetPositions = self.targets.positions
        cobraCenters = self.bench.cobras.centers

        # Set the fiber positions to their associated target positions, leaving
        # unused cobras at their home positions
        fiberPositions = self.bench.cobras.home0.copy()
        usedCobras = self.assignedTargetIndices != NULL_TARGET_INDEX
        fiberPositions[usedCobras] = targetPositions[
            self.assignedTargetIndices[usedCobras]]

        # Get the indices of the cobra associations where we have a collision
        problematicCobraAssociations = self.bench.getProblematicCobraAssociations(
            fiberPositions)
        problematicCobras = problematicCobraAssociations[0]
        nearbyProblematicCobras = problematicCobraAssociations[1]

        # Try to solve the collisions one by one
        freeTargets = np.full(nTargets, True)
        freeTargets[self.assignedTargetIndices[usedCobras]] = False

        for c, nc in zip(problematicCobras, nearbyProblematicCobras):
            # Check if one of the colliding cobras is unused
            if self.assignedTargetIndices[c] == NULL_TARGET_INDEX or self.assignedTargetIndices[nc] == NULL_TARGET_INDEX:
                # The unused cobra is the cobra that we are going to move
                cobraToMove = c if (self.assignedTargetIndices[c] == NULL_TARGET_INDEX) else nc

                # Move the cobra until we find a position with zero collisions
                cobraCenter = cobraCenters[cobraToMove]
                initialPosition = fiberPositions[cobraToMove]
                bestPosition = initialPosition

                for ang in np.linspace(0, 2 * np.pi, 7)[1:-1]:
                    # Rotate the cobra around its center
                    fiberPositions[cobraToMove] = (initialPosition - cobraCenter) * np.exp(1j * ang) + cobraCenter

                    # Exit the loop if we found a fiber position with zero
                    # collisions
                    if self.bench.getCollisionsForCobra(cobraToMove, fiberPositions) == 0:
                        bestPosition = fiberPositions[cobraToMove]
                        break

                # Use the best fiber position
                fiberPositions[cobraToMove] = bestPosition
            else:
                # Calculate the initial number of collisions associated with
                # the two cobras
                collisions = self.bench.getCollisionsForCobra(c, fiberPositions)
                collisions += self.bench.getCollisionsForCobra(nc, fiberPositions)

                # Free the current targets
                initialTarget1 = self.assignedTargetIndices[c]
                initialTarget2 = self.assignedTargetIndices[nc]
                freeTargets[initialTarget1] = True
                freeTargets[initialTarget2] = True

                # Get the targets that can be reached by each cobra
                targets1 = self.accessibleTargetIndices[c, self.accessibleTargetIndices[c] != NULL_TARGET_INDEX]
                targets2 = self.accessibleTargetIndices[nc, self.accessibleTargetIndices[nc] != NULL_TARGET_INDEX]

                # Select only the free targets
                targets1 = targets1[freeTargets[targets1]]
                targets2 = targets2[freeTargets[targets2]]

                # Shuffle the arrays to remove the ordering by distance
                np.random.shuffle(targets1)
                np.random.shuffle(targets2)

                # Create two arrays reflecting all the possible target
                # combinations
                targetsCombination1 = np.repeat(targets1, len(targets2))
                targetsCombination2 = np.tile(targets2, len(targets1))

                # Exclude the current target combination and combinations that
                # use the same target for the two cobras
                validCombinations = np.logical_or(targetsCombination1 != initialTarget1, targetsCombination2 != initialTarget2)
                validCombinations = np.logical_and(validCombinations, targetsCombination1 != targetsCombination2)
                targetsCombination1 = targetsCombination1[validCombinations]
                targetsCombination2 = targetsCombination2[validCombinations]

                # Loop over all the possible combinations until we find the
                # minimum number of collisions
                bestTarget1 = initialTarget1
                bestTarget2 = initialTarget2

                for newTarget1, newTarget2 in zip(targetsCombination1, targetsCombination2):
                    # Assign the new fiber positions
                    fiberPositions[c] = targetPositions[newTarget1]
                    fiberPositions[nc] = targetPositions[newTarget2]

                    # Calculate the number of collisions at the current
                    # positions
                    currentCollisions = self.bench.getCollisionsForCobra(c, fiberPositions)
                    currentCollisions += self.bench.getCollisionsForCobra(nc, fiberPositions)

                    # Check if the number of collisions decreased
                    # significantly. A decrease of one means that we solved the
                    # current collision, but we created a new collision with
                    # another nearby cobra.
                    if currentCollisions <= collisions - 2:
                        # Save the information from this target combination
                        bestTarget1 = newTarget1
                        bestTarget2 = newTarget2
                        collisions = currentCollisions

                    # Exit the loop if the number of collisions is already zero
                    if collisions == 0:
                        break

                # Use the target combination where we had less collisions
                self.assignedTargetIndices[c] = bestTarget1
                self.assignedTargetIndices[nc] = bestTarget2
                fiberPositions[c] = targetPositions[bestTarget1]
                fiberPositions[nc] = targetPositions[bestTarget2]
                freeTargets[bestTarget1] = False
                freeTargets[bestTarget2] = False

    def getSelectedTargets(self):
        """Returns a new target group with a selected target for each cobra.

        This method should not be run before the selecTargets method.

        Returns
        -------
        object
            A new target group instance with a selected target for each cobra.

        """
        return self.targets.select(self.assignedTargetIndices)
