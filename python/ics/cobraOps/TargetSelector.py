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
            If True, the selector will try to solve cobra end-point collisions
            assigning them alternative targets. Default is True.

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
        blackDotsPositions = self.bench.blackDots.centers
        blackDotsRadius = self.bench.blackDots.radius

        # Obtain the cobra-target associations
        associations = []
        maxTargetsPerCobra = 0

        for i in range(nCobras):
            # Get the targets that fall inside the cobra patrol area
            indices, positions, distances = self.getTargetsInsidePatrolArea(
                i, maximumDistance)

            # Invalidate all the targets if the cobra has a problem
            if self.bench.cobras.hasProblem[i]:
                indices = indices[[]]
                positions = positions[[]]
                distances = distances[[]]

            # Invalidate targets falling in the nearby black dots
            blackDotsIndices = np.append([i], self.bench.getCobraNeighbors(i))
            blackDotdistances = np.abs(
                positions[:, np.newaxis] - blackDotsPositions[blackDotsIndices])
            validTargets = np.all(
                blackDotdistances > blackDotsRadius[blackDotsIndices], axis=1)
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
        """Detects and solves cobra end-point collisions assigning them
        alternative targets.

        This method should always be run after the selectTargets method.

        """
        # Get the indices of the targets that are currently assigned to cobras
        indices = self.assignedTargetIndices

        # Check which cobras are used (i.e. have a target assigned)
        usedCobras = indices != NULL_TARGET_INDEX

        # Check which targets are currently free (i.e. not assigned to a cobra)
        freeTargets = np.full(self.targets.nTargets, True)
        freeTargets[indices[usedCobras]] = False

        # Set the cobra fiber positions to the assigned target positions,
        # leaving unused cobras at their home positions
        positions = self.bench.cobras.home0.copy()
        positions[usedCobras] = self.targets.positions[indices[usedCobras]]

        # Get the cobra associations where we have an end-point collision
        problematicAssociations = self.bench.getProblematicCobraAssociations(
            positions).T

        # Try to solve the cobra collisions one by one
        for c, nc in problematicAssociations:
            # Check if one of the colliding cobras is not used
            if not usedCobras[c] or not usedCobras[nc]:
                # The unused cobra is the cobra that we are going to move
                cobraToMove = c if not usedCobras[c] else nc

                # Calculate the initial number of collisions for that cobra
                initialCollisions = self.bench.getCollisionsForCobra(
                    cobraToMove, positions)

                # Move to the next association if the number of collisions is
                # already zero (it could have been solved in a previous step)
                if initialCollisions == 0:
                    continue

                # Move the cobra until we find the position with the minimum
                # number of collisions
                cobraCenter = self.bench.cobras.centers[cobraToMove]
                initialPosition = positions[cobraToMove]
                bestPosition = initialPosition
                bestCollisions = initialCollisions

                for ang in np.linspace(0, 2 * np.pi, 7)[1:-1]:
                    # Rotate the cobra around its center
                    positions[cobraToMove] = cobraCenter + (
                        initialPosition - cobraCenter) * np.exp(1j * ang)

                    # Calculate the number of collisions at the rotated position
                    collisions = self.bench.getCollisionsForCobra(
                        cobraToMove, positions)

                    # Check if the number of collisions decreased
                    if collisions < bestCollisions:
                        # Save the information from this cobra position
                        bestPosition = positions[cobraToMove]
                        bestCollisions = collisions

                        # Exit the loop if the number of collisions is zero
                        if collisions == 0:
                            break

                # Use the best fiber position
                positions[cobraToMove] = bestPosition
            else:
                # Calculate the initial number of collisions associated with
                # the two cobras
                initialCollisions = self.bench.getCollisionsForCobra(
                    c, positions)
                initialCollisions += self.bench.getCollisionsForCobra(
                    nc, positions)

                # Free the current targets
                initialTarget1 = indices[c]
                initialTarget2 = indices[nc]
                freeTargets[initialTarget1] = True
                freeTargets[initialTarget2] = True

                # Get the targets that can be reached by each cobra
                targets1 = self.accessibleTargetIndices[c]
                targets1 = targets1[targets1 != NULL_TARGET_INDEX]
                targets2 = self.accessibleTargetIndices[nc]
                targets2 = targets2[targets2 != NULL_TARGET_INDEX]

                # Select only the free targets
                targets1 = targets1[freeTargets[targets1]]
                targets2 = targets2[freeTargets[targets2]]

                # Create an array with all the possible target combinations
                combinations = np.column_stack((
                    np.repeat(targets1, len(targets2)),
                    np.tile(targets2, len(targets1))))

                # Exclude the current target combination and combinations that
                # use the same target for the two cobras
                validCombinations = np.logical_or(
                    combinations[:, 0] != initialTarget1,
                    combinations[:, 1] != initialTarget2)
                validCombinations = np.logical_and(
                    validCombinations,
                    combinations[:, 0] != combinations[:, 1])
                combinations = combinations[validCombinations]

                # Loop over all the possible combinations until we find the
                # minimum number of collisions
                bestTarget1 = initialTarget1
                bestTarget2 = initialTarget2
                bestCollisions = initialCollisions

                for newTarget1, newTarget2 in combinations:
                    # Assign the new fiber positions
                    positions[c] = self.targets.positions[newTarget1]
                    positions[nc] = self.targets.positions[newTarget2]

                    # Calculate the number of collisions at the new positions
                    collisions = self.bench.getCollisionsForCobra(
                        c, positions)
                    collisions += self.bench.getCollisionsForCobra(
                        nc, positions)

                    # Check if the number of collisions decreased
                    if collisions < bestCollisions:
                        # Save the information from this target combination
                        bestTarget1 = newTarget1
                        bestTarget2 = newTarget2
                        bestCollisions = collisions

                        # Exit the loop if the number of collisions is zero
                        if bestCollisions == 0:
                            break

                # Do not use the best target combination if the decrease in the
                # number of collisions is only 1, because this means that we
                # solved the current collision, but we created a new collision
                # with another nearby cobra
                if (initialCollisions - bestCollisions) == 1:
                    bestTarget1 = initialTarget1
                    bestTarget2 = initialTarget2

                # Use the target combination where we had less collisions
                indices[c] = bestTarget1
                indices[nc] = bestTarget2
                positions[c] = self.targets.positions[bestTarget1]
                positions[nc] = self.targets.positions[bestTarget2]
                freeTargets[bestTarget1] = False
                freeTargets[bestTarget2] = False

    def getSelectedTargets(self):
        """Returns a new target group with the selected target for each cobra.

        This method should not be run before the selecTargets method.

        Returns
        -------
        object
            A new target group instance with the selected target for each cobra.

        """
        return self.targets.select(self.assignedTargetIndices)
