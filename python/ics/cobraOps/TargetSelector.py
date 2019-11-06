"""

TargetSelector class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from .cobraConstants import NULL_TARGET_INDEX


class TargetSelector():
    """

    Class used to select optimal targets for a given PFI bench.

    """

    def __init__(self, bench, targets):
        """Constructs a new target selector instance.

        Parameters
        ----------
        bench: object
            The PFI bench instance.
        targets: object
            The target group instance.

        Returns
        -------
        object
            The target selector instance.

        """
        # Save the bench and target group instances
        self.bench = bench
        self.targets = targets

        # Define some internal arrays that will be used and filled by the
        # computeAccessibleTargets and selectTargets methods
        self.accessibleTargetIndices = None
        self.accessibleTargetDistances = None
        self.accessibleTargetElbows = None
        self.assignedTargetIndices = None


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
        # Do nothing. This method should be implemented by the subclass
        return


    def calculateAccessibleTargets(self, maximumDistance=np.Inf):
        """Calculates the targets that each cobra can reach.

        The results are saved in the accessibleTargetIndices and
        accesssibleTargetDistances internal arrays.

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
        rMin = self.bench.cobras.rMin
        rMax = self.bench.cobras.rMax
        targetPositions = self.targets.positions
        notNullTargets = self.targets.notNull
        cobraBlackDotPosition = self.bench.cobras.blackDotPosition
        cobraBlackDotRadius = self.bench.cobras.blackDotRadius

        # Calculate the maximum target distance allowed for each cobra
        rMax = rMax.copy()
        rMax[rMax > maximumDistance] = maximumDistance

        # Obtain the cobra-target associations: select first by the x axis
        # distance and then by the y axis distance. Remember to mask any
        # possible NULL target that might exist.
        xDistanceMatrix = np.abs(cobraCenters.real[:, np.newaxis] - targetPositions.real)
        xDistanceMatrix[:, notNullTargets == False] = np.Inf
        (cobraIndices, targetIndices) = np.where(xDistanceMatrix < rMax[:, np.newaxis])
        yDistance = np.abs(cobraCenters[cobraIndices].imag - targetPositions[targetIndices].imag)
        validIndices = yDistance < rMax[cobraIndices]
        cobraIndices = cobraIndices[validIndices]
        targetIndices = targetIndices[validIndices]

        # Select only those targets that can be reached by each cobra and avoid
        # targets falling in the cobra black dots
        distances = np.abs(cobraCenters[cobraIndices] - targetPositions[targetIndices])
        blackDotdistances = np.abs(cobraCenters[cobraIndices] + cobraBlackDotPosition[cobraIndices] - targetPositions[targetIndices])
        validIndices = np.logical_and(distances > rMin[cobraIndices], distances < rMax[cobraIndices])
        validIndices = np.logical_and(validIndices, blackDotdistances > cobraBlackDotRadius[cobraIndices])
        cobraIndices = cobraIndices[validIndices]
        targetIndices = targetIndices[validIndices]
        distances = distances[validIndices]

        # Calculate the elbow positions for all the cobra-target associations
        elbows = self.bench.cobras.calculateMultipleElbowPositions(targetPositions, cobraIndices, targetIndices)

        # Calculate the total number of targets that each cobra can reach
        nTargetsPerCobra = np.bincount(cobraIndices)

        # Create the accessible target arrays
        maxTagetsPerCobra = nTargetsPerCobra.max()
        self.accessibleTargetIndices = np.full((nCobras, maxTagetsPerCobra), NULL_TARGET_INDEX, dtype="int")
        self.accessibleTargetDistances = np.zeros((nCobras, maxTagetsPerCobra))
        self.accessibleTargetElbows = np.zeros((nCobras, maxTagetsPerCobra), dtype="complex")

        # Fill the arrays ordering the targets by their distance to the cobra
        counter = 0

        for i in range(len(nTargetsPerCobra)):
            # Get the target indices and distances for this cobra
            nTargetsForThisCobra = nTargetsPerCobra[i]
            targetsForThisCobra = targetIndices[counter:counter + nTargetsForThisCobra]
            distancesForThisCobra = distances[counter:counter + nTargetsForThisCobra]
            elbowsForThisCobra = elbows[counter:counter + nTargetsForThisCobra]

            # Sort the targets by their distance to the cobra and fill the arrays
            sortedIndices = distancesForThisCobra.argsort()
            self.accessibleTargetIndices[i, :nTargetsForThisCobra] = targetsForThisCobra[sortedIndices]
            self.accessibleTargetDistances[i, :nTargetsForThisCobra] = distancesForThisCobra[sortedIndices]
            self.accessibleTargetElbows[i, :nTargetsForThisCobra] = elbowsForThisCobra[sortedIndices]

            # Increase the counter
            counter += nTargetsForThisCobra


    def selectTargets(self):
        """Selects a single target for each cobra.

        This method should always be run after the calculateAccessibleTargets
        method.

        """
        # Do nothing. This method should be implemented by the subclass
        return


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
        fiberPositions[usedCobras] = targetPositions[self.assignedTargetIndices[usedCobras]]

        # Get the indices of the cobra associations where we have a collision
        problematicCobraAssociations = self.bench.getProblematicCobraAssociations(fiberPositions)
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
