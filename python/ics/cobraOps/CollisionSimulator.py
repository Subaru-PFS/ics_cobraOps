"""

CollisionSimulator class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from . import plotUtils

from .TrajectoryGroup import TrajectoryGroup

#from IPython.core.debugger import Tracer

class CollisionSimulator():
    """

    Class used to simulate a PFS observation.

    """

    def __init__(self, bench, targets, trajectorySteps=170, trajectoryStepWidth=50):
        """Constructs a new collision simulator instance.

        Parameters
        ----------
        bench: object
            The PFI bench instance.
        targets: object
            The target group instance.
        nSteps: int, optional
            The total number of steps in the cobra trajectories. Default is
            150.
        stepWidth: int, optional
            The trajectory step width in units of motor steps. Default is 50.

        Returns
        -------
        object
            The collision simulator instance.

        """
        # Save the bench and target group instances
        self.bench = bench
        self.targets = targets

        # Save the trajectory parameters
        self.trajectorySteps = trajectorySteps
        self.trajectoryStepWidth = trajectoryStepWidth

        # Check which cobras are assigned to a target
        self.assignedCobras = self.targets.notNull.copy()

        # Define some internal variables that will filled by the run method
        self.finalFiberPositions = None
        self.posSteps = None
        self.negSteps = None
        self.movementDirections = None
        self.movementStrategies = None
        self.trajectories = None
        self.associationCollisions = None
        self.associationEndPointCollisions = None
        self.collisions = None
        self.endPointCollisions = None
        self.nCollisions = None
        self.nEndPointCollisions = None


    def run(self, solveCollisions=True):
        """Runs the collisions simulator.

        Parameters
        ----------
        solveCollisions: bool, optional
            If True, the simulator will try to solve trajectory collisions,
            changing the movement directions and strategies of the affected
            cobras. Default is True.

        """
        # Calculate the final fiber positions
        self.calculateFinalFiberPositions()

        # Define the theta and phi movement directions
        self.defineMovementDirections()

        # Define the theta and phi movement strategies
        self.defineMovementStrategies()

        # Calculate the cobra trajectories
        self.calculateTrajectories()

        # Detect cobra collisions during the trajectory
        self.detectTrajectoryCollisions()

        # Solve trajectory collisions if requested
        if solveCollisions:
            self.solveTrajectoryCollisions(True)
            self.solveTrajectoryCollisions(False)
            self.solveTrajectoryCollisions(True)
            self.solveTrajectoryCollisions(False)


    def calculateFinalFiberPositions(self):
        """Calculates the cobras final fiber positions.

        """
        # Set the final fiber positions to their associated target positions,
        # leaving unassigned cobras at their home positive positions
        self.finalFiberPositions = self.bench.cobras.home0.copy()
        self.finalFiberPositions[self.assignedCobras] = self.targets.positions[self.assignedCobras]

        # Optimize the unassigned cobra positions to minimize their possible
        # collisions with other cobras
        self.optimizeUnassignedCobraPositions()


    def optimizeUnassignedCobraPositions(self):
        """Finds the unassigned cobras final fiber positions that minimize
        their collisions with other cobras.

        """
        # Calculate the cobras elbow positions at their current final fiber
        # positions
        elbowPositions = self.bench.cobras.calculateElbowPositions(self.finalFiberPositions)

        # Get the unassigned cobra indices
        (unassigendCobraIndices,) = np.where(self.assignedCobras == False)

        # Find the optimal position for each unassigned cobra
        for c in unassigendCobraIndices:
            # Get the cobra nearest neighbors
            cobraNeighbors = self.bench.getCobraNeighbors(c)

            # Select only those that are assigned to a target
            cobraNeighbors = cobraNeighbors[self.assignedCobras[cobraNeighbors]]

            # Jump to the next cobra if all the neighbors are unassigned
            if len(cobraNeighbors) == 0:
                continue

            # Calculate all the possible cobra elbow rotations
            rotationAngles = np.arange(0, 2 * np.pi, 0.1 * np.pi)
            center = self.bench.cobras.centers[c]
            rotatedElbowPositions = (elbowPositions[c] - center) * np.exp(1j * rotationAngles) + center

            # Obtain the angle that maximizes the closer distance to a neighbor
            distances = np.abs(self.finalFiberPositions[cobraNeighbors] - rotatedElbowPositions[:, np.newaxis])
            minDistances = np.min(distances, axis=1)
            optimalAngle = rotationAngles[np.argmax(minDistances)]

            # Update the cobra final fiber position
            self.finalFiberPositions[c] = (self.finalFiberPositions[c] - center) * np.exp(1j * optimalAngle) + center


    def defineMovementDirections(self):
        """Defines the theta and phi movement directions that the cobras should
        follow.

        The movement directions are encoded in a boolean array with 2 rows and
        as many columns as cobras in the array:
            - The first row indicates the theta movement direction: positive
            for True values, and negative for False values.
            - The second row indicates the phi movement direction: positive
            for True values, and negative for False values.

        """
        # Get the cobra rotation angles for the starting positive and negative
        # home positions and the final fiber positions
        (posStartTht, posStartPhi) = self.bench.cobras.calculateRotationAngles(self.bench.cobras.home0)
        (negStartTht, negStartPhi) = self.bench.cobras.calculateRotationAngles(self.bench.cobras.home1)
        (finalTht, finalPhi) = self.bench.cobras.calculateRotationAngles(self.finalFiberPositions)

        # Calculate the required theta and phi delta offsets to move from the
        # positive and negative home positions to the final positions
        posDeltaTht = np.mod(finalTht - posStartTht, 2 * np.pi)
        negDeltaTht = -np.mod(negStartTht - finalTht, 2 * np.pi)
        posDeltaPhi = finalPhi - posStartPhi
        negDeltaPhi = finalPhi - negStartPhi

        # Calculate the total number of motor steps required to reach the final
        # positions from the positive and the negative home positions
        (posThtMotorSteps, posPhiMotorSteps) = self.bench.cobras.motorMaps.calculateSteps(posDeltaTht, posStartPhi, posDeltaPhi)
        (negThtMotorSteps, negPhiMotorSteps) = self.bench.cobras.motorMaps.calculateSteps(negDeltaTht, negStartPhi, negDeltaPhi)

        # Calculate the number of trajectory steps required to reach the final
        # positions from the positive and the negative starting positions
        self.posSteps = np.ceil(np.max((posThtMotorSteps, posPhiMotorSteps), axis=0) / self.trajectoryStepWidth).astype("int") + 1
        self.negSteps = np.ceil(np.max((negThtMotorSteps, negPhiMotorSteps), axis=0) / self.trajectoryStepWidth).astype("int") + 1

        # Make sure that at least one of the movements requires less steps than
        # the maximum number of steps allowed
        if np.any(np.min((self.posSteps, self.negSteps), axis=0) > self.trajectorySteps):
            raise Exception("Some cobras cannot reach their assigned targets "
                            "because the trajectorySteps parameter value is "
                            "too low. Please set it to a higher value.")

        # Decide if the cobras should follow a positive theta movement
        # direction:
        #    - The positive movement should require less steps than the
        #    negative movement.
        posThtMovement = self.posSteps < self.negSteps

        # Calculate the phi movement direction
        posPhiMovement = negDeltaPhi > 0
        posPhiMovement[posThtMovement] = posDeltaPhi[posThtMovement] > 0

        # Save the results in a single array
        self.movementDirections = np.vstack((posThtMovement, posPhiMovement))


    def defineMovementStrategies(self):
        """Defines the theta and phi movement strategies that the cobras should
        follow.

        There are only three strategies that really make sense:
            - Early theta and early phi movements, when phi is moving towards
            the cobra center.
            - Early theta and late phi movements, when the theta movement is in
            the positive direction.
            - Late theta and late phi movements, when the theta movement is in
            the negative direction.

        The movement strategies are encoded in a boolean array with 2 rows and
        as many columns as cobras in the array:
            - The first row indicates the theta movement strategy: "early" for
            True values, and "late" for False values.
            - The second row indicates the phi movement strategy: "early" for
            True values, and "late" for False values.

        """
        # By default always do late theta and late phi movements
        thtEarly = np.full(self.bench.cobras.nCobras, False)
        phiEarly = thtEarly.copy()

        # Do early theta movements when moving in the theta positive direction
        thtEarly[self.movementDirections[0]] = True

        # If the cobra is moving towards the center, do the theta and phi
        # movements as early as possible
        towardsTheCenter = np.logical_not(self.movementDirections[1])
        thtEarly[towardsTheCenter] = True
        phiEarly[towardsTheCenter] = True

        # Do early movements for unassigned cobras
        thtEarly[self.assignedCobras == False] = True
        phiEarly[self.assignedCobras == False] = True

        # Save the results in a single array
        self.movementStrategies = np.vstack((thtEarly, phiEarly))


    def calculateTrajectories(self):
        """Calculates the cobra trajectories.

        """
        self.trajectories = TrajectoryGroup(nSteps=self.trajectorySteps,
                                            stepWidth=self.trajectoryStepWidth,
                                            bench=self.bench,
                                            finalFiberPositions=self.finalFiberPositions,
                                            movementDirections=self.movementDirections,
                                            movementStrategies=self.movementStrategies)


    def detectTrajectoryCollisions(self):
        """Detects collisions in the cobra trajectories.

        """
        # Detect trajectory collisions between cobra associations
        trajectoryCollisions, self.distances = self.trajectories.calculateCobraAssociationCollisions()

        # Check which are the cobra associations affected by collisions
        self.associationCollisions = np.any(trajectoryCollisions, axis=1)
        self.associationEndPointCollisions = trajectoryCollisions[:, -1]

        # Check which cobras are involved in collisions
        collidingCobras = np.unique(self.bench.cobraAssociations[:, self.associationCollisions])
        self.collisions = np.full(self.bench.cobras.nCobras, False)
        self.collisions[collidingCobras] = True
        self.nCollisions = np.sum(self.collisions)

        # Check which cobras are involved in end point collisions
        collidingCobras = np.unique(self.bench.cobraAssociations[:, self.associationEndPointCollisions])
        self.endPointCollisions = np.full(self.bench.cobras.nCobras, False)
        self.endPointCollisions[collidingCobras] = True
        self.nEndPointCollisions = np.sum(self.endPointCollisions)


    def solveTrajectoryCollisions(self, selectLowerIndices):
        """Solves trajectory collisions changing the cobras theta movement
        directions.

        Parameters
        ----------
        selectLowerIndices: bool
            If True, the cobras selected to be changed in the association will
            be the ones with the lower indices values.

        """
        # Get the indices of the cobras involved in a mid point trajectory
        # collision
        associationMidPointCollisions = np.logical_and(self.associationCollisions, self.associationEndPointCollisions == False)
        collidingAssociations = self.bench.cobraAssociations[:, associationMidPointCollisions]

        # Select the indices of the cobras whose movement should be changed
        cobraIndices = np.unique(collidingAssociations[0] if selectLowerIndices else collidingAssociations[1])

        # Make sure that there is at least one cobra to change
        if len(cobraIndices) > 0:
            # Change the cobras theta movement directions
            self.movementDirections[0, cobraIndices] = np.logical_not(self.movementDirections[0, cobraIndices])

            # Make sure that the new movement directions don't require too many
            # steps
            self.movementDirections[0, self.posSteps > self.trajectories.nSteps] = False
            self.movementDirections[0, self.negSteps > self.trajectories.nSteps] = True

            # Define the theta and phi movement strategies
            self.defineMovementStrategies()

            # Calculate the cobra trajectories
            self.calculateTrajectories()

            # Recalculate the cobra collisions during the trajectory
            self.recalculateTrajectoryCollisions(cobraIndices)


    def recalculateTrajectoryCollisions(self, cobraIndices):
        """Recalculates the trajectory collision information for the given
        cobras.

        Parameters
        ----------
        cobraIndices: object
            A numpy array with the indices of the cobras whose movement has
            changed.

        """
        # Get the cobra associations for the given cobras
        cobraAssociationIndices = np.in1d(self.bench.cobraAssociations[0], cobraIndices)
        cobraAssociationIndices = np.logical_or(cobraAssociationIndices, np.in1d(self.bench.cobraAssociations[1], cobraIndices))

        # Detect trajectory collisions between these cobra associations
        trajectoryCollisions, self.distances = self.trajectories.calculateCobraAssociationCollisions(cobraAssociationIndices)

        # Update the cobra associations affected by collisions
        self.associationCollisions[cobraAssociationIndices] = np.any(trajectoryCollisions, axis=1)
        self.associationEndPointCollisions[cobraAssociationIndices] = trajectoryCollisions[:, -1]

        # Check which cobras are involved in collisions
        collidingCobras = np.unique(self.bench.cobraAssociations[:, self.associationCollisions])
        self.collisions[:] = False
        self.collisions[collidingCobras] = True
        self.nCollisions = np.sum(self.collisions)

        # Check which cobras are involved in end point collisions
        collidingCobras = np.unique(self.bench.cobraAssociations[:, self.associationEndPointCollisions])
        self.endPointCollisions[:] = False
        self.endPointCollisions[collidingCobras] = True
        self.nEndPointCollisions = np.sum(self.endPointCollisions)


    def plotResults(self, extraTargets=None, paintFootprints=False):
        """Plots the collision simulator results in a new figure.

        Parameters
        ----------
        extraTargets: object, optional
            Extra targets that should also be plotted in the figure, in
            addition to the targets that were used in the simulation. Default
            is None.
        paintFootprints: bool, optional
            If True, the cobra trajectory footprints will be painted. Default
            is False.

        """
        # Create a new figure
        plotUtils.createNewFigure("Collision simulation results", "x position (mm)", "y position (mm)")

        # Set the axes limits
        limRange = 1.05 * self.bench.radius * np.array([-1, 1])
        xLim = self.bench.center.real + limRange
        yLim = self.bench.center.imag + limRange
        plotUtils.setAxesLimits(xLim, yLim)

        # Draw the cobra patrol areas
        patrolAreaColors = np.full((self.bench.cobras.nCobras, 4), [0.0, 0.0, 1.0, 0.15])
        patrolAreaColors[self.collisions] = [1.0, 0.0, 0.0, 0.3]
        patrolAreaColors[self.endPointCollisions] = [0.0, 1.0, 0.0, 0.5]
        self.bench.cobras.addPatrolAreasToFigure(colors=patrolAreaColors)

        # Draw the cobra links at the final fiber positions
        linkColors = np.full(patrolAreaColors.shape, [0.0, 0.0, 1.0, 0.5])
        linkColors[self.assignedCobras == False] = [1.0, 0.0, 0.0, 0.25]
        self.bench.cobras.addLinksToFigure(self.finalFiberPositions, colors=linkColors)

        # Draw the cobra trajectories and the trajectory footprints of those
        # that have a collision
        footprintColors = np.zeros((self.bench.cobras.nCobras, self.trajectories.nSteps, 4))
        footprintColors[self.collisions, :] = [0.0, 0.0, 1.0, 0.05]
        self.trajectories.addToFigure(paintFootprints=paintFootprints, footprintColors=footprintColors)

        # Draw the targets assigned to the cobras
        self.targets.addToFigure(colors=np.array([1.0, 0.0, 0.0, 1.0]))

        # Draw the extra targets if necessary
        if extraTargets is not None:
            # Draw only those targets that are not part of the simulation
            unusedTargets = np.logical_not(np.in1d(extraTargets.ids, self.targets.ids))
            extraTargets.addToFigure(indices=unusedTargets)


    def animateCobraTrajectory(self, cobraIndex, extraTargets=None, fileName=None):
        """Animates the trajectory of a given cobra and its nearest neighbors.

        Parameters
        ----------
        cobraIndex: int
            The index of the cobra to animate.
        extraTargets: object, optional
            Extra targets that should also be plotted in the figure, in
            addition to the targets that were used in the simulation. Default
            is None.
        fileName: object, optional
            The file name path where a video of the animation should be saved.
            If it is set to None, no video will be saved. Default is None.

        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        cobraCenters = self.bench.cobras.centers
        rMax = self.bench.cobras.rMax

        # Create a new figure
        plotUtils.createNewFigure("Trajectory animation for cobra " + str(cobraIndex), "x position (mm)", "y position (mm)")

        # Get the cobra neighbors
        cobraNeighbors = self.bench.getCobraNeighbors(cobraIndex)

        # Set the axes limits
        distances = np.abs(cobraCenters[cobraNeighbors] - cobraCenters[cobraIndex])
        limRange = 1.05 * np.max(distances + rMax[cobraNeighbors]) * np.array([-1, 1])
        xLim = cobraCenters[cobraIndex].real + limRange
        yLim = cobraCenters[cobraIndex].imag + limRange
        plotUtils.setAxesLimits(xLim, yLim)

        # Draw the cobra patrol areas
        patrolAreaColors = np.full((nCobras, 4), [0.0, 0.0, 1.0, 0.15])
        patrolAreaColors[self.collisions] = [1.0, 0.0, 0.0, 0.3]
        patrolAreaColors[self.endPointCollisions] = [0.0, 1.0, 0.0, 0.5]
        self.bench.cobras.addPatrolAreasToFigure(colors=patrolAreaColors)

        # Select which cobras should be animated: only those that fall inside
        # the displayed area
        toAnimate = np.full(nCobras, False)
        toAnimate[cobraIndex] = True
        toAnimate[cobraNeighbors] = True
        toAnimate[self.bench.getCobrasNeighbors(cobraNeighbors)] = True

        # Draw the cobras that should not be animated at their final positions
        linkColors = np.full((nCobras, 4), [0.0, 0.0, 1.0, 0.5])
        linkColors[self.assignedCobras == False] = [1.0, 0.0, 0.0, 0.25]
        self.bench.cobras.addLinksToFigure(self.finalFiberPositions, colors=linkColors, indices=np.logical_not(toAnimate))

        # Draw every point in the elbow and fiber trajectories
        plotUtils.addPoints(self.trajectories.elbowPositions.ravel(), s=2, facecolor=[1.0, 1.0, 1.0, 1.0])
        plotUtils.addPoints(self.trajectories.fiberPositions.ravel(), s=2, facecolor=[1.0, 1.0, 1.0, 1.0])

        # Draw the targets assigned to the cobras
        self.targets.addToFigure(colors=np.array([1.0, 0.0, 0.0, 1.0]))

        # Draw the extra targets if necessary
        if extraTargets is not None:
            # Draw only those targets that are not part of the simulation
            unusedTargets = np.logical_not(np.in1d(extraTargets.ids, self.targets.ids))
            extraTargets.addToFigure(indices=unusedTargets)

        # Add the animation
        self.trajectories.addAnimationToFigure(linkColors=linkColors, indices=toAnimate, fileName=fileName)
