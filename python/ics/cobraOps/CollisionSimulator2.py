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

from procedures.moduleTest import engineer
from . import plotUtils


class CollisionSimulator2():
    """Class used to simulate a PFS observation.

    """

    def __init__(self, bench, cobraCoach, targets):
        """Constructs a new collision simulator instance.

        Parameters
        ----------
        bench: object
            The PFI bench instance.
        cobraCoach: object
            The cobra coach instance.
        targets: object
            The target group instance with the targets to observe.

        Returns
        -------
        object
            The collision simulator instance.

        """
        # Save the bench, the cobra coach and the target group instances
        self.bench = bench
        self.cobraCoach = cobraCoach
        self.targets = targets

        # Check which cobras are working fine and are not bad or broken
        self.nCobras = self.cobraCoach.nCobras
        self.goodCobras = np.full(self.nCobras, False)
        self.goodCobras[self.cobraCoach.goodIdx] = True

        # Check which cobras are assigned to a target
        self.assignedCobras = self.targets.notNull.copy()

        # Set the cobras that we are going to move
        self.movingCobras = np.logical_and(self.goodCobras, self.assignedCobras)

        # Define some internal variables that will filled by the run method
        self.finalFiberPositions = None
        self.trajectories = None
        self.fiberPositions = None
        self.elbowPositions = None
        self.nSteps = None
        self.associationCollisions = None
        self.associationEndPointCollisions = None
        self.collisions = None
        self.endPointCollisions = None
        self.nCollisions = None
        self.nEndPointCollisions = None

    def run(self, timeStep=20, maxSteps=2000):
        """Runs the collisions simulator.

        Parameters
        ----------
        timeStep: int, optional
            The trajectories time step resolution in steps. Default is 20
            steps.
        maxSteps: int, optional
            The trajectories maximum number of steps. Default is 2000 steps.

        """
        # Calculate the final fiber positions
        self.calculateFinalFiberPositions()

        # Calculate the cobra trajectories
        self.calculateTrajectories(timeStep, maxSteps)

        # Detect cobra collisions during the trajectory
        self.detectTrajectoryCollisions()

    def calculateFinalFiberPositions(self):
        """Calculates the cobras final fiber positions.

        """
        # Set the final fiber positions to their associated target positions,
        # with the exception of bad and/or unassigned cobras, which will be
        # moved to a position where they cannot collide with other cobras
        self.finalFiberPositions = self.cobraCoach.pfi.anglesToPositions(
            self.cobraCoach.allCobras, np.zeros(self.nCobras) - 0.00001,
            np.deg2rad(-170) - self.cobraCoach.calibModel.phiIn)
        self.finalFiberPositions[self.movingCobras] = self.targets.positions[
            self.movingCobras]

    def calculateTrajectories(self, timeStep, maxSteps):
        """Calculates the cobra trajectories.

        Parameters
        ----------
        timeStep: int
            The trajectories time step resolution in steps.
        maxSteps: int
            The trajectories maximum number of steps.

        """
        # Calculate the final theta and phi angles for the good cobras
        thetaAngles, phiAngles, _ = self.cobraCoach.pfi.positionsToAngles(
            self.cobraCoach.allCobras[self.goodCobras],
            self.finalFiberPositions[self.goodCobras])

        # Select the first angles solution
        thetaAngles = thetaAngles[:, 0]
        phiAngles = phiAngles[:, 0]

        # Initialize the engineer module
        engineer.setCobraCoach(self.cobraCoach)
        engineer.setConstantOntimeMode(maxSteps=maxSteps)

        # Calculate the cobra trajectories
        self.trajectories, _ = engineer.createTrajectory(
            np.where(self.goodCobras)[0], thetaAngles, phiAngles,
            tries=8, twoSteps=True, threshold=20.0, timeStep=timeStep)

        # Calculate the fiber and elbow positions along the cobra trajectories
        self.fiberPositions = self.trajectories.calculateFiberPositions(
            self.cobraCoach)
        self.elbowPositions = self.trajectories.calculateElbowPositions(
            self.cobraCoach)
        self.nSteps = self.fiberPositions.shape[1]

    def detectTrajectoryCollisions(self):
        """Detects collisions in the cobra trajectories.

        """
        # Extract some useful information from the bench instance
        cobraAssociations = self.bench.cobraAssociations
        linkRadius = self.bench.cobras.linkRadius

        # Calculate the distances between the cobras links for each step in the
        # trajectory
        startPoints1 = self.fiberPositions[cobraAssociations[0]].ravel()
        endPoints1 = self.elbowPositions[cobraAssociations[0]].ravel()
        startPoints2 = self.fiberPositions[cobraAssociations[1]].ravel()
        endPoints2 = self.elbowPositions[cobraAssociations[1]].ravel()
        distances = self.bench.distancesBetweenLineSegments(
            startPoints1, endPoints1, startPoints2, endPoints2)

        # Reshape the distances array
        distances = distances.reshape((len(cobraAssociations[0]), self.nSteps))

        # Detect trajectory collisions between cobra associations
        minimumSeparation = linkRadius[cobraAssociations[0]] + linkRadius[
            cobraAssociations[1]]
        trajectoryCollisions = distances < minimumSeparation[:, np.newaxis]

        # Check which cobra associations are affected by collisions
        self.associationCollisions = np.any(trajectoryCollisions, axis=1)
        self.associationEndPointCollisions = trajectoryCollisions[:, -1]

        # Check which cobras are involved in collisions
        collidingCobras = np.unique(
            self.bench.cobraAssociations[:, self.associationCollisions])
        self.collisions = np.full(self.nCobras, False)
        self.collisions[collidingCobras] = True
        self.nCollisions = np.sum(self.collisions)

        # Check which cobras are involved in end point collisions
        collidingCobras = np.unique(
            self.bench.cobraAssociations[:, self.associationEndPointCollisions])
        self.endPointCollisions = np.full(self.nCobras, False)
        self.endPointCollisions[collidingCobras] = True
        self.nEndPointCollisions = np.sum(self.endPointCollisions)

    def addTrajectories(self, colors=np.array([0.4, 0.4, 0.4, 1.0]),
                        indices=None, paintFootprints=False,
                        footprintColors=np.array([0.0, 0.0, 1.0, 0.05])):
        """Draws the cobra trajectories on top of an existing figure.

        Parameters
        ----------
        colors: object, optional
            The trajectory colors. Default is dark grey.
        indices: object, optional
            A numpy array with the cobra trajectory indices to use. If it is
            set to None, all the trajectories will be used. Default is None.
        paintFootprints: bool, optional
            If True, the cobra footprints will be painted. Default is False.
        footprintColors: object, optional
            The cobra footprints colors. Default is very light blue.

        """
        # Calculate some useful information
        fiberPositions = self.fiberPositions.copy()
        elbowPositions = self.elbowPositions.copy()
        linkRadius = self.bench.cobras.linkRadius

        # Select a subset of the trajectories if necessary
        if indices is not None:
            elbowPositions = elbowPositions[indices]
            fiberPositions = fiberPositions[indices]
            linkRadius = linkRadius[indices]

            if colors.ndim >= 2:
                colors = colors[indices]

            if footprintColors.ndim >= 2:
                footprintColors = footprintColors[indices]

        # Plot the elbow and fiber trajectories as continuous lines
        plotUtils.addTrajectories(
            np.vstack((elbowPositions, fiberPositions)),
            color=np.vstack((colors, colors)), linewidth=1)

        # Paint the cobra trajectory footprints if necessary
        if paintFootprints:
            # Calculate the line thicknesses
            thiknesses = np.empty(elbowPositions.shape)
            thiknesses[:] = linkRadius[:, np.newaxis]

            # Only use the elbow and fiber positions where the cobra is moving
            isMoving = np.empty(elbowPositions.shape, dtype="bool")
            isMoving[:, :-1] = (fiberPositions[:, 1:] - fiberPositions[:, :-1]) != 0
            isMoving[:, -1] = isMoving[:, -2]
            elbowPositions = elbowPositions[isMoving]
            fiberPositions = fiberPositions[isMoving]
            thiknesses = thiknesses[isMoving]

            # Update the colors if necessary
            if footprintColors.ndim > 2 and footprintColors.shape[:2] == isMoving.shape:
                # Set the colors for the moving positions
                footprintColors = footprintColors[isMoving]

                # Only use positions where the alpha color is not exactly zero
                visible = footprintColors[:, 3] != 0
                elbowPositions = elbowPositions[visible]
                fiberPositions = fiberPositions[visible]
                thiknesses = thiknesses[visible]
                footprintColors = footprintColors[visible]

            # Represent the trajectory footprint as a combination of thick
            # lines
            plotUtils.addThickLines(
                fiberPositions, elbowPositions, thiknesses,
                facecolor=footprintColors)

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
        plotUtils.createNewFigure(
            "Collision simulation results", "x position (mm)",
            "y position (mm)")

        # Set the axes limits
        limRange = 1.05 * self.bench.radius * np.array([-1, 1])
        xLim = self.bench.center.real + limRange
        yLim = self.bench.center.imag + limRange
        plotUtils.setAxesLimits(xLim, yLim)

        # Draw the cobra patrol areas
        patrolAreaColors = np.full((self.nCobras, 4), [0.0, 0.0, 1.0, 0.15])
        patrolAreaColors[self.collisions] = [1.0, 0.0, 0.0, 0.3]
        patrolAreaColors[self.endPointCollisions] = [0.0, 1.0, 0.0, 0.5]
        patrolAreaColors[~self.goodCobras] = [1.0, 0.0, 1.0, 0.2]
        self.bench.cobras.addPatrolAreasToFigure(colors=patrolAreaColors)

        # Draw the black dots
        self.bench.blackDots.addToFigure(colors=[0.0, 0.0, 0.0, 0.15])

        # Draw the cobra links at the final fiber positions
        linkColors = np.full(patrolAreaColors.shape, [0.0, 0.0, 1.0, 0.5])
        linkColors[~self.movingCobras] = [1.0, 0.0, 0.0, 0.25]
        linkColors[~self.goodCobras] = [1.0, 1.0, 1.0, 0.25]
        self.bench.cobras.addLinksToFigure(
            self.finalFiberPositions, colors=linkColors)

        # Draw the cobra trajectories and the trajectory footprints of those
        # that have a collision
        footprintColors = np.zeros((self.nCobras, self.nSteps, 4))
        footprintColors[self.collisions, :] = [0.0, 0.0, 1.0, 0.05]
        linkColors[self.collisions & ~self.goodCobras, :] = [1.0, 1.0, 1.0, 0.05]
        self.addTrajectories(
            paintFootprints=paintFootprints, footprintColors=footprintColors)

        # Draw the targets assigned to the cobras
        self.targets.addToFigure(colors=np.array([1.0, 0.0, 0.0, 1.0]))

        # Draw the extra targets if necessary
        if extraTargets is not None:
            # Draw only those targets that are not part of the simulation
            unusedTargets = np.logical_not(
                np.in1d(extraTargets.ids, self.targets.ids))
            extraTargets.addToFigure(indices=unusedTargets)
