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

from ics.cobraCharmer.cobraCoach import engineer

from . import plotUtils


class CollisionSimulator():
    """Class used to simulate a PFS observation.

    """

    def __init__(self, bench, targets):
        """Constructs a new collision simulator instance.

        Parameters
        ----------
        bench: object
            The PFI bench instance.
        targets: object
            The target group instance with the targets to observe.

        Returns
        -------
        object
            The collision simulator instance.

        """
        # Save the bench and the target group instances
        self.bench = bench
        self.targets = targets

        # Check which cobras are working correctly and are not bad or broken
        self.nCobras = self.bench.cobras.nCobras
        self.goodCobras = self.bench.cobras.isGood

        # Set the cobras that we are going to move
        self.movingCobras = np.logical_and(
            self.goodCobras, self.targets.notNull)

        # Define some internal variables that will be filled by the run method
        self.fiberPositions = None
        self.elbowPositions = None
        self.nSteps = None
        self.associationCollisions = None
        self.associationEndPointCollisions = None
        self.collisions = None
        self.endPointCollisions = None
        self.nCollisions = None
        self.nEndPointCollisions = None
        self.interferences = None
        self.nInterferences = None

    def run(self, timeStep=20, maxSteps=2000):
        """Runs the collisions simulator.

        Parameters
        ----------
        timeStep: int, optional
            The trajectories time step resolution in steps. Default is 20 steps.
        maxSteps: int, optional
            The trajectories maximum number of steps. Default is 2000 steps.

        """
        # Calculate the cobra trajectories
        self.calculateTrajectories(timeStep, maxSteps)

        # Detect cobra collisions while they move
        self.detectCollisions()

    def calculateTrajectories(self, timeStep, maxSteps):
        """Calculates the cobra trajectories.

        Parameters
        ----------
        timeStep: int
            The trajectories time step resolution in steps.
        maxSteps: int
            The trajectories maximum number of steps.

        """
        # Get the cobra coach instance from the bench
        cobraCoach = self.bench.cobras.cobraCoach

        # Calculate the theta and phi angles at the target positions
        thetaAngles, phiAngles, _ = cobraCoach.pfi.positionsToAngles(
            cobraCoach.allCobras[self.movingCobras],
            self.targets.positions[self.movingCobras])

        # Select the first angles solution
        thetaAngles = thetaAngles[:, 0]
        phiAngles = phiAngles[:, 0]

        # Initialize the engineer module
        engineer.setCobraCoach(cobraCoach)
        engineer.setConstantOntimeMode(maxSteps=maxSteps)

        # Calculate the cobra trajectories
        trajectories, _ = engineer.createTrajectory(
            np.where(self.movingCobras)[0], thetaAngles, phiAngles,
            tries=8, twoSteps=True, threshold=20.0, timeStep=timeStep)

        # Calculate the fiber and elbow positions along the cobra trajectories
        self.fiberPositions = trajectories.calculateFiberPositions(cobraCoach)
        self.elbowPositions = trajectories.calculateElbowPositions(cobraCoach)
        self.nSteps = self.fiberPositions.shape[1]

    def detectCollisions(self):
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
        distances = CollisionSimulator.distancesBetweenLineSegments(
            startPoints1, endPoints1, startPoints2, endPoints2)

        # Reshape the distances array
        distances = distances.reshape((cobraAssociations.shape[1], self.nSteps))

        # Detect trajectory collisions between cobra associations
        minimumSeparation = linkRadius[cobraAssociations[0]] + linkRadius[
            cobraAssociations[1]]
        trajectoryCollisions = distances < minimumSeparation[:, np.newaxis]

        # Check which cobra associations are affected by collisions
        self.associationCollisions = np.any(trajectoryCollisions, axis=1)
        self.associationEndPointCollisions = trajectoryCollisions[:, -1]

        # Check which cobras are involved in collisions
        collidingCobras = np.unique(
            cobraAssociations[:, self.associationCollisions])
        self.collisions = np.full(self.nCobras, False)
        self.collisions[collidingCobras] = True
        self.nCollisions = len(collidingCobras)

        # Check which cobras are involved in end point collisions
        collidingCobras = np.unique(
            cobraAssociations[:, self.associationEndPointCollisions])
        self.endPointCollisions = np.full(self.nCobras, False)
        self.endPointCollisions[collidingCobras] = True
        self.nEndPointCollisions = len(collidingCobras)

        # Check which good cobras could collide with fiducial fibers
        cobraCoach = self.bench.cobras.cobraCoach
        thetaAngles, phiAngles, _ = cobraCoach.pfi.positionsToAngles(
            cobraCoach.goodCobras, self.fiberPositions[self.goodCobras, -1])
        thetaAngles = thetaAngles[:, 0]
        phiAngles = phiAngles[:, 0]
        interferenceCobrasIndices = np.array(
            cobraCoach.checkFiducialInterference(
                thetaAngles, phiAngles), dtype=int)
        self.interferences = np.full(self.nCobras, False)
        self.interferences[interferenceCobrasIndices] = True
        self.nInterferences = self.interferences.sum()

    @staticmethod
    def distancesBetweenLineSegments(startPoints1, endPoints1, startPoints2,
                                     endPoints2):
        """Calculates the minimum distances between two sets of line segments.

        Parameters
        ----------
        startPoints1: object
            A complex numpy array with the first set of line segments start
            coordinates.
        endPoints1: object
            A complex numpy array with the first set of line segments end
            coordinates.
        startPoints2: object
            A complex numpy array with the second set of line segments start
            coordinates.
        endPoints2: object
            A complex numpy array with the second set of line segments end
            coordinates.

        Returns
        -------
        object
            A numpy array with the minimum distance between the two sets of line
            segments.

        """
        # Calculate the minimum distances for each point to segment combination
        distances1 = CollisionSimulator.distancesToLineSegments(
            startPoints1, startPoints2, endPoints2)
        distances2 = CollisionSimulator.distancesToLineSegments(
            endPoints1, startPoints2, endPoints2)
        distances3 = CollisionSimulator.distancesToLineSegments(
            startPoints2, startPoints1, endPoints1)
        distances4 = CollisionSimulator.distancesToLineSegments(
            endPoints2, startPoints1, endPoints1)

        # Return the minimum distances
        return np.min((distances1, distances2, distances3, distances4), axis=0)

    @staticmethod
    def distancesToLineSegments(points, startPoints, endPoints):
        """Calculates the minimum distances between a set of points and a set of
        line segments.

        Parameters
        ----------
        points: object
            A complex numpy array with the point coordinates.
        startPoints: object
            A complex numpy array with the line segments start coordinates.
        endPoints: object
            A complex numpy array with the line segments end coordinates.

        Returns
        -------
        object
            A numpy array with the minimum distances between the set of points
            and the set of line segments.

        """
        # Translate the points and the line segment end points to the line
        # segment starting points
        translatedPoints = points - startPoints
        translatedEndPoints = endPoints - startPoints

        # Rotate the translated points to have the line segment on the x axis
        rotatedPoints = translatedPoints * np.exp(-1j * np.angle(
            translatedEndPoints))

        # Define 3 regions for the points: left of the origin, over the line
        # segments, and right of the line segments
        x = rotatedPoints.real
        lineLengths = np.abs(translatedEndPoints)
        (region1,) = np.where(x <= 0)
        (region2,) = np.where(np.logical_and(x > 0 , x < lineLengths))
        (region3,) = np.where(x >= lineLengths)

        # Calculate the minimum distances in each region
        distances = np.empty(len(points))
        distances[region1] = np.abs(rotatedPoints[region1])
        distances[region2] = np.abs(rotatedPoints[region2].imag)
        distances[region3] = np.abs(
            rotatedPoints[region3] - lineLengths[region3])

        return distances

    def addTrajectoriesToFigure(self, colors=np.array([0.4, 0.4, 0.4, 1.0]),
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
        # Extract some useful information
        fiberPositions = self.fiberPositions
        elbowPositions = self.elbowPositions
        linkRadius = self.bench.cobras.linkRadius

        # Select a subset of the trajectories if necessary
        if indices is not None:
            fiberPositions = fiberPositions[indices]
            elbowPositions = elbowPositions[indices]
            linkRadius = linkRadius[indices]

            if colors.ndim >= 2:
                colors = colors[indices]

            if footprintColors.ndim >= 2:
                footprintColors = footprintColors[indices]

        # Plot the elbow and fiber trajectories as continuous lines
        plotUtils.addTrajectories(
            np.vstack((fiberPositions, elbowPositions)),
            color=np.vstack((colors, colors)), linewidth=1)

        # Paint the cobra trajectory footprints if necessary
        if paintFootprints:
            # Calculate the line thicknesses
            thiknesses = np.empty(fiberPositions.shape)
            thiknesses[:] = linkRadius[:, np.newaxis]

            # Only use the elbow and fiber positions where the cobra is moving
            isMoving = np.empty(fiberPositions.shape, dtype="bool")
            isMoving[:, :-1] = (
                fiberPositions[:, 1:] - fiberPositions[:, :-1]) != 0
            isMoving[:, -1] = isMoving[:, -2]
            fiberPositions = fiberPositions[isMoving]
            elbowPositions = elbowPositions[isMoving]
            thiknesses = thiknesses[isMoving]

            # Update the colors if necessary
            if footprintColors.ndim > 2 and footprintColors.shape[:2] == isMoving.shape:
                # Set the colors for the moving positions
                footprintColors = footprintColors[isMoving]

                # Only use positions where the alpha color is not exactly zero
                isVisible = footprintColors[:, 3] != 0
                fiberPositions = fiberPositions[isVisible]
                elbowPositions = elbowPositions[isVisible]
                thiknesses = thiknesses[isVisible]
                footprintColors = footprintColors[isVisible]

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

        # Draw the fiducial fibers if necessary
        if self.bench.fiducials is not None:
            fiducials = self.bench.fiducials
            nFiducials = len(fiducials)
            fiducialPositions = (
                fiducials["x_mm"] + 1j * fiducials["y_mm"]).to_numpy()
            rMin = np.full(nFiducials, 1)
            rMax = np.full(nFiducials, 2)
            plotUtils.addRings(
                fiducialPositions, rMin, rMax, facecolors=[1.0, 0.0, 0.0, 0.75])

        # Draw the cobra patrol areas
        patrolAreaColors = np.full((self.nCobras, 4), [0.0, 0.0, 1.0, 0.15])
        patrolAreaColors[self.collisions] = [1.0, 0.0, 0.0, 0.3]
        patrolAreaColors[self.endPointCollisions] = [0.0, 1.0, 0.0, 0.5]
        patrolAreaColors[~self.goodCobras] = [0.0, 0.0, 0.0, 0.25]
        patrolAreaColors[self.interferences] = [0.0, 0.0, 0.0, 0.5]
        self.bench.cobras.addPatrolAreasToFigure(colors=patrolAreaColors)

        # Draw the black dots
        self.bench.blackDots.addToFigure(colors=[0.0, 0.0, 0.0, 0.15])

        # Draw the cobra links at the final fiber positions
        linkColors = np.full((self.nCobras, 4), [0.0, 0.0, 1.0, 0.5])
        linkColors[~self.movingCobras] = [1.0, 0.0, 0.0, 0.25]
        linkColors[~self.goodCobras] = [0.0, 0.0, 0.0, 0.25]
        self.bench.cobras.addLinksToFigure(
            self.fiberPositions[:, -1], self.elbowPositions[:,-1],
            colors=linkColors)

        # Draw the cobra trajectories and the trajectory footprints of those
        # that have a collision
        footprintColors = np.zeros((self.nCobras, self.nSteps, 4))
        footprintColors[self.collisions, :] = [0.0, 0.0, 1.0, 0.05]
        self.addTrajectoriesToFigure(
            paintFootprints=paintFootprints, footprintColors=footprintColors)

        # Draw the targets assigned to the cobras
        self.targets.addToFigure(colors=np.array([1.0, 0.0, 0.0, 1.0]))

        # Draw the extra targets if necessary
        if extraTargets is not None:
            # Draw only those targets that are not part of the simulation
            unusedTargets = ~np.in1d(extraTargets.ids, self.targets.ids)
            extraTargets.addToFigure(indices=unusedTargets)
