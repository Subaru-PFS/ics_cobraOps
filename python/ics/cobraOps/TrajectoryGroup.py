"""

Trajectory group class.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

import ics.cobraOps.plotUtils as plotUtils

from ics.cobraOps.AttributePrinter import AttributePrinter
from ics.cobraOps.cobraConstants import (TRAJECTORY_BIN_WIDTH,
                                         TRAJECTORY_MAX_STEPS)


class TrajectoryGroup(AttributePrinter):
    """
    
    Class describing the properties of a group of cobra trajectories.
    
    """
    
    def __init__(self, bench, finalFiberPositions, movementDirections, movementStrategies):
        """Constructs a new trajectory group instance.
        
        Parameters
        ----------
        bench: object
            The PFI bench instance.
        finalFiberPositions: object
            A complex numpy array with the trajectory final fiber positions for
            each cobra.
        movementDirections: object
            A boolean numpy array with the theta and phi movement directions to
            use. True values indicate that the cobras should move in the
            positive theta or phi directions, while False values indicate that
            the movement should be in the negative direction.
        movementStrategies: object
            A boolean numpy array with the theta and phi movement strategies to
            use. True values indicate that the cobras should move in those
            angles as soon as possible, while False values indicate that the
            angle movement should be as late as possible.
        
        Returns
        -------
        object
            The trajectory group instance.
        
        """
        # Save the bench instance, the final positions and the movement arrays
        self.bench = bench
        self.finalFiberPositions = finalFiberPositions.copy()
        self.movementDirections = movementDirections.copy()
        self.movementStrategies = movementStrategies.copy()
        
        # Fix the number of trajectory steps to the maximum allowed
        self.nSteps = TRAJECTORY_MAX_STEPS
        
        # Calculate the trajectory stating fiber positions
        self.calculateStartingFiberPositions()
        
        # Calculate the cobra trajectories
        self.calculateCobraTrajectories()
    
    
    def calculateStartingFiberPositions(self):
        """Calculates the trajectories starting fiber positions.
        
        """
        # Set the start positions according to the specified movement direction
        self.startFiberPositions = self.bench.cobras.home1.copy()
        self.startFiberPositions[self.movementDirections[0]] = self.bench.cobras.home0[self.movementDirections[0]]
    
    
    def calculateCobraTrajectories(self):
        """Calculates the cobra trajectories.
        
        This method assumes perfect cobras (i.e., it doesn't use the motor
        maps).
        
        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        cobraCenters = self.bench.cobras.centers
        L1 = self.bench.cobras.L1
        L2 = self.bench.cobras.L2
        posThtMovement = self.movementDirections[0]
        thtEarly = self.movementStrategies[0]
        phiEarly = self.movementStrategies[1]
        
        # Get the cobra rotation angles for the starting and the final fiber
        # positions
        (startTht, startPhi) = self.bench.cobras.calculateRotationAngles(self.startFiberPositions)
        (finalTht, finalPhi) = self.bench.cobras.calculateRotationAngles(self.finalFiberPositions)
        
        # Calculate the required theta and phi delta offsets
        deltaTht = -np.mod(startTht - finalTht, 2 * np.pi)
        deltaTht[posThtMovement] = np.mod(finalTht[posThtMovement] - startTht[posThtMovement], 2 * np.pi)
        deltaPhi = finalPhi - startPhi
        
        # Reassign the final theta positions to be sure that they are
        # consistent with the deltaTht values
        finalTht = startTht + deltaTht
        
        # Calculate the theta and phi bin step widths
        binTht = np.sign(deltaTht) * TRAJECTORY_BIN_WIDTH
        binPhi = np.sign(deltaPhi) * TRAJECTORY_BIN_WIDTH
        
        # Calculate theta and phi angle values along the trajectories
        tht = np.empty((nCobras, self.nSteps))
        phi = np.empty((nCobras, self.nSteps))
        tht[:] = finalTht[:, np.newaxis]
        phi[:] = finalPhi[:, np.newaxis]
        
        for c in range(nCobras):
            # Jump to the next cobra if the two deltas are zero
            if deltaTht[c] == 0 and deltaPhi[c] == 0:
                continue
            
            # Get the theta and phi moves from the starting to the final
            # position
            movesTht = np.arange(startTht[c], finalTht[c], binTht[c]) if deltaTht[c] != 0 else np.array([])
            movesPhi = np.arange(startPhi[c], finalPhi[c], binPhi[c]) if deltaPhi[c] != 0 else np.array([])
            nMovesTht = len(movesTht)
            nMovesPhi = len(movesPhi)
            
            # Fill the rotation angles according to the movement strategies
            if thtEarly[c]:
                if phiEarly[c]:
                    # Early-early movement
                    tht[c, :nMovesTht] = movesTht
                    phi[c, :nMovesPhi] = movesPhi
                else:
                    # Early-late movement
                    tht[c, :nMovesTht] = movesTht
                    phi[c, :-nMovesPhi - 1] = startPhi[c]
                    phi[c, -nMovesPhi - 1:-1] = movesPhi
            else:
                if phiEarly[c]:
                    # Late-early movement. This should never happen!!
                    tht[c, :-nMovesTht - 1] = startTht[c]
                    tht[c, -nMovesTht - 1:-1] = movesTht
                    phi[c, :nMovesPhi] = movesPhi
                else:
                    # Late-late movement
                    if nMovesTht > nMovesPhi:
                        # If we have more theta moves, do the extra moves
                        # early, because the other cobras are still close to
                        # the center and the collision probability is smaller
                        extraMoves = nMovesTht - nMovesPhi
                        tht[c, :extraMoves] = movesTht[:extraMoves]
                        
                        # Keep the theta position before phi and theta start to
                        # move together
                        tht[c, extraMoves:-nMovesPhi - 1] = movesTht[extraMoves - 1]
                        
                        # Execute the rest of the theta movement together with
                        # the phi movement
                        tht[c, -nMovesPhi - 1:-1] = movesTht[extraMoves:]
                        phi[c, :-nMovesPhi - 1] = startPhi[c]
                        phi[c, -nMovesPhi - 1:-1] = movesPhi
                    else:
                        tht[c, :-nMovesTht - 1] = startTht[c]
                        tht[c, -nMovesTht - 1:-1] = movesTht
                        phi[c, :-nMovesPhi - 1] = startPhi[c]
                        phi[c, -nMovesPhi - 1:-1] = movesPhi
        
        # Calculate the elbow and fiber positions along the trajectory
        self.elbowPositions = cobraCenters[:, np.newaxis] + L1[:, np.newaxis] * np.exp(1j * tht)
        self.fiberPositions = self.elbowPositions + L2[:, np.newaxis] * np.exp(1j * (tht + phi))
    
    
    def calculateCobraAssociationCollisions(self, associationIndices=None):
        """Calculates which cobra associations are involved in a collision for
        each step in the trajectory.
        
        Parameters
        ----------
        associationIndices: object, optional
            A numpy array with the cobra associations indices to use. If it is
            set to None, all the cobra associations will be used. Default is
            None.
        
        Returns
        -------
        object
            A boolean numpy array indicating which cobra associations are
            involved in a collision for each step in the trajectory.
        
        """
        # Extract some useful information
        cobraAssociations = self.bench.cobraAssociations
        linkRadius = self.bench.cobras.linkRadius
        
        # Select a subset of the cobra associations if necessary
        if associationIndices is not None:
            cobraAssociations = cobraAssociations[:, associationIndices]
        
        # Calculate the distances between the cobras links for each step in the
        # trajectory
        startPoints1 = self.fiberPositions[cobraAssociations[0]].ravel()
        endPoints1 = self.elbowPositions[cobraAssociations[0]].ravel()
        startPoints2 = self.fiberPositions[cobraAssociations[1]].ravel()
        endPoints2 = self.elbowPositions[cobraAssociations[1]].ravel()
        distances = self.bench.distancesBetweenLineSegments(startPoints1, endPoints1, startPoints2, endPoints2)
        
        # Reshape the distances array
        distances = distances.reshape((len(cobraAssociations[0]), self.nSteps))
        
        # Return the cobra association collisions along the trajectory
        return distances < (linkRadius[cobraAssociations[0], np.newaxis] + linkRadius[cobraAssociations[1], np.newaxis])
    
    
    def addToFigure(self, colors=np.array([0.4, 0.4, 0.4, 1.0]), indices=None, paintFootprints=False, footprintColors=np.array([0.0, 0.0, 1.0, 0.05])):
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
        elbowPositions = self.elbowPositions
        fiberPositions = self.fiberPositions
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
        plotUtils.addTrajectories(np.vstack((elbowPositions, fiberPositions)), color=np.vstack((colors, colors)), linewidth=1)
        
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
            plotUtils.addThickLines(fiberPositions, elbowPositions, thiknesses, facecolor=footprintColors)
    
    
    def addAnimationToFigure(self, colors=np.array([0.4, 0.4, 0.4, 1.0]), linkColors=np.array([0.0, 0.0, 1.0, 0.5]), indices=None, fileName=None):
        """Animates the cobra trajectories on top of an existing figure.
        
        Parameters
        ----------
        colors: object, optional
            The trajectory colors. Default is dark grey.
        linkColors: object, optional
            The cobra link colors. Default is light blue.
        indices: object, optional
            A numpy array with the cobra trajectory indices to animate. If it
            is set to None, all the trajectories will be animated. Default is
            None.
        fileName: object, optional
            The file name path where a video of the animation should be saved.
            If it is set to None, no video will be saved. Default is None.
        
        """
        # Extract some useful information
        elbowPositions = self.elbowPositions
        fiberPositions = self.fiberPositions
        cobraCenters = self.bench.cobras.centers
        linkRadius = self.bench.cobras.linkRadius
        
        # Select a subset of the trajectories if necessary
        if indices is not None:
            elbowPositions = elbowPositions[indices]
            fiberPositions = fiberPositions[indices]
            cobraCenters = cobraCenters[indices]
            linkRadius = linkRadius[indices]
            
            if colors.ndim >= 2:
                colors = colors[indices]
            
            if linkColors.ndim >= 2:
                linkColors = linkColors[indices]
        
        # Define the animation update function
        lineCollection = None
        thickLineCollection = None
        trajectoryCollection = None
        
        def update(frame):
            # The function should be able to modify these variables
            nonlocal lineCollection
            nonlocal thickLineCollection
            nonlocal trajectoryCollection
            
            # Remove the cobras line collections painted in the previous step
            if lineCollection is not None:
                plotUtils.plt.gca().collections.remove(lineCollection)
                plotUtils.plt.gca().collections.remove(thickLineCollection)
                plotUtils.plt.gca().collections.remove(trajectoryCollection)
            
            # Draw the cobras using a combination of thin and thick lines
            lineCollection = plotUtils.addLines(cobraCenters, elbowPositions[:, frame], edgecolor=linkColors, linewidths=2)
            thickLineCollection = plotUtils.addThickLines(elbowPositions[:, frame], fiberPositions[:, frame], linkRadius, facecolors=linkColors)
            
            # Draw also their line trajectories
            combinedTrajectories = np.vstack((elbowPositions[:, :frame + 1], fiberPositions[:, :frame + 1]))
            trajectoryCollection = plotUtils.addTrajectories(combinedTrajectories, color=np.vstack((colors, colors)), linewidth=1)
        
        # Add the animation to the current figure
        plotUtils.addAnimation(update, elbowPositions.shape[1], fileName=fileName)
