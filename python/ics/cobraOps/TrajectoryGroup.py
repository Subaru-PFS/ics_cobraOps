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
from ics.cobraOps.cobraConstants import TRAJECTORY_BIN_WIDTH


class TrajectoryGroup(AttributePrinter):
    """
    
    Class describing the properties of a group of cobra trajectories.
    
    """
    
    def __init__(self, bench, finalFiberPositions, thtMovementDirection, movementStrategies):
        """Constructs a new trajectory group instance.
        
        Parameters
        ----------
        bench: object
            The PFI bench instance.
        finalFiberPositions: object
            A complex numpy array with the trajectory final fiber positions for
            each cobra.
        thtMovementDirection: object
            A boolean numpy array with the theta movement direction to use.
            True values indicate that the cobras should move in the positive
            theta direction, while False values indicate that the movement
            should be in the negative theta direction.
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
        self.thtMovementDirection = thtMovementDirection.copy()
        self.movementStrategies = movementStrategies.copy()
        
        # Calculate the cobra trajectories
        self.calculateCobraTrajectories()
    
    
    def calculateCobraTrajectories(self):
        """Calculates the cobra trajectories starting from their home
        positions.
        
        This method assumes perfect cobras (i.e., it doesn't use the motor
        maps).
        
        """
        # Extract some useful information
        nCobras = self.bench.cobras.nCobras
        cobraCenters = self.bench.cobras.centers
        L1 = self.bench.cobras.L1
        L2 = self.bench.cobras.L2
        home0 = self.bench.cobras.home0
        home1 = self.bench.cobras.home1
        earlyThtMovement = self.movementStrategies[0]
        earlyPhiMovement = self.movementStrategies[1]
        
        # Set the start positions according to the specified movement direction
        self.startFiberPositions = home1.copy()
        self.startFiberPositions[self.thtMovementDirection] = home0[self.thtMovementDirection]
        
        # Get the cobra rotation angles from the starting to the final
        # positions
        (startTht, startPhi) = self.bench.cobras.calculateRotationAngles(self.startFiberPositions)
        (finalTht, finalPhi) = self.bench.cobras.calculateRotationAngles(self.finalFiberPositions)
        
        # Calculate the required theta and phi delta offsets
        deltaTht = -np.mod(startTht - finalTht, 2 * np.pi)
        deltaTht[self.thtMovementDirection] = np.mod(finalTht[self.thtMovementDirection] - startTht[self.thtMovementDirection], 2 * np.pi)
        deltaPhi = finalPhi - startPhi
        
        # Reassign the final theta positions to be sure that they are
        # consistent with the deltaTht values
        finalTht = startTht + deltaTht
        
        # Calculate the theta and phi bin step widths
        binTht = np.full(nCobras, -TRAJECTORY_BIN_WIDTH)
        binTht[deltaTht > 0] = TRAJECTORY_BIN_WIDTH
        binPhi = np.full(nCobras, -TRAJECTORY_BIN_WIDTH)
        binPhi[deltaPhi > 0] = TRAJECTORY_BIN_WIDTH
        
        # Calculate the number of trajectory steps
        self.nSteps = int(np.ceil(np.max((deltaTht / binTht, deltaPhi / binPhi))) + 1)
        
        # Calculate theta and phi trajectories
        trajectoriesTht = np.empty((nCobras, self.nSteps))
        trajectoriesPhi = np.empty((nCobras, self.nSteps))
        trajectoriesTht[:] = finalTht[:, np.newaxis]
        trajectoriesPhi[:] = finalPhi[:, np.newaxis]
        
        for c in range(nCobras):
            # Jump to the next cobra if the two deltas are zero (unassigned
            # cobras)
            if deltaTht[c] == 0 and deltaPhi[c] == 0:
                continue
            
            # Get the theta and phi moves from the starting position to the
            # final positions
            movesTht = np.arange(startTht[c], finalTht[c], binTht[c])
            movesPhi = np.arange(startPhi[c], finalPhi[c], binPhi[c])
            nMovesTht = len(movesTht)
            nMovesPhi = len(movesPhi)
            
            # Fill the trajectories according to the movement strategies
            if earlyThtMovement[c]:
                if earlyPhiMovement[c]:
                    # Early-early movement
                    trajectoriesTht[c, :nMovesTht] = movesTht
                    trajectoriesPhi[c, :nMovesPhi] = movesPhi
                else:
                    # Early-late movement
                    trajectoriesTht[c, :nMovesTht] = movesTht
                    trajectoriesPhi[c, :-nMovesPhi - 1] = startPhi[c]
                    trajectoriesPhi[c, -nMovesPhi - 1:-1] = movesPhi
            else:
                if earlyPhiMovement[c]:
                    # Late-early movement (this should never happen!!)
                    trajectoriesTht[c, :-nMovesTht - 1] = startTht[c]
                    trajectoriesTht[c, -nMovesTht - 1:-1] = movesTht
                    trajectoriesPhi[c, :nMovesPhi] = movesPhi
                else:
                    # Late-late movement
                    if nMovesTht > nMovesPhi:
                        # If we have more moves in theta, do the extra theta
                        # moves early, because the other cobras are still close
                        # to the center and that decreases the collision
                        # probability
                        extraMoves = nMovesTht - nMovesPhi
                        trajectoriesTht[c, :extraMoves] = movesTht[:extraMoves]
                        
                        # Keep the theta position before phi and theta start to
                        # move together
                        trajectoriesTht[c, extraMoves:-nMovesPhi - 1] = movesTht[extraMoves - 1]
                        
                        # Execute the rest of the theta movement together with
                        # the phi movement
                        trajectoriesTht[c, -nMovesPhi - 1:-1] = movesTht[extraMoves:]
                        trajectoriesPhi[c, :-nMovesPhi - 1] = startPhi[c]
                        trajectoriesPhi[c, -nMovesPhi - 1:-1] = movesPhi
                    else:
                        trajectoriesTht[c, :-nMovesTht - 1] = startTht[c]
                        trajectoriesTht[c, -nMovesTht - 1:-1] = movesTht
                        trajectoriesPhi[c, :-nMovesPhi - 1] = startPhi[c]
                        trajectoriesPhi[c, -nMovesPhi - 1:-1] = movesPhi
        
        # Calculate the elbow and fiber positions along the trajectory
        self.elbowPositions = cobraCenters[:, np.newaxis] + L1[:, np.newaxis] * np.exp(1j * trajectoriesTht)
        self.fiberPositions = self.elbowPositions + L2[:, np.newaxis] * np.exp(1j * (trajectoriesTht + trajectoriesPhi))
    
    
    def addToFigure(self, colors=[0.4, 0.4, 0.4, 1.0], indices=None, paintFootprints=False, footprintColors=[0.0, 0.0, 1.0, 0.05]):
        """Draws the cobra trajectories on top off an existing figure.
        
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
