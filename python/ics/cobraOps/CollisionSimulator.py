"""

Collision simulator class.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

from ics.cobraOps.TrajectoryGroup import TrajectoryGroup
from ics.cobraOps.cobraConstants import (TRAJECTORY_BIN_WIDTH,
                                         TRAJECTORY_MAX_STEPS,
                                         NULL_TARGET_ID)


class CollisionSimulator():
    """
    
    Class used to simulate a PFS observation.
    
    """
    
    def __init__(self, bench, targets):
        """Constructs a new collision simulator instance.
        
        Parameters
        ----------
        bench: object
            The PFI bench instance.
        targets: object
            The target group instance.
        
        Returns
        -------
        object
            The collision simulator instance.
        
        """
        # Save the bench and target group instances
        self.bench = bench
        self.targets = targets
        
        # Define some internal arrays that will filled by the run method
        self.assignedCobras = None
        self.finalFiberPositions = None
        self.nStepsP = None
        self.nStepsN = None
        self.positiveThtMovement = None
        self.positivePhiMovement = None
        self.movementStrategies = None
        self.trajectories = None
        self.trajectoryCollisions = None
        self.associationCollisions = None
        self.collisionDetected = None
        self.endPointCollisionDetected = None
    
    
    def run(self, solveCollisions=True):
        """Runs the collisions simulator.
        
        Parameters
        ----------
        solveCollisions: bool, optional
            If True, the simulator will try to solve trajectory collisions,
            changing the theta movement direction of the affected cobras.
            Default is True.
        
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
        self.detectCollisions()
        
        # Solve trajectory collisions if requested
        if solveCollisions:
            self.solveTrajectoryCollisions(selectLowerIndices=True)
            self.solveTrajectoryCollisions(selectLowerIndices=False)
            self.solveTrajectoryCollisions(selectLowerIndices=True)
    
    
    def solveTrajectoryCollisions(self, selectLowerIndices=True):
        """Solves trajectory collisions changing the cobras theta movement
        directions.
        
        Parameters
        ----------
        selectLowerIndices: bool, optional
            If True, the cobras selected to be changed will be the ones with
            the lower index values. Default is True.
        
        """
        # Get the cobras whose theta movement direction should be swapped
        cobrasToSwap = self.getCobrasToSwap(selectLowerIndices=selectLowerIndices)
        
        # Make sure that there is at least one cobra to change
        if len(cobrasToSwap) > 0:
            # Swapp the cobras theta movement direction
            self.swapThetaMovementDirection(cobrasToSwap)
            
            # Define the theta and phi movement strategies
            self.defineMovementStrategies()
            
            # Calculate the cobra trajectories
            self.calculateTrajectories()
            
            # Detect cobra collisions during the trajectory
            self.detectCollisions()
    
    
    def calculateFinalFiberPositions(self):
        """Calculates the cobras final fiber positions.
        
        """
        # Check which cobras are assigned to a target
        self.assignedCobras = self.targets.notNull.copy()
        
        # Set the final fiber positions to their associated target positions,
        # leaving unassigned cobras at their home positive positions
        self.finalFiberPositions = self.bench.cobras.home0.copy()
        self.finalFiberPositions[self.assignedCobras] = self.targets.positions[self.assignedCobras]
    
    
    def defineMovementDirections(self):
        """Defines the cobras theta and phi movement directions.
        
        """
        # Get the cobra rotation angles for the starting home positions and the
        # final fiber positions
        (startThtP, startPhiP) = self.bench.cobras.calculateRotationAngles(self.bench.cobras.home0)
        (startThtN, startPhiN) = self.bench.cobras.calculateRotationAngles(self.bench.cobras.home1)
        (finalTht, finalPhi) = self.bench.cobras.calculateRotationAngles(self.finalFiberPositions)
        
        # Calculate the required theta and phi delta offsets to move from the
        # positive and negative starting positions to the final positions
        deltaThtP = np.mod(finalTht - startThtP, 2 * np.pi)
        deltaThtN = -np.mod(startThtN - finalTht, 2 * np.pi)
        deltaPhiP = finalPhi - startPhiP
        deltaPhiN = finalPhi - startPhiN
        
        # Calculate the number of steps in the positive and the negative
        # directions
        self.nStepsP = np.ceil(np.max((np.abs(deltaThtP), np.abs(deltaPhiP)), axis=0) / TRAJECTORY_BIN_WIDTH).astype("int") + 1
        self.nStepsN = np.ceil(np.max((np.abs(deltaThtN), np.abs(deltaPhiN)), axis=0) / TRAJECTORY_BIN_WIDTH).astype("int") + 1
        
        # Make sure that at least one of the movements requires less steps than
        # the maximum number of steps allowed
        if np.any(np.min((self.nStepsP, self.nStepsN), axis=0) > TRAJECTORY_MAX_STEPS):
            raise Exception("TRAJECTORY_MAX_STEPS value should be set to a higher value.")
        
        # Decide if the cobras should follow a positive theta movement
        # direction:
        #
        #  - The positive movement should require less steps than the negative
        #    movement.
        #  - The cobra fibers should not go too far in phi (it is easier to
        #    have collisions when moving in the positive direction).
        self.positiveThtMovement = np.logical_and(self.nStepsP < self.nStepsN, finalPhi < -0.3 * np.pi)
        
        # Select the positive movement if the negative movement would require
        # too many steps
        self.positiveThtMovement[self.nStepsN > TRAJECTORY_MAX_STEPS] = True
        
        # Calculate the phi movement direction
        self.positivePhiMovement = deltaPhiN > 0
        self.positivePhiMovement[self.positiveThtMovement] = deltaPhiP[self.positiveThtMovement] > 0
    
    
    def defineMovementStrategies(self):
        """Defines the movement strategies that the cobras should follow.
        
        There are only three options that really make sense:
        
            - Early theta and early phi movements
            - Early theta and late phi movements
            - Late theta and late phi movements
        
        """
        # By default always do late-late movements
        thtEarly = np.full(self.bench.cobras.nCobras, False)
        phiEarly = thtEarly.copy()
        
        # Do early theta movements when moving in the theta positive direction
        thtEarly[self.positiveThtMovement] = True
        
        # If the cobra is moving towards the center, do the theta and phi
        # movements as early as possible
        thtEarly[self.positivePhiMovement == False] = True
        phiEarly[self.positivePhiMovement == False] = True
        
        # Save the results in a single array
        self.movementStrategies = np.vstack((thtEarly, phiEarly))
    
    
    def calculateTrajectories(self):
        """Obtains the cobra trajectories.
        
        """
        self.trajectories = TrajectoryGroup(self.bench, self.finalFiberPositions, self.positiveThtMovement, self.movementStrategies)
    
    
    def detectCollisions(self):
        """Detects collisions in the cobra trajectories.
        
        """
        # Detect cobra to nearby cobra collisions at the same trajectory step
        fiberPositions = self.trajectories.fiberPositions
        self.trajectoryCollisions = np.full(fiberPositions.shape, False)
        self.associationCollisions = np.full(self.bench.nearestNeighbors.shape[1], False)
        
        for i in range(self.trajectories.nSteps):
            # Calculate the nearest neighbors cobra association collisions
            collisions = self.bench.calculateNearestNeighborsCollisions(fiberPositions[:, i])
            
            # Get the indices of the cobras involved in the collisions
            problematicCobras = self.bench.nearestNeighbors[0, collisions]
            nearbyProblematicCobras = self.bench.nearestNeighbors[1, collisions]
            
            # Update the collisions arrays
            self.trajectoryCollisions[problematicCobras, i] = True
            self.trajectoryCollisions[nearbyProblematicCobras, i] = True
            self.associationCollisions[collisions] = True
        
        # Check which cobras are involved in collisions
        self.collisionDetected = np.any(self.trajectoryCollisions, axis=1)
        self.endPointCollisionDetected = self.trajectoryCollisions[:, -1]
    
    
    def getCobrasToSwap(self, selectLowerIndices=True):
        """Returns the indices of the cobras that should be swapped.
        
        Parameters
        ----------
        selectLowerIndices: bool, optional
            If True, the cobras selected to be swapped will be the ones with
            the lower index values. Default is True.
        
        Returns
        -------
        object
            A numpy array with the indices of the cobras that should be swapped.
        
        """
        # Get the indices of the cobra associations where we have a collision
        # at some point in the trajectory
        problematicCobras = self.bench.nearestNeighbors[0, self.associationCollisions]
        nearbyProblematicCobras = self.bench.nearestNeighbors[1, self.associationCollisions]
        
        # Check which cobra collision association corresponds to an end collision
        endCollisions = np.logical_and(self.trajectoryCollisions[problematicCobras, -1], self.trajectoryCollisions[nearbyProblematicCobras, -1])
        
        # Get the cobra association with a collision that is not an end collision
        problematicCobras = problematicCobras[endCollisions == False]
        nearbyProblematicCobras = nearbyProblematicCobras[endCollisions == False]
        
        # Get the indices of the cobras that should be swapped
        cobrasToSwap = problematicCobras if selectLowerIndices else nearbyProblematicCobras
        
        return np.unique(cobrasToSwap)
    
    
    def swapThetaMovementDirection(self, cobrasToSwap):
        """Swaps the theta movement direction for the given cobras.
        
        Parameters
        ----------
        cobrasToSwap: object
            A numpy array with the indices of the cobras to swap.
        
        """
        # Swap the given cobras movement direction
        self.positiveThtMovement[cobrasToSwap] = np.logical_not(self.positiveThtMovement[cobrasToSwap])
    
        # Make sure that the new selected movement doesn't require too many steps
        self.positiveThtMovement[self.nStepsP > TRAJECTORY_MAX_STEPS] = False
        self.positiveThtMovement[self.nStepsN > TRAJECTORY_MAX_STEPS] = True
