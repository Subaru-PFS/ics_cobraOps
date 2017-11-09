"""

Some utility methods related with the cobra trajectories.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

import ics.cobraOps.benchUtils as benchUtils
import ics.cobraOps.plotUtils as plotUtils
import ics.cobraOps.targetUtils as targetUtils


MAX_TRAJECTORY_STEPS = 80
"""The maximum number of steps that a cobra trajectory is allowed to have."""


def defineThetaMovementDirection(finalPositions, bench):
    """Defines the cobras theta movement directions starting from the bench
    home positions.

    Parameters
    ----------
    finalPositions: object
        A complex numpy array with the cobras fiber final positions.
    bench: object
        The bench geometry to use.

    Returns
    -------
    tuple
        A python tuple with the cobras theta movement direction information
        (True values indicate that the cobras should move in the positive theta
        direction, while False values indicate that the movement should be in
        the negative direction), and two arrays with the total number of steps
        required for the positive and the negative movements.

    """
    # Extract some useful information from the bench geometry
    cobraCenters = bench["center"]
    L1 = bench["L1"]
    L2 = bench["L2"]
    home0 = bench["home0"]
    home1 = bench["home1"]
    binWidth = bench["binWidth"]

    # Get the cobra rotation angles for the starting home positions and the
    # final positions
    (startThtP, startPhiP) = benchUtils.getCobraRotationAngles(home0 - cobraCenters, L1, L2)
    (startThtN, startPhiN) = benchUtils.getCobraRotationAngles(home1 - cobraCenters, L1, L2)
    (finalTht, finalPhi) = benchUtils.getCobraRotationAngles(finalPositions - cobraCenters, L1, L2)
    
    # Calculate the required theta and phi delta offsets to move from the
    # positive and negative starting positions to the final positions
    deltaThtP = np.mod(finalTht - startThtP, 2 * np.pi)
    deltaThtN = -np.mod(startThtN - finalTht, 2 * np.pi)
    deltaPhiP = finalPhi - startPhiP
    deltaPhiN = finalPhi - startPhiN
    
    # Calculate the number of steps in the positive and the negative directions
    nStepsP = np.ceil(np.max((np.abs(deltaThtP / binWidth), np.abs(deltaPhiP / binWidth)), axis=0)).astype("int") + 1
    nStepsN = np.ceil(np.max((np.abs(deltaThtN / binWidth), np.abs(deltaPhiN / binWidth)), axis=0)).astype("int") + 1

    # Make sure that at least one of the movements requires less steps than the
    # maximum number of steps allowed
    if np.any(np.min((nStepsP, nStepsN), axis=0) > MAX_TRAJECTORY_STEPS):
        raise Exception("MAX_TRAJECTORY_STEPS value should be set to a higher value.")

    # Decide if the cobras should follow a positive theta movement direction:
    #  - The positive movement should require less steps than the negative
    #    movement.
    #  - The cobra fibers should not go too far in phi (it is easier to have
    #    collisions when moving in the positive direction).
    positiveMovement = np.logical_and(nStepsP < nStepsN, finalPhi < -0.3 * np.pi)   

    # Select the positive movement if the negative movement would require too
    # many steps
    positiveMovement[nStepsN > MAX_TRAJECTORY_STEPS] = True
 
    return (positiveMovement, nStepsP, nStepsN)


def calculateTrajectories(finalPositions, positiveMovement, bench):
    """Calculates the cobra trajectories starting from their home positions.
    
    This method assumes perfect cobras (i.e., it doesn't use the motor maps).

    Parameters
    ----------
    finalPositions: object
        A complex numpy array with the cobras fiber final positions.
    positiveMovement: object
        A boolean numpy array with the cobras theta movement direction. True
        values indicate that the cobras should move in the positive theta
        direction, while False values indicate that the movement should be in
        the negative direction.
    bench: object
        The bench geometry to use.

    Returns
    -------
    tuple
        A python tuple containing two numpy complex arrays with the cobras
        elbow and fiber trajectories.

    """
    # Extract some useful information from the bench geometry
    cobraCenters = bench["center"]
    L1 = bench["L1"]
    L2 = bench["L2"]
    home0 = bench["home0"]
    home1 = bench["home1"]
    binWidth = bench["binWidth"]

    # Set the start positions according to the specified movement direction
    startPositions = home1.copy()
    startPositions[positiveMovement] = home0[positiveMovement]
    
    # Get the cobra rotation angles for the starting to the final positions
    (startTht, startPhi) = benchUtils.getCobraRotationAngles(startPositions - cobraCenters, L1, L2)
    (finalTht, finalPhi) = benchUtils.getCobraRotationAngles(finalPositions - cobraCenters, L1, L2)
    
    # Calculate the required theta and phi delta offsets
    deltaTht = -np.mod(startTht - finalTht, 2 * np.pi)
    deltaTht[positiveMovement] = np.mod(finalTht[positiveMovement] - startTht[positiveMovement], 2 * np.pi)
    deltaPhi = finalPhi - startPhi
    
    # Reassign the final theta positions to be sure that they are consistent
    # with the deltaTht values
    finalTht = startTht + deltaTht
    
    # Calculate the theta and phi bin step widths
    nCobras = len(cobraCenters)
    binTht = np.full(nCobras, -binWidth)
    binTht[deltaTht > 0] = binWidth
    binPhi = np.full(nCobras, -binWidth)
    binPhi[deltaPhi > 0] = binWidth
    
    # Calculate the maximum number of steps that a trajectory could have
    maxSteps = int(np.ceil(np.max((deltaTht / binTht, deltaPhi / binPhi))) + 1) 
        
    # Calculate theta and phi trajectories
    trajectoriesTht = np.empty((nCobras, maxSteps))
    trajectoriesTht[:] = finalTht[:, np.newaxis]
    trajectoriesPhi = np.empty((nCobras, maxSteps))
    trajectoriesPhi[:] = finalPhi[:, np.newaxis]

    for c in range(nCobras):
        # Jump to the next cobra if the two deltas are zero (unassigned cobras)
        if deltaTht[c] == 0 and deltaPhi[c] == 0:
            continue
        
        # Get the theta and phi moves from the starting position to the final
        # positions
        movesTht = np.arange(startTht[c], finalTht[c], binTht[c])
        movesPhi = np.arange(startPhi[c], finalPhi[c], binPhi[c])
        nMovesTht = len(movesTht)
        nMovesPhi = len(movesPhi)
   
        # Check if the phi move is going towards the center
        if binPhi[c] < 0:
            # Moving towards the center: do the theta and phi movements early,
            # because the other cobras are still close to the center and that
            # decreases the collision probability
            trajectoriesTht[c, :nMovesTht] = movesTht
            trajectoriesPhi[c, :nMovesPhi] = movesPhi
        else:
            # Moving outwards: we should differentiate between positive and
            # negative theta movements
            if positiveMovement[c]:
                # Always make early theta movements, because the other cobras
                # are still close to the center and that decreases the
                # collision probability
                trajectoriesTht[c, :nMovesTht] = movesTht
            else:
                # Negative movement: check which angle requires more moves to
                # complete
                if nMovesTht > nMovesPhi:
                    # If we have more moves in theta, do the extra theta moves
                    # early, because the other cobras are still close to the
                    # center and that decreases the collision probability
                    extraMoves = nMovesTht - nMovesPhi
                    trajectoriesTht[c, :extraMoves] = movesTht[:extraMoves]

                    # Keep the theta position before phi and theta start to
                    # move together
                    trajectoriesTht[c, extraMoves:-nMovesPhi - 1] = movesTht[extraMoves - 1]
                
                    # Execute the rest of the theta movement together with the
                    # phi movement: as late as possible to decrease the
                    # collision probability
                    trajectoriesTht[c, -nMovesPhi - 1:-1] = movesTht[-nMovesPhi:]
                else:
                    # Move theta as late as possible
                    trajectoriesTht[c, :-nMovesTht - 1] = startTht[c]
                    trajectoriesTht[c, -nMovesTht - 1:-1] = movesTht

            # Do the phi movements as late as possible to decrease the
            # collision probability
            trajectoriesPhi[c, :-nMovesPhi - 1] = startPhi[c]
            trajectoriesPhi[c, -nMovesPhi - 1:-1] = movesPhi

    # Calculate the elbow and fiber trajectories
    elbowTrajectories = cobraCenters[:, np.newaxis] + L1[:, np.newaxis] * np.exp(1j * trajectoriesTht)
    fiberTrajectories = elbowTrajectories + L2[:, np.newaxis] * np.exp(1j * (trajectoriesTht + trajectoriesPhi))

    return (elbowTrajectories, fiberTrajectories)


def detectCollisions(elbowTrajectories, fiberTrajectories, bench):
    """Detects collisions in the cobra trajectories.
    
    Parameters
    ----------
    elbowTrajectories: object
        A complex numpy array with the cobra elbow trajectories.
    fiberTrajectories: object
        A complex numpy array with the cobra fiber trajectories.
    bench: object
        The bench geometry to use.

    Returns
    -------
    tuple
        A python tuple containing two numpy array with the indices of the
        cobras involved in a collision and a boolean numpy array indicating the
        points in the trajectory where a collision is detected.

    """
    # Extract some useful information from the bench geometry
    minDist = bench["minDist"]

    # Get the bench precalculated nearest neighbors information 
    cobras = bench["NN"]["row"]
    nearbyCobras = bench["NN"]["col"]

    # We only need to test half of the cobra associations
    uniqueAssociations = cobras < nearbyCobras
    cobras = cobras[uniqueAssociations]
    nearbyCobras = nearbyCobras[uniqueAssociations]
    
    # Detect cobra to nearby cobra collisions at the same trajectory time step    
    (nCobras, nSteps) = elbowTrajectories.shape
    cobraCollisions = np.full(len(cobras), False)
    trajectoryCollisions = np.full((nCobras, nSteps), False)
    
    for i in range(nSteps): 
        # Get the cobra elbow and fiber positions at the given time
        elbowPositions = elbowTrajectories[:, i]
        fiberPositions = fiberTrajectories[:, i]

        # Calculate the distances between the cobras links and the nearby cobras links
        startPoints1 = elbowPositions[cobras]
        endPoints1 = fiberPositions[cobras]
        startPoints2 = elbowPositions[nearbyCobras]
        endPoints2 = fiberPositions[nearbyCobras]
        distances = targetUtils.distanceBetweenLineSegments(startPoints1, endPoints1, startPoints2, endPoints2)

        # Get the collisions
        collisions = distances < (minDist[cobras] + minDist[nearbyCobras]) / 2
        
        # Update the cobra collisions and the trajectory collisions arrays
        cobraCollisions[collisions] = True
        trajectoryCollisions[cobras[collisions], i] = True
        trajectoryCollisions[nearbyCobras[collisions], i] = True

    # Obtain the indices of the problematic cobras and their associations and
    # don't forget to include the other half of the associations
    problematicCobras = np.concatenate((cobras[cobraCollisions], nearbyCobras[cobraCollisions]))
    nearbyProblematicCobras = np.concatenate((nearbyCobras[cobraCollisions], cobras[cobraCollisions]))
   
    return (problematicCobras, nearbyProblematicCobras, trajectoryCollisions)


def getCobrasToSwap(cobras, nearbyCobras, trajectoryCollisions, selectLowerIndices=True):
    """Returns the indices of the cobras that should be swapped.

    Parameters
    ----------
    cobras: object
        A numpy array with the indices of the cobras involved in a collision.
    nearbyCobras: object
        A numpy array with the indices of the secondary cobras involved in a
        collision.
    trajectoryCollisions: object
        A boolean numpy array indicating the points in the cobras trajectories
        where a collision has been detected.
    selectLowerIndices: bool, optional
        If True, the cobra selected to be swapped will be the one with the
        lower index value. Default is True.

    Returns
    -------
    object
        A numpy array with the indices of the cobras that should be swapped.

    """
    # Check which cobra collision association corresponds to an end collision
    endCollisions = np.logical_and(trajectoryCollisions[cobras, -1], trajectoryCollisions[nearbyCobras, -1])
    
    # Get the cobra collision association that should be swapped
    swapMovement = np.logical_and(endCollisions == False, (cobras < nearbyCobras) == selectLowerIndices)
    
    return np.unique(cobras[swapMovement])


def swapThetaMovementDirection(cobrasToSwap, movementDirection, nStepsP, nStepsN):
    """Swaps the theta movement direction for the given cobras.

    Parameters
    ----------
    cobrasToSwap: object
        A numpy array with the indices of the cobras to swap.
    movementDirection: object
        A boolean numpy array with the current cobras theta movement direction.
        True values indicate that the cobras should move in the positive theta
        direction, while False values indicate that the movement should be in
        the negative direction.
    nStepsP: object
        A numpy array with the total number of steps required for the positive
        movement.
    nStepsN: object
        A numpy array with the total number of steps required for the negative
        movement.

    Returns
    -------
    object
        A boolean numpy array with the update cobras theta movement direction.
        True values indicate that the cobras should move in the positive theta
        direction, while False values indicate that the movement should be in
        the negative direction.

    """
    # Swap the given cobras movement direction
    newMovementDirection = movementDirection.copy()
    newMovementDirection[cobrasToSwap] = np.logical_not(movementDirection[cobrasToSwap])
    
    # Make sure that the new selected movement doesn't require too many steps
    newMovementDirection[nStepsP > MAX_TRAJECTORY_STEPS] = False
    newMovementDirection[nStepsN > MAX_TRAJECTORY_STEPS] = True

    return newMovementDirection


def plotTrajectories(elbowTrajectories, fiberTrajectories, bench, paintFootprints=False, footprintColors=[0.0, 0.0, 1.0, 0.05]):
    """Plots cobras trajectories.

    Parameters
    ----------
    elbowTrajectories: object
        A complex numpy array with the elbow trajectory positions for each
        cobra.
    fiberTrajectories: object
        A complex numpy array with the fiber trajectory positions for each
        cobra.
    bench: object
        The bench geometry to use.
    paintFootprints: bool, optional
        If True, the cobra footprints will be painted. Default is False.
    footprintColors: object, optional
        The cobra footprints colors. Default is very light blue.
    
    """
    # Plot the elbow and fiber trajectories as continuous lines
    plotUtils.addTrajectories(np.vstack((elbowTrajectories, fiberTrajectories)), color="0.4", linewidth=1)
   
    # Paint the cobra trajectory footprints if necessary
    if paintFootprints:
        # Calculate the line thicknesses
        (nCobras, nSteps) = elbowTrajectories.shape
        thiknesses = np.empty((nCobras, nSteps)) 
        thiknesses[:] = 0.5 * bench["minDist"][:, np.newaxis]

        # Use only those elbow and fiber positions where the cobra is moving
        isMoving = np.empty((nCobras, nSteps), dtype="bool")
        isMoving[:, :-1] = (fiberTrajectories[:, 1:] - fiberTrajectories[:, :-1]) != 0
        isMoving[:, -1] = isMoving[:, -2]
        elbowPositions = elbowTrajectories[isMoving]
        fiberPositions = fiberTrajectories[isMoving]
        thiknesses = thiknesses[isMoving]
        
        # Update the colors if necessary
        if footprintColors.ndim > 2 and footprintColors.shape[:2] == isMoving.shape:
            # Set the colors for the moving positions
            footprintColors = footprintColors[isMoving]
            
            # Only use positions where the alpha color is not exactly zero
            visible = footprintColors[:, 3] != 0
            elbowPositions = elbowPositions[visible]
            fiberPositions = fiberPositions[visible]
            footprintColors = footprintColors[visible]
            thiknesses = thiknesses[visible]
            
        # Represent the trajectory footprint as a combination of thick lines
        plotUtils.addThickLines(fiberPositions, elbowPositions, thiknesses, facecolor=footprintColors)


def animateCobraTrajecty(centralCobraIndex, elbowTrajectories, fiberTrajectories, bench, cobraColors=[0.0, 0.0, 1.0, 0.5], fileName=False):
    """Animates the trajectory of a given cobra.

    Parameters
    ----------
    centraCobraIndex: int
        The index of the cobra that should be animated.
    elbowTrajectories: object
        A complex numpy array with the elbow trajectory positions for each
        cobra.
    fiberTrajectories: object
        A complex numpy array with the fiber trajectory positions for each
        cobra.
    bench: object
        The bench geometry to use.
    cobraColors: object, optional
        The cobra footprints colors. Default is very light blue.
    fileName: object, optional
        The video file name path. If it is set to None, no video will be saved.
        Default is None.
    
    """
    # Extract some useful information from the bench geometry
    cobraCenters = bench["center"]
    rMax = bench["rMax"]
    minDist = bench["minDist"]

    # Get the bench precalculated nearest neighbors information 
    cobras = bench["NN"]["row"]
    nearbyCobras = bench["NN"]["col"]

    # Get the indices of the cobras that are near the central cobra
    nearCentralCobra = nearbyCobras[cobras == centralCobraIndex]
 
    # Set the figure axes limits using the distance from those cobras
    radius = np.max(np.abs(cobraCenters[nearCentralCobra] - cobraCenters[centralCobraIndex]) + rMax[nearCentralCobra])
    limRange = radius * np.array([-1, 1]) 
    xLim = cobraCenters[centralCobraIndex].real + limRange
    yLim = cobraCenters[centralCobraIndex].imag + limRange
    plotUtils.setAxesLimits(xLim, yLim)

    # Check which cobras should be animated: only those that fall inside the figure
    toAnimate = np.full(len(cobraCenters), False)
    toAnimate[centralCobraIndex] = True
    toAnimate[nearCentralCobra] = True
    toAnimate[nearbyCobras[np.in1d(cobras, nearCentralCobra)]] = True   
    
    # Plot the cobras that should not be animated at their final position
    benchUtils.plotCobras(bench, fiberTrajectories[:, -1], cobraColors=cobraColors, cobraIndices=np.logical_not(toAnimate))

    # Limit the relevant arrays to those cobras
    cobraCenters = cobraCenters[toAnimate]
    elbowTrajectories = elbowTrajectories[toAnimate]
    fiberTrajectories = fiberTrajectories[toAnimate]
    minDist = minDist[toAnimate]

    if cobraColors.ndim == 2:
        cobraColors = cobraColors[toAnimate]

    # Define the update function
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
        
        # Paint the cobras that should be animated
        lineCollection = plotUtils.addLines(cobraCenters, elbowTrajectories[:, frame], edgecolor=cobraColors, linewidths=2)
        thickLineCollection = plotUtils.addThickLines(elbowTrajectories[:, frame], fiberTrajectories[:, frame], 0.5 * minDist, facecolors=cobraColors)  
    
        # Plot also their line trajectories
        trajectoryCollection = plotUtils.addTrajectories(np.vstack((elbowTrajectories[:, :frame + 1], fiberTrajectories[:, :frame + 1])), color="0.4", linewidth=1)

        # Log some animation information
        animationPercentage = int(100 * (frame + 1) / elbowTrajectories.shape[1])
        
        if animationPercentage % 20 == 0:
            print("Trajectory percentage: {}%".format(animationPercentage))

    # Add the animation to the current figure
    plotUtils.addAnimation(update, elbowTrajectories.shape[1], fileName=fileName)


def animateTrajectories(elbowTrajectories, fiberTrajectories, bench, cobraColors=[0.0, 0.0, 1.0, 0.5]):
    """Animates the cobra trajectories.

    Parameters
    ----------
    elbowTrajectories: object
        A complex numpy array with the elbow trajectory positions for each
        cobra.
    fiberTrajectories: object
        A complex numpy array with the fiber trajectory positions for each
        cobra.
    bench: object
        The bench geometry to use.
    cobraColors: object, optional
        The cobra footprints colors. Default is very light blue.
    
    """
    # Extract some useful information from the bench geometry
    cobraCenters = bench["center"]
    minDist = bench["minDist"]

    # Define the update function
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
        
        # Paint the cobras that should be animated
        lineCollection = plotUtils.addLines(cobraCenters, elbowTrajectories[:, frame], edgecolor=cobraColors, linewidths=2)
        thickLineCollection = plotUtils.addThickLines(elbowTrajectories[:, frame], fiberTrajectories[:, frame], 0.5 * minDist, facecolors=cobraColors)  
    
        # Plot also their line trajectories
        trajectoryCollection = plotUtils.addTrajectories(np.vstack((elbowTrajectories[:, :frame + 1], fiberTrajectories[:, :frame + 1])), color="0.4", linewidth=1)

        # Log some animation information
        animationPercentage = int(100 * (frame + 1) / elbowTrajectories.shape[1])
        
        if animationPercentage % 20 == 0:
            print("Trajectory percentage: {}%".format(animationPercentage))

    # Add the animation to the current figure
    plotUtils.addAnimation(update, elbowTrajectories.shape[1], fileName=None)


if __name__ == "__main__":
    # Import the necessary modules
    import time as time
    import ics.cobraOps.cobraUtils as cobraUtils

    # Define the target density to use
    targetDensity = 2.0

    # Get the cobras central positions for the full PFI
    start = time.time()
    cobraCenters = cobraUtils.getCobrasCenters("full")
    print("Number of cobras:", len(cobraCenters))

    # Define the bench geometry
    bench = benchUtils.defineBenchGeometry(cobraCenters, True, True)

    # Create a random sample of targets
    targetPositions = targetUtils.generateTargets(targetDensity, bench)
    print("Number of simulated targets:", len(targetPositions))

    # Assign the targets to the cobras and get the fiber positions
    (assignedTargets, fiberPositions) = targetUtils.assignTargets(targetPositions, bench)

    # Calculate the cobra trajectories
    (positiveMovement, nStepsP, nStepsN) = defineThetaMovementDirection(fiberPositions, bench)
    (elbowTrajectories, fiberTrajectories) = calculateTrajectories(fiberPositions, positiveMovement, bench)

    # Calculate the cobra trajectory collisions
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = detectCollisions(elbowTrajectories, fiberTrajectories, bench)
    print("Number of cobras affected by a collision (  I):", len(np.unique(problematicCobras)))

    # Swap the theta movement direction for some of the problematic cobras
    cobrasToSwap = getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, True)
    positiveMovement = swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
    (elbowTrajectories, fiberTrajectories) = calculateTrajectories(fiberPositions, positiveMovement, bench)
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = detectCollisions(elbowTrajectories, fiberTrajectories, bench)
    print("Number of cobras affected by a collision ( II):", len(np.unique(problematicCobras)))

    # Swap the theta movement direction for some of the problematic cobras
    cobrasToSwap = getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, False)
    positiveMovement = swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
    (elbowTrajectories, fiberTrajectories) = calculateTrajectories(fiberPositions, positiveMovement, bench)
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = detectCollisions(elbowTrajectories, fiberTrajectories, bench)
    print("Number of cobras affected by a collision (III):", len(np.unique(problematicCobras)))
 
    # Swap the theta movement direction for some of the problematic cobras
    cobrasToSwap = getCobrasToSwap(problematicCobras, nearbyProblematicCobras, trajectoryCollisions, True)
    positiveMovement = swapThetaMovementDirection(cobrasToSwap, positiveMovement, nStepsP, nStepsN)
    (elbowTrajectories, fiberTrajectories) = calculateTrajectories(fiberPositions, positiveMovement, bench)
    (problematicCobras, nearbyProblematicCobras, trajectoryCollisions) = detectCollisions(elbowTrajectories, fiberTrajectories, bench)
    print("Number of cobras affected by a collision ( IV):", len(np.unique(problematicCobras)))

    # Check how many of the collisions are not end point collisions
    endCollisions = np.logical_and(trajectoryCollisions[problematicCobras, -1], trajectoryCollisions[nearbyProblematicCobras, -1])
    print("Number of cobras unaffected by end collisions: ", len(np.unique(problematicCobras)) - len(np.unique(problematicCobras[endCollisions])))
    print("Total computation time (s):", time.time() - start)

    # Plot the cobra trajectories
    start = time.time()
    plotUtils.createNewFigure("Cobra trajectories", "x position", "y position")

    patrolAreaColors = np.full((len(cobraCenters), 4), [0.0, 0.0, 1.0, 0.15])
    patrolAreaColors[problematicCobras] = [1.0, 0.0, 0.0, 0.3]
    patrolAreaColors[problematicCobras[endCollisions]] = [0.0, 1.0, 0.0, 0.5]
    benchUtils.plotBenchGeometry(bench, patrolAreaColors)

    footprintColors = np.zeros((elbowTrajectories.shape[0], fiberTrajectories.shape[1], 4))
    footprintColors[problematicCobras, :] = [0.0, 0.0, 1.0, 0.05]
    plotTrajectories(elbowTrajectories, fiberTrajectories, bench, paintFootprints=True, footprintColors=footprintColors)

    cobraColors = np.full((len(cobraCenters), 4), [0.0, 0.0, 1.0, 0.5])
    cobraColors[assignedTargets == targetUtils.NULL_TARGET_INDEX] = [1.0, 0.0, 0.0, 0.25]
    benchUtils.plotCobras(bench, fiberPositions, cobraColors)

    targetColors = np.full((len(targetPositions), 4), [0.4, 0.4, 0.4, 1.0])
    targetColors[assignedTargets[assignedTargets != targetUtils.NULL_TARGET_INDEX]] = [1.0, 0.0, 0.0, 1.0]
    targetUtils.plotTargets(targetPositions, targetColors)

    print("Plotting time (s):", time.time() - start)
    plotUtils.pauseExecution()

