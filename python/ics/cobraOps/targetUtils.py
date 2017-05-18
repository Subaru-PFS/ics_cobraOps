"""

Some utility methods related with the science targets.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np
import time as time

import cobraUtils as cobraUtils
import benchUtils as benchUtils
import plotUtils as plotUtils


def generateTargets(density, bench):
    """Generates a set of targets uniformly distributed over the bench field of view.

    Parameters
    ----------
    density: float
        The average number of targets per patrol area.
    bench: object
        The bench geometry to use.

    Returns
    -------
    Object
        Complex numpy array with the targets positions.

    """
    # Calculate the total number of targets based on the bench properties
    benchCenter = np.mean(bench["center"])
    benchRadius = np.max(np.abs(bench["center"] - benchCenter) + bench["rMax"])
    medianPatrolRadius = np.median(bench["rMax"])
    nTargets = int(np.ceil(density * (benchRadius / medianPatrolRadius) ** 2))
 
    # Calculate the uniformly distributed target positions
    ang = 2 * np.pi * np.random.random(nTargets)
    radius = benchRadius * np.sqrt(np.random.random(nTargets)) 
    targetPositions = radius * np.exp(1j * ang)

    # Move the targets positions to the bench central position
    targetPositions += benchCenter

    return targetPositions


def assignTargets(targetPositions, bench):
    """Assigns a set of targets to the cobras in the bench.

    Parameters
    ----------
    targetPositions: object
        A complex numpy array with the targets coordinates.
    bench: object
        The bench geometry to use.
    
    Returns
    -------
    tuple
        A python tuple with the assigned target indices and the cobra positions.

    """
    # Get the indices and distances of those targets that can be reached by each cobra
    (targetIndices, targetDistances) = getAccesibleTargets(targetPositions, bench)
    
    # Assign a single target to each cobra based of the target distances
    assignedTargets = assignTargetsByDistance(targetIndices, targetDistances)
    
    # Calculate the cobra positions solving possible collisions between cobras
    cobraPositions = solveCobraCollisions(assignedTargets, targetIndices, targetDistances, targetPositions, bench)

    return (assignedTargets, cobraPositions)


def getAccesibleTargets(targetPositions, bench):
    """Returns the targets that each cobra can reach ordered by distance.

    Parameters
    ----------
    targetPositions: object
        A complex numpy array with the targets coordinates.
    bench: object
        The bench geometry to use.
    
    Returns
    -------
    tuple
        A python tuple with the indices and distances of each target that can
        be reached by the cobra.

    """
    # Calculate the distance matrix between the cobra centers and the target
    # positions: select first by the x axis distance, then by the y axis 
    # distance, and finally by the radial distance.
    distanceMatrix = np.abs(bench["center"].real[:, np.newaxis] - targetPositions.real)
    (cobras, targets) = np.where(distanceMatrix < bench["rMax"][:, np.newaxis])
    yDistance = np.abs(bench["center"][cobras].imag - targetPositions[targets].imag)
    distanceMatrix[cobras, targets] = yDistance
    (indx,) = np.where(yDistance < bench["rMax"][cobras])
    cobras = cobras[indx]
    targets = targets[indx]
    distanceMatrix[cobras, targets] = np.abs(bench["center"][cobras] - targetPositions[targets])

    # This is simpler method, but it's two to three times slower.
    # distanceMatrixOld = np.abs(bench["center"][:, np.newaxis] - targetPositions)
    
    # Check which targets can be reached by each cobra
    reachableTargets = np.logical_and(distanceMatrix > bench["rMin"][:, np.newaxis], distanceMatrix < bench["rMax"][:, np.newaxis])

    # Calculate the total number of targets that each cobra can reach
    nTargetsPerCobra = np.sum(reachableTargets, axis=1)

    # Order the targets by their distance to the cobra
    nCobras = len(bench["center"])
    maxTagetsPerCobra = np.max(nTargetsPerCobra)
    targetIndices = np.full((nCobras, maxTagetsPerCobra), -1, dtype="int")
    targetDistances = np.zeros((nCobras, maxTagetsPerCobra))

    for i in range(nCobras):
        # Get the indices of the targets that can be reached by the cobra
        (validTargetIndices,) = np.where(reachableTargets[i])

        # Sort the valid target indices by their distance to the cobra
        distances = distanceMatrix[i]
        orderedTargetIndices = validTargetIndices[distances[validTargetIndices].argsort()]
        
        # Fill the target arrays
        targetIndices[i, :nTargetsPerCobra[i]] = orderedTargetIndices
        targetDistances[i, :nTargetsPerCobra[i]] = distances[orderedTargetIndices]

    return (targetIndices, targetDistances)


def assignTargetsByDistance(targetIndices, targetDistances):
    """Assigns a single target to each cobra based on the target distance.
    
    This method assumes that the input arrays are ordered by the target 
    distance to the center of the cobra.
    
    This method changes the input arrays, removing references to the assigned 
    cobras.

    Parameters
    ----------
    targetIndices: object
        A numpy array with the indices of the targets that can be reached by
        a given cobra.
    targetDistances: object
        A numpy array with the distances of the targets that can be reached by
        a given cobra.
    
    Returns
    -------
    tuple
        A numpy array with the indices of the targets assigned to each cobra.

    """
    # Assign targets to cobras looping from the closest to the more far away ones
    nCobras = targetIndices.shape[0]
    maxTargetsPerCobra = targetIndices.shape[1]
    assignedTargets = np.full(nCobras, -1, dtype="int")       
    freeCobras = assignedTargets < 0
    
    for i in range(maxTargetsPerCobra):
        # Get a list with the unique targets in the given column
        columnTargetIndices = targetIndices[:, i]
        uniqueTargetIndices = np.unique(columnTargetIndices[freeCobras])

        # Loop over the unique target indices
        for targetIndex in uniqueTargetIndices:    
            # Jump to the next target if the index does not represents a real target
            if targetIndex == -1:
                continue

            # Get the free cobras for which this target is the closest in the current column
            (associatedCobras,) = np.where(np.logical_and(columnTargetIndices == targetIndex, freeCobras))
            
            # Check how many associated cobras we have
            if len(associatedCobras) == 1:
                # Use this single cobra for this target
                cobraToUse = associatedCobras[0]
            else:
                # Select the cobras for which this is the only target
                nAvailableTargets = np.sum(targetIndices[associatedCobras] != -1, axis=1)
                singleTargetCobras = associatedCobras[nAvailableTargets == 1]
                
                # Decide depending on how many of these cobras we have
                if len(singleTargetCobras) == 0:
                    # All cobras have multiple targets: select the closest cobra
                    distances = targetDistances[associatedCobras, i]
                    cobraToUse = associatedCobras[distances.argmin()]
                elif len(singleTargetCobras) == 1:
                    # Assign the target to the cobra that can only reach this target
                    cobraToUse = singleTargetCobras[0]
                else:
                    # Assign the target to the closest single target cobra
                    distances = targetDistances[singleTargetCobras, i]
                    cobraToUse = singleTargetCobras[distances.argmin()]

            
            # Assign the target to the correct cobra
            assignedTargets[cobraToUse] = targetIndex
            freeCobras[cobraToUse] = False
            
            # Remove the target from the target indices array
            targetIndices[targetIndices == targetIndex] = -1

    return assignedTargets


def solveCobraCollisions(assignedTargets, targetIndices, targetDistances, targetPositions, bench):
    """Detects and solves cobra collisions assigning them alternative targets.

    This method assumes that the targetIndices and targetDistances arrays are 
    ordered by the target distance to the center of the cobra.

    This method changes the assignedTargets, targetIndices and targetDistances 
    input arrays.

    Parameters
    ----------
    assignedTargets: object
        A numpy array with the indices of the targets assigned to each cobra.
    targetIndices: object
        A numpy array with the indices of all other alternative targets that
        can be assigned to each cobra.
    targetDistances: object
        A numpy array with the distances of all the alternative targets that
        can be assigned to each cobra.
    targetPositions: object
        A complex numpy array with the targets coordinates.
    bench: object
        The bench geometry to use.
    
    Returns
    -------
    tuple
        A complex numpy array with the cobras positions.

    """
    # Set the cobra positions to their associated target positions, leaving 
    # unused cobras at their home positions 
    cobraPositions = bench["home0"].copy()
    usedCobras = assignedTargets >= 0
    cobraPositions[usedCobras] = targetPositions[assignedTargets[usedCobras]]  

    # Get the indices of the cobras where we have a collision
    (problematicCobras, nearbyProblematicCobras) = getProblematicCobras(cobraPositions, bench)
 
    # Try to solve the collisions one by one
    for c, nc in zip(problematicCobras, nearbyProblematicCobras):
        # We only need to solve the first half of the problematic cases
        if nc > c:
            # Check if one of the colliding cobras is unused
            if assignedTargets[c] < 0 or assignedTargets[nc] < 0:
                # The unused cobra is the cobra that we are going to move
                if assignedTargets[c] < 0:
                    cobraToMove = c
                else:
                    cobraToMove = nc

                # Calculate the initial number of collisions associated with that cobra 
                collisions = getCollisionsForCobra(cobraToMove, cobraPositions, bench)

                # Rotate the cobra until we find the minimum number of collisions
                bestPosition = cobraPositions[cobraToMove]
            
                for ang in range(5):
                    # Rotate the cobra 60 degrees around its center
                    cobraCenter = bench["center"][cobraToMove]
                    cobraPositions[cobraToMove] = (cobraPositions[cobraToMove] - cobraCenter) * np.exp(1j * np.pi / 3) + cobraCenter

                    # Calculate the number of collisions at the current position 
                    currentCollisions = getCollisionsForCobra(cobraToMove, cobraPositions, bench)
                    
                    # Check if the number of collisions decreased
                    if currentCollisions < collisions:
                        # Save the information from this position
                        bestPosition = cobraPositions[cobraToMove]
                        collisions = currentCollisions
 
                    # Exit the loop if the number of collisions is already zero
                    if collisions == 0:
                        break

                # Use the cobra position where we had less collisions
                cobraPositions[cobraToMove] = bestPosition
            else:
                # Calculate the initial number of collisions associated with the two cobras
                collisions = getCollisionsForCobra(c, cobraPositions, bench)
                collisions += getCollisionsForCobra(nc, cobraPositions, bench)

                # Get the targets that can be reached by the each cobra
                initialTarget1 = assignedTargets[c]
                initialTarget2 = assignedTargets[nc]
                targets1 = np.append(initialTarget1, targetIndices[c][targetIndices[c] >= 0])
                targets2 = np.append(initialTarget2, targetIndices[nc][targetIndices[nc] >= 0])
                
                # Create two arrays reflecting all the possible target combinations
                targetsCombination1 = np.repeat(targets1, len(targets2))
                targetsCombination2 = np.tile(targets2, len(targets1))
                
                # Exclude the current target combination
                targetsCombination1 = targetsCombination1[1:]
                targetsCombination2 = targetsCombination2[1:]
                
                # Loop over all the possible combinations until we find the minimum 
                # number of collisions
                bestTarget1 = initialTarget1    
                bestTarget2 = initialTarget2    
                bestPosition1 = cobraPositions[c]
                bestPosition2 = cobraPositions[nc]
 
                for newTarget1, newTarget2 in zip(targetsCombination1, targetsCombination2):
                    # Assign the new cobra positions
                    cobraPositions[c] = targetPositions[newTarget1]
                    cobraPositions[nc] = targetPositions[newTarget2]
                    
                    # Calculate the number of collisions at the current positions 
                    currentCollisions = getCollisionsForCobra(c, cobraPositions, bench)
                    currentCollisions += getCollisionsForCobra(nc, cobraPositions, bench)
                               
                    # Check if the number of collisions decreased significantly
                    # A decrease of one means that we solved the current collision,
                    # but we created a new collision with another nearby cobra.
                    if currentCollisions <= collisions - 2:
                        # Save the information from these targets
                        bestTarget1 = newTarget1
                        bestTarget2 = newTarget2
                        bestPosition1 = cobraPositions[c]
                        bestPosition2 = cobraPositions[nc]
                        collisions = currentCollisions
                    
                    # Exit the loop if the number of collisions is already zero
                    if collisions == 0:
                        break
                
                # Use the targets where we had less collisions
                assignedTargets[c] = bestTarget1
                assignedTargets[nc] = bestTarget2
                cobraPositions[c] = bestPosition1
                cobraPositions[nc] = bestPosition2
                
                # Swap the initial and best targets information in the targets arrays
                (index,) = np.where(targetIndices[c] == bestTarget1)
                targetIndices[c, index] = initialTarget1
                targetDistances[c, index] = np.abs(bench["center"][c] - targetPositions[initialTarget1])
                (index,) = np.where(targetIndices[nc] == bestTarget2)
                targetIndices[nc, index] = initialTarget2
                targetDistances[nc, index] = np.abs(bench["center"][nc] - targetPositions[initialTarget2])
               
                # Remove the new assigned targets from the target indices array
                # to avoid that other cobras can use them
                targetIndices[np.logical_or(targetIndices == bestTarget1, targetIndices == bestTarget2)] = -1
    
    return cobraPositions


def getCollisionsForCobra(cobraIndex, cobraPositions, bench):
    """Calculates the number of collisions of a given cobra with its neighbors.
    
    Parameters
    ----------
    cobraIndex: int
        The index of the cobra for which we want to count the cobra collisions.
    cobraPositions: object
        A complex numpy array with the cobras positions.
    bench: object
        The bench geometry to use.

    Returns
    -------
    int
        The number of collisions between this cobra and its neighbors.
    
    """  
    # Get the cobra neighbors from the bench precalculated information
    nearbyCobras = bench["NN"]["col"][bench["NN"]["row"] == cobraIndex]

    # Calculate the cobras rotation angles to reach the given positions
    allIndices = np.append(cobraIndex, nearbyCobras)
    cobraCenters = bench["center"][allIndices]
    L1 = bench["L1"][allIndices]
    L2 = bench["L2"][allIndices]
    (tht, phi) = getCobraRotationAngles(cobraPositions[allIndices] - cobraCenters, L1, L2)

    # Compute the cobras elbow positions
    cobraElbows = cobraCenters + L1 * np.exp(1j * tht)

    # Calculate the distances between the cobra link and the nearby cobras links
    startPoints1 = np.repeat(cobraPositions[cobraIndex], len(nearbyCobras))
    endPoints1 = np.repeat(cobraElbows[0], len(nearbyCobras))
    startPoints2 = cobraPositions[nearbyCobras]
    endPoints2 = cobraElbows[1:]
    distances = distanceBetweenLineSegments(startPoints1, endPoints1, startPoints2, endPoints2)

    # Get the cobra collisions for the current configuration
    cobraCollisions = distances < (bench["minDist"][cobraIndex] + bench["minDist"][nearbyCobras]) / 2

    return np.sum(cobraCollisions)


def getProblematicCobras(cobraPositions, bench):
    """Obtains the indices of the cobras involved in collisions.
    
    Parameters
    ----------
    cobraPositions: object
        A complex numpy array with the cobras positions.
    bench: object
        The bench geometry to use.

    Returns
    -------
    tuple
        A python tuple with the indices of the cobras involved in collisions
        and the indices of the nearby cobras with which they collide.
    
    """
    # Calculate the cobras rotation angles to reach the given positions
    cobraCenters = bench["center"]
    L1 = bench["L1"]
    L2 = bench["L2"]
    (tht, phi) = getCobraRotationAngles(cobraPositions - cobraCenters, L1, L2)
    
    # Compute the cobras elbow positions
    cobraElbows = cobraCenters + L1 * np.exp(1j * tht)
  
    # Get the bench precalculated nearest neighbors information
    cobras = bench["NN"]["row"]
    nearbyCobras = bench["NN"]["col"]
    
    # Calculate the distances between the cobras links and the nearby cobras links
    startPoints1 = cobraPositions[cobras]
    endPoints1 = cobraElbows[cobras]
    startPoints2 = cobraPositions[nearbyCobras]
    endPoints2 = cobraElbows[nearbyCobras]
    distances = distanceBetweenLineSegments(startPoints1, endPoints1, startPoints2, endPoints2)

    # Get the cobra collisions for the current configuration
    cobraCollisions = distances < (bench["minDist"][cobras] + bench["minDist"][nearbyCobras]) / 2

    return (cobras[cobraCollisions], nearbyCobras[cobraCollisions])


def getCobraRotationAngles(relativePositions, link1, link2, useNegativePhi=True):
    """Calculates the cobra rotation angles to reach some relative positions. 
    
    The code assumes that the cobras can reach the given positions.
    
    Parameters
    ----------
    relativePositions: object
        A complex numpy array with the position coordinates relative to the
        cobras centers.
    link1: object
        A numpy array or constant with the links1 lengths.
    link2: object
        A numpy array or constant with the links2 lengths.
    useNegativePhi: bool, optional
        If True the phi angle values will be negative. If False, the phi 
        angles will be positive. Default is True.

    Returns
    -------
    tuple
        A python tuple with the cobras rotation angles (tht, phi). 
    
    """  
    # Calculate the cobra angles applying the law of cosines
    distance = np.abs(relativePositions)
    distanceSq = distance ** 2
    link1Sq = link1 ** 2
    link2Sq = link2 ** 2    
    phiSign = 1 - 2 * useNegativePhi
    phi = phiSign * np.arccos((distanceSq - link1Sq - link2Sq) / (2 * link1 * link2))
    tht = np.angle(relativePositions) - phiSign * np.arccos(-(link2Sq - link1Sq - distanceSq) / (2 * link1 * distance))
    
    # Force tht to go from -pi to pi, instead of from 0 to 2pi
    tht = (tht - np.pi) % (2 * np.pi) - np.pi

    return (tht, phi)


def distanceBetweenLineSegments(startPoints1, endPoints1, startPoints2, endPoints2):
    """Calculates the minimum distances between line segments.
    
    Parameters
    ----------
    startPoints1: object
        A complex numpy array with the first line segment start coordinates.
    endPoints1: object
        A complex numpy array with the first line segment end coordinates. 
    startPoints2: object
        A complex numpy array with the second line segment start coordinates. 
    endPoints2: object
        A complex numpy array with the second line segment end coordinates. 

    Returns
    -------
    object
        A numpy array with the minimum distance between the line segments.
        
    """
    # Calculate the minimum distances for each point to segment combination
    distances1 = distanceToLineSegment(startPoints1, startPoints2, endPoints2)
    distances2 = distanceToLineSegment(endPoints1, startPoints2, endPoints2)
    distances3 = distanceToLineSegment(startPoints2, startPoints1, endPoints1)
    distances4 = distanceToLineSegment(endPoints2, startPoints1, endPoints1)

    # Return the minimum distances
    return np.min((distances1, distances2, distances3, distances4), axis=0)


def distanceToLineSegment(points, startPoints, endPoints):
    """Calculates the minimum distances between points and line segments.
    
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
        A numpy array with the minimum distances between the points and the 
        line segments.
        
    """
    # Translate the points and the line segment end points to the line segment
    # starting points
    translatedPoints = points - startPoints
    translatedEndPoints = endPoints - startPoints
    
    # Rotate the translated points to have the line segment on the x axis
    rotatedPoints = translatedPoints * np.exp(-1j * np.angle(translatedEndPoints))
    
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
    distances[region3] = np.abs(rotatedPoints[region3] - lineLengths[region3])

    return distances


def plotCobraTargetAssociations(cobraPositions, problematicCobras, assignedTargets, targetPositions, bench):
    """Plots the cobra-target associations.

    Parameters
    ----------
    cobraPositions
        A complex numpy array with the cobra positions.
    problematicCobras: object
        A numpy array with the indices of the cobras involved on a collision.
    assignedTargets: object
        A numpy array with the indices of the targets assigned to each cobra.
    targetPositions: object
        A complex numpy array with the targets coordinates.
    bench: object
        The bench geometry to use.
    
    """
    # Create the figure
    plotUtils.createNewFigure("Cobra-target associations", "x position", "y position")
      
    # Set the axes limits
    benchCenter = np.mean(bench["center"])
    benchRadius = np.max(np.abs(bench["center"] - benchCenter) + bench["rMax"])
    limRange = 1.05 * benchRadius * np.array([-1, 1]) 
    xLim = benchCenter.real + limRange
    yLim = benchCenter.imag + limRange
    plotUtils.setAxesLimits(xLim, yLim)
    
    # Plot the cobra patrol areas using ring shapes and use a 
    # different color for those with detected collisions
    cobraCenters = bench["center"]
    rMin = bench["rMin"]
    rMax = bench["rMax"]
    colors = np.full((len(cobraCenters), 4), [0.0, 0.0, 1.0, 0.15])
    colors[problematicCobras] = [0.0, 1.0, 0.0, 0.5]
    plotUtils.addRings(cobraCenters, rMin, rMax, facecolors=colors)
 
    # Draw the cobras positions using a different color for those
    # that do not have an assigned target
    L1 = bench["L1"]
    L2 = bench["L2"]
    (tht, phi) = getCobraRotationAngles(cobraPositions - cobraCenters, L1, L2)
    link1 = cobraCenters + L1 * np.exp(1j * tht)
    link2 = link1 + L2 * np.exp(1j * (tht + phi))
    colors = np.full((len(cobraPositions), 4), [1.0, 0.0, 0.0, 0.25])
    colors[assignedTargets >= 0] = [0.0, 0.0, 1.0, 0.5]
    plotUtils.addLines(cobraCenters, link1, edgecolor=colors, linewidths=2)
    plotUtils.addThickLines(link1, link2, 0.5 * bench["minDist"], facecolors=colors)  

    # Draw the target positions and highlight those that are assigned to a cobra
    plotUtils.addPoints(targetPositions, s=2, facecolor="0.4")
    plotUtils.addPoints(targetPositions[assignedTargets[assignedTargets >= 0]], s=2, facecolor="red")


if __name__ == "__main__":
    # Define the target density to use
    targetDensity = 5
    
    # Get the cobras central positions for the full PFI
    start = time.time()
    centers = cobraUtils.getPFICenters()
    print("Number of cobras: " + str(len(centers)))

    # Define the bench geometry
    bench = benchUtils.defineBenchGeometry(centers, True, True)

    # Create a random sample of targets
    targetPositions = generateTargets(targetDensity, bench)
    print("Number of simulated targets: " + str(len(targetPositions)))

    # Assign the target positions
    (assignedTargets, cobraPositions) = assignTargets(targetPositions, bench)
    
    # Get the cobra for which the collision could not solved
    (problematicCobras, nearbyProblematicCobras) = getProblematicCobras(cobraPositions, bench)
    print("Number of unsolved collisions: " + str(len(problematicCobras)/2))
    print("Total computation time (s): " + str(time.time() - start))

    # Plot the cobra-target associations
    start = time.time()
    plotCobraTargetAssociations(cobraPositions, problematicCobras, assignedTargets, targetPositions, bench)
    print("Plotting time (s): " + str(time.time() - start))
    plotUtils.pauseExecution()

