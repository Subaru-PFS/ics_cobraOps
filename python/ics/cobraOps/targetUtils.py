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
        The number of targets per patrol area.
    bench: object
        The bench to use.
    
    Returns
    -------
    Object
        Complex numpy array with the targets positions.
            
    """
    # Calculate the total number of targets based on the bench properties
    benchRadius = bench["field"]["R"]
    medianPatrolRadius = np.median(bench["rMax"])
    nTargets = int(np.ceil(density * (benchRadius / medianPatrolRadius) ** 2))
    
    # Calculate the uniformly distributed target positions
    ang = 2 * np.pi * np.random.random(nTargets)
    radius = benchRadius * np.sqrt(np.random.random(nTargets)) 
    targets = radius * np.exp(1j * ang)

    # Move the targets positions to the bench central position
    targets += bench["field"]["cm"]

    return targets


def assignTargets(targets, bench=None):
    """Assigns targets to cobras.

    Parameters
    ----------
    targets: object
        Target list (complex numpy array) or density (number).
    bench: object, optional
        The bench to use. If None, a full bench will be generated.
        Default is None.
    
    Returns
    -------
    Object
        ...
            
    """
    # Define the bench geometry if necessary
    if bench is None:
        # Get the cobras central positions for the full PFI
        centers = cobraUtils.getPFICenters()

        # Define the bench geometry
        bench = benchUtils.defineBenchGeometry(centers, True, True)

        # Set the bench alpha value to zero
        bench["alpha"] = 0

    # Generate the targets if a density has been provided
    if isinstance(targets, int) or isinstance(targets, float):
        targetDensity = targets
        targets = generateTargets(targetDensity, bench)

    # Calculate the distance matrix between the cobras center and the targets positions
    distanceMatrix = np.abs(bench["center"][:, np.newaxis] - targets)
    
    # Set to zero the distances from those targets that cannot be reached by the cobras
    unreachableTargets = np.logical_or(distanceMatrix < bench["rMin"][:, np.newaxis], distanceMatrix > bench["rMax"][:, np.newaxis])
    distanceMatrix[unreachableTargets] = 0
    
    # Calculate the number of targets associated to each cobra
    nTargetsPerCobra = np.sum(distanceMatrix > 0, axis=1)
    
    # Reorganize the distance information to have them ordered by distance to the cobra
    nCobras = len(bench["center"]) 
    maxTagetsPerCobra = np.max(nTargetsPerCobra)
    targetIndices = np.full((nCobras, maxTagetsPerCobra), -1, dtype="int")
    targetDistances = np.zeros((nCobras, maxTagetsPerCobra))

    for i in range(nCobras):
        # Get the target distances to this cobra
        distances = distanceMatrix[i]

        # Select only those target indices where the distance is not zero
        (validTargetIndices,) = distances.nonzero()

        # Sort the valid target indices by their distance to the cobra
        sortedTargetIndices = validTargetIndices[distances[validTargetIndices].argsort()]
        
        # Fill the target arrays
        targetIndices[i, :len(sortedTargetIndices)] = sortedTargetIndices
        targetDistances[i, :len(sortedTargetIndices)] = distances[sortedTargetIndices]

    # Assign targets to cobras looping from the closest to the more far away ones
    assignedTarget = np.full(nCobras, -1, dtype="int")
        
    for i in range(maxTagetsPerCobra):
        # Consider only cobras that do not have an assigned target
        freeCobras = assignedTarget < 0

        # Get a list with the unique targets in the given column 
        uniqueTargetIndices = np.unique(targetIndices[freeCobras, i])

        # Loop over the unique target indices
        for targetIndex in uniqueTargetIndices:    
            # Jump to the next target if the index does not represents a real target
            if targetIndex == -1:
                continue

            # Get the free cobras for which this target is the closest in the current column
            (associatedCobras,) = np.where(np.logical_and(targetIndices[:, i] == targetIndex, freeCobras))
            
            # Check how many associated cobras we have
            if len(associatedCobras) == 1:
                # Use this single cobra for this target
                cobraToUse = associatedCobras[0]
            else:
                # Select the cobras for which this is the only target
                singleTargetCobras = associatedCobras[nTargetsPerCobra[associatedCobras] == 1]
                
                # Decide depending on home many of these cobras we have
                if len(singleTargetCobras) == 0:
                    # All cobras have multiple targets
                    # HACK SOLUTION: assign the target to the first cobra in the list
                    cobraToUse = associatedCobras[0]
                elif len(singleTargetCobras) == 1:
                    # Assign the target to the cobra that can only reach this target
                    cobraToUse = singleTargetCobras[0]
                else:
                    # Assign the target to the closest cobra
                    distances = targetDistances[singleTargetCobras, i]
                    cobraToUse = singleTargetCobras[distances == np.min(distances)]
            
            # Assign the target to the correct cobra
            assignedTarget[cobraToUse] = targetIndex
            
            # Remove the target from the target arrays
            indicesToClean = np.where(targetIndices == targetIndex)
            targetIndices[indicesToClean] = -1
            targetDistances[indicesToClean] = 0
            
            # Update the number of targets per cobra array
            nTargetsPerCobra = np.sum(targetDistances > 0, axis=1)
    
    # Calculate the cobra collisions
    collisions = nCobras
    
    while collisions > 0:
        # Calculate the expected cobra positions (leave at home unused cobras) 
        cobraPositions = bench["home0"].copy()
        usedCobras = assignedTarget >= 0
        cobraPositions[usedCobras] = targets[assignedTarget[usedCobras]]  
    
        # Calculate the collision distance matrix between the cobra positions 
        # and the nearby cobra arm positions.
        collisionMatrix = calculateCollisionMatrix(cobraPositions, bench)

        # Get the cobra collisions for the current configuration
        tooClose = np.logical_and(collisionMatrix > 0, collisionMatrix < bench["minDist"])
        (problematicCobras, neabyCobras) = np.where(tooClose)

        # Solve the binary collisions
        collisions = 0
    
    
    # Save the output data
    output = {}
    output["bench"] = bench
    output["tgt"] = targets
    output["nTargetsForCobra"] = nTargetsPerCobra 
    output["targetIndices"] = targetIndices 
    output["targetDistances"] = targetDistances 
    output["assignedTarget"] = assignedTarget
    output["collisions"] = problematicCobras

    return output


def calculateCollisionMatrix(positions, bench):
    """Calculates the collision distance matrix between the cobra positions 
    and the nearby cobra arm positions.
    
    Parameters
    ----------
    positions: object
        Complex numpy array with the cobras commanded positions.
    bench: object
        The bench geometry

    Returns
    -------
    object
        The collision distance matrix.
    
    """  
    # Calculate the cobras rotation angles to reach the given position
    (tht, phi) = getCobraRotationAngles(positions - bench["center"], bench["L1"], bench["L2"])
    
    # Compute the cobras elbow positions
    elbows = positions - bench["L2"] * np.exp(1j * (tht + phi))
  
    # Use the bench precalculated nearest neighbors information
    distanceMatrix = np.zeros(bench["nnMap"].shape)
    cobras = bench["NN"]["row"]
    nearbyCobras = bench["NN"]["col"]
    
    targetPositions = positions[cobras]
    nearbyTargetPositions = positions[nearbyCobras]
    nearbyElbowPositions = elbows[nearbyCobras]

    # 
    distances = distanceToLineSegment(targetPositions, nearbyElbowPositions, nearbyTargetPositions)
    distanceMatrix[cobras, nearbyCobras] = distances
    
    return distanceMatrix


def getCobraRotationAngles(targetPostions, link1, link2, useNegativePhi=True):
    """Calculates the cobra rotation angles to reach the given target position. 
    
    The code assumes that the cobras can reach the given positions.
    
    Parameters
    ----------
    targetPositions: object
        Complex numpy array with the target coordinates centered at the cobra
        position.
    link1: object
        Numpy array or constant with the links1 lengths.
    link2: object
        Numpy array or constant with the links2 lengths.
    useNegativePhi: bool, optional
        If True the phi values will be negative. If False, they will be
        positive. Default is True.

    Returns
    -------
    tuple
        A python tuple with two numpy arrays containing the (tht, phi) cobra 
        rotation angles. 
    
    """  
    # Calculate the cobra angles applying the law of cosines
    distance = np.abs(targetPostions)
    distanceSq = distance ** 2
    link1Sq = link1 ** 2
    link2Sq = link2 ** 2    
    phiSign = 1 - 2 * useNegativePhi
    phi = phiSign * np.arccos((distanceSq - link1Sq - link2Sq) / (2 * link1 * link2))
    tht = np.angle(targetPostions) - phiSign * np.arccos(-(link2Sq - link1Sq - distanceSq) / (2 * link1 * distance))
    
    # Force tht go from -pi to pi instead of from 0 to 2pi
    tht = (tht - np.pi) % (2 * np.pi) - np.pi

    return (tht, phi)


def distanceToLineSegment(points, startPoints, endPoints):
    """Calculates the minimum distances between points and line segments.
    
    Parameters
    ----------
    points: object
        Complex numpy array with the point coordinates. 
    startPoints: object
        Complex numpy array with the lines start coordinates. 
    endPoints: object
        Complex numpy array with the lines end coordinates. 

    Returns
    -------
    object
        A numpy array with the point to line distances.
        
    """
    # Translate the points and the line end points to the lines starting points
    translatedPoints = points - startPoints
    translatedEndPoints = endPoints - startPoints
    
    # Rotate the translated points to have the lines on the x axis
    rotatedPoints = translatedPoints * np.exp(-1j * np.angle(translatedEndPoints))
    
    # Define 3 regions for the points: left of origin, over the lines, right of the lines
    x = np.real(rotatedPoints)
    lineLengths = np.abs(translatedEndPoints)
    (region1,) = np.where(x <= 0)
    (region2,) = np.where(np.logical_and(x > 0 , x < lineLengths))
    (region3,) = np.where(x >= lineLengths)

    # Calculate the distances in each region
    distances = np.empty(len(points))
    distances[region1] = np.abs(rotatedPoints[region1])
    distances[region2] = np.abs(np.imag(rotatedPoints[region2]))
    distances[region3] = np.abs(rotatedPoints[region3] - lineLengths[region3])

    return distances


def plotCollisions(targets, indices, colIndices, bench):
    """Plots the cobra collisions.

    Parameters
    ----------
    targets: object
        A numpy complex array with the target positions. 
    
    """
    # Create the figure
    plotUtils.createNewFigure("Cobra collisions", "x position", "y position")
      
    # Set the axes limits
    limRange = 1.05 * bench["field"]["R"] * np.array([-1, 1]) 
    xLim = bench["field"]["cm"].real + limRange
    yLim = bench["field"]["cm"].imag + limRange
    plotUtils.setAxesLimits(xLim, yLim)
    
    # Plot the cobra patrol areas using ring shapes and use a 
    # different color for those with detected collisions
    centers = bench["center"]
    rMin = bench["rMin"]
    rMax = bench["rMax"]
    colors = np.full((len(centers), 4), [0.0, 0.0, 1.0, 0.15])
    colors[colIndices] = [0.0, 1.0, 0.0, 0.5]
    plotUtils.addRings(centers, rMin, rMax, facecolors=colors)
 
    # Draw the cobras positions using a different color for those
    # that do not have an assigned target
    (assignedCobras,) = np.where(indices >= 0)
    positions = bench["home0"]
    positions[assignedCobras] = targets[indices[assignedCobras]]
    cobraCenters = bench["center"]
    L1 = bench["L1"]
    L2 = bench["L2"]
    (tht, phi) = getCobraRotationAngles(positions - cobraCenters, L1, L2)
    link1 = cobraCenters + L1 * np.exp(1j * tht)
    link2 = link1 + L2 * np.exp(1j * (tht + phi))
    colors = np.full((len(positions), 4), [1.0, 0.0, 0.0, 0.25])
    colors[assignedCobras] = [0.0, 0.0, 1.0, 0.5]
    plotUtils.addLines(cobraCenters, link1, edgecolor=colors, linewidths=2)
    plotUtils.addThickLines(link1, link2, 0.5 * bench["minDist"], facecolors=colors)  

    # Draw the target positions and highlight those that are assigned to a cobra
    plotUtils.addPoints(targets, s=2, facecolor="0.4")
    plotUtils.addPoints(targets[indices[assignedCobras]], s=2, facecolor="red")


if __name__ == "__main__":
    # Assign the target positions
    start = time.time()   
    out = assignTargets(1.0, None)
    print("Target calculation time:" + str(time.time() - start))

    # Plot the target positions
    start = time.time()
    plotCollisions(out["tgt"], out["assignedTarget"], out["collisions"], out["bench"])
    print("Ploting time:" + str(time.time() - start))
    plotUtils.pauseExecution()

