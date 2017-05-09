"""

Some utility methods related with the science targets.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T
  
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as collections
import time as time

import cobraUtils as cobraUtils
import benchUtils as benchUtils


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
        Target list (complex) or density (real scalar).
    bench: object
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

    # Calculate the distance matrix between the cobras and the targets positions
    distanceMatrix = np.abs(bench["center"][:, np.newaxis] - targets)
    
    # Set to zero the distances from those targets that cannot be reached by the cobras
    unreachableTargets = np.logical_or(distanceMatrix < bench["rMin"][:, np.newaxis], distanceMatrix > bench["rMax"][:, np.newaxis])
    distanceMatrix[unreachableTargets] = 0
    
    # Calculate the number of targets associated to eacg cobra
    nTargetsForCobra = np.sum(distanceMatrix > 0, axis=1)
    
    # Reorganize the distance information to have them ordered by distance to the cobra
    nCobras = len(bench["center"]) 
    maxTagetsForCobra = np.max(nTargetsForCobra)
    targetIndices = np.full((nCobras, maxTagetsForCobra), -1, dtype="int")
    targetDistances = np.zeros((nCobras, maxTagetsForCobra))

    for i in range(nCobras):
        # Get only those indices where the distance is not zero
        (validIndices,) = distanceMatrix[i].nonzero()

        # Sort the valid indices by their distance
        sortedIndices = validIndices[distanceMatrix[i, validIndices].argsort()]
        
        # Fill the target arrays
        targetIndices[i, :len(sortedIndices)] = sortedIndices
        targetDistances[i, :len(sortedIndices)] = distanceMatrix[i, sortedIndices]

    # Assign a targets to cobras looping from the closest targets to the more 
    # far away ones
    assignedTarget = np.full(nCobras, -1, dtype="int")
    
    for i in range(targetIndices.shape[1]):
        # Get a list with the unique targets in the given column 
        uniqueTargetIndices = np.unique(targetIndices[:, i])

        # Loop over the unique target indices
        for targetIndex in uniqueTargetIndices:    
            # Jump to the next target if the index does not represents a real target
            if targetIndex == -1:
                continue
            
            # Get the cobras for which this target is the closest in the current column
            (associatedCobras,) = np.where(targetIndices[:, i] == targetIndex)
            
            # Check how many associated cobras we have
            if len(associatedCobras) == 1:
                # Assign this cobra to the target
                assignedCobra = associatedCobras[0]
            else:
                # Get the number of target that can be reached by each associated cobra
                nTargets = nTargetsForCobra[associatedCobras]
            
                # Get the cobras for which this is the only target
                singleTargetCobras = associatedCobras[np.where(nTargets == 1)[0]]
   
                #              
                if len(singleTargetCobras) == 0:
                    # All cobras have multiple targets
                    # HACK SOLUTION: assign the target to the first cobra in the list
                    assignedCobra = associatedCobras[0]
                elif len(singleTargetCobras) == 1:
                    # Assign the target to the cobra that can only reach this target
                    assignedCobra = singleTargetCobras[0]
                else:
                    # Multiple cobras have only one choice.
                    # Assign the target to the closest cobra.
                    distances = targetDistances[singleTargetCobras, i]
                    assignedCobra = singleTargetCobras[np.where(distances == np.min(distances))[0]]
            
            # Assign the target to the correct cobra
            assignedTarget[assignedCobra] = targetIndex
            
            # Remove the target from the target arrays
            indicesToClean = np.where(targetIndices == targetIndex)
            targetIndices[indicesToClean] = -1
            targetDistances[indicesToClean] = 0
            
            # Make sure that we don't use the assign cobra anymore
            targetIndices[assignedCobra, :] = -1
            targetDistances[assignedCobra, :] = 0
            
            # Update the number of targets per cobra array
            nTargetsForCobra = np.sum(targetDistances > 0, axis=1)
    
    # Save the output data
    output = {}
    output["bench"] = bench
    output["tgt"] = targets
    output["nTargetsForCobra"] = nTargetsForCobra 
    output["targetIndices"] = targetIndices 
    output["targetDistances"] = targetDistances 
    output["assignedTarget"] = assignedTarget

    return output


def getCobraTargetConnections(centers, targets, indices):
    """Calculates the patrol areas for the bench cobras.

    Parameters
    ----------
    bench: Object
        The bench object. 
    
    Returns
    -------
    Object
        Pyhton tuple with the cobras inner and outer patrol area collections
         
    """
    delta = targets[indices] - centers
    arrowPatches = [patches.Arrow(c.real, c.imag, d.real, d.imag) for c, d, i in zip(centers, delta, indices) if i != -1]
    collection = collections.PatchCollection(arrowPatches, color="black")

    return collection

def plotTargets(targets, indices, bench):
    """Plots the target positions.

    Parameters
    ----------
    targets: Object
        A numpy complex array with the target positions. 
    
    """
    # Create the figure
    plt.figure("Target positions", facecolor="white", tight_layout=True, figsize=(7, 7))
    plt.title("Target positions")
    plt.xlabel("x position")
    plt.ylabel("y position")
    plt.show(block=False)

    # Fix the axes aspect ratio
    ax = plt.gca()
    ax.set_aspect("equal")
      
    # Set the axes limits
    limRange = 1.05 * bench["field"]["R"] * np.array([-1, 1]) 
    ax.set_xlim(bench["field"]["cm"].real + limRange)
    ax.set_ylim(bench["field"]["cm"].imag + limRange)

    # Obtain the cobras inner and outer patrol area collections
    (innerCollection, outerCollection) = benchUtils.getCobraPatrolAreas(bench)
    arrowCollection = getCobraTargetConnections(bench["center"], targets, indices)
            
    # Add the collections to the figure
    ax.add_collection(outerCollection)
    ax.add_collection(innerCollection)
    ax.add_collection(arrowCollection)

    # Draw the target positions and highlight those that can be reached by a cobra
    ax.scatter(np.real(targets), np.imag(targets), s=2)
    ax.scatter(np.real(targets[indices]), np.imag(targets[indices]), s=2, color="red")


if __name__ == "__main__":
    # Assign the target positions
    out = assignTargets(3.0, None)
        
    # Plot the target positions
    plotTargets(out["tgt"], out["assignedTarget"], out["bench"])
    plt.show()

