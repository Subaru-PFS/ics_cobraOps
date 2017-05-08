"""

Some utility methods related with the science targets.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T
  
"""

import numpy as np
import matplotlib.pyplot as plt

import cobraUtils as cobraUtils
import benchUtils as benchUtils


def generateTargets(density, bench):
    """Generates a set of targets uniformly distributed over the bench field of view.

    Parameters
    ----------
    density: float
        The number of targets per patrol area.
    bench: object
        The bench to use. If None, a full bench will be generated using a 
        calibration file.
    
    Returns
    -------
    Object
        Complex numpy array with the targets positions.
            
    """
    # Define the bench geometry if necessary
    if bench is None:
        # Get the cobras central positions for the full PFI
        centers = cobraUtils.getPFICenters()

        # Define the bench geometry
        bench = benchUtils.defineBenchGeometry(centers, True, True)

        # Set the bench alpha value to zero
        bench["alpha"] = 0
    
    # Get the total number of cobras
    nCobras = len(bench["center"]) 

    # Calculate the total number of targets based on the bench properties
    benchRadius = bench["field"]["R"]
    medianPatrolRadius = np.median(bench["rMax"])
    nTargets = int(np.ceil(density * (benchRadius / medianPatrolRadius) ** 2))
    
    # Calculate the uniformly distributed target positions
    print("Generating " + str(nTargets) + " targets for " + str(nCobras) + " positioners")
    ang = 2 * np.pi * np.random.random(nTargets)
    radius = benchRadius * np.sqrt(np.random.random(nTargets)) 
    targets = radius * np.exp(1j * ang)

    # Move the targets positions to the bench central position
    targets += bench["field"]["cm"]

    return targets


def assignTargets(tgt, bench=None, makefigs=False):
    """Assigns targets to cobras.

    Parameters
    ----------
    tgt: object
        Target list (complex) or density (real scalar).
    bench: object
        The bench to use. If None, a full bench will be generated.
        Default is None.
    makefigs: bool
        If true, some figures will be made. Default is False.
    
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
    if isinstance(tgt, int) or isinstance(tgt, float):
        targetDensity = tgt
        tgt = generateTargets(targetDensity, bench)

    # Get the total number of cobras and targets
    nCobras = len(bench["center"]) 
    nTargets = len(tgt)

    # Calculate the distance matrix between the cobras and the targets positions
    distanceMatrix = np.zeros((nCobras, nTargets))

    for i in range(nCobras):
        # Calculate the distance between the cobra and all the target positions
        distanceMatrix[i] = np.abs(tgt - bench["center"][i])
        
        # Set to zero the distances from those targets that cannot be reached by the cobra
        unreachableTargets = np.logical_or(distanceMatrix[i] < bench["rMin"][i], distanceMatrix[i] > bench["rMax"][i])
        distanceMatrix[i, unreachableTargets] = 0   
    
    # Calculate the number of targets per cobra and the number of cobras per target
    nCobrasForTarget = np.sum(distanceMatrix > 0, axis=0)
    nTargetsForCobra = np.sum(distanceMatrix > 0, axis=1)
    
    # Reorganize the distance information to have them ordered by distance to the cobra
    maxTagetsForCobra = np.max(nTargetsForCobra)
    targetIndices = np.zeros((nCobras, maxTagetsForCobra), dtype="int")
    targetDistances = np.zeros((nCobras, maxTagetsForCobra))
    
    for i in range(nCobras):
        # Get only those indices where the distance is not zero
        validIndices = distanceMatrix[i].nonzero()[0]
        
        # Get the sorted valid indices
        sortedIndices = validIndices[np.argsort(distanceMatrix[i, validIndices])]
        
        # Fill the arrays
        targetIndices[i, :len(sortedIndices)] = sortedIndices
        targetDistances[i, :len(sortedIndices)] = distanceMatrix[i, sortedIndices]
    
    output = {}
    output["tgt"] = tgt
    output["nCobrasForTarget"] = nCobrasForTarget
    output["nTargetsForCobra"] = nTargetsForCobra 
    output["targetIndices"] = targetIndices 
    output["targetDistances"] = targetDistances 

    return output


def plotTargets(targets, indices):
    """Plots the target positions.

    Parameters
    ----------
    targets: Object
        A numpy complex array with the target positions. 
    
    """
    if targets is not None:
        plt.figure("Target positions", facecolor="white")
        plt.scatter(np.real(targets[indices]), np.imag(targets[indices]), s=8, color="red")
        plt.scatter(np.real(targets), np.imag(targets), s=2)
        plt.xlabel("x position")
        plt.ylabel("y position")
        plt.title("Target positions")
        plt.show(block=False)


if __name__ == "__main__":
    # Assign the target positions
    out = assignTargets(1, None, False)
    
    # Plot the target positions
    plotTargets(out["tgt"], out["targetIndices"])
    plt.show()

