
import numpy as np


def getCentersRails(numOfRails):
    # Set the number of cobras per rail
    cobrasPerRail = 29 + 28
    
    # Create the centers array
    centers = np.zeros(numOfRails * cobrasPerRail, dtype="complex")
    
    # Fill the first rail
    firstRail = centers[:cobrasPerRail]
    firstRail[:29] = 8 * np.arange(29)
    firstRail[29:] = 8 * np.arange(28) + 8 * np.exp(1j * np.pi / 3)
    firstRail += 8 * np.exp(1j * 2 * np.pi / 3)
    
    # Order the first rail center by the x coordinate
    orderedIndexes = firstRail.sort()
    firstRail[:] = firstRail[orderedIndexes]
    
    # Fill the rest of the rails
    for i in xrange(1, numOfRails):
        centers[i * cobrasPerRail:(i + 1) * cobrasPerRail] = firstRail + i * 16 * np.exp(1j * 2 * np.pi / 3)
    
    return centers
