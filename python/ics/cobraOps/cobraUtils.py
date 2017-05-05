"""

Some utility methods related with the PFS cobras.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T
  
"""

import numpy as np
import matplotlib.pyplot as plt


COBRAS_SEPARATION = 8.0
"""The separation in mm between two consecutive cobras."""

MODULE_FIRST_LINE_LENGTH = 29
"""The number of cobras in the first line of a module."""

MODULE_SECOND_LINE_LENGTH = 28
"""The number of cobras in the second line of a module."""

MODULES_PER_SECTOR = 14
"""The number of modules in one PFI sector."""


def getCobrasCenters(cobraLayout):
    """Calculates the cobras central positions for a given cobra layout.

    Parameters
    ----------
    cobraLayout: str
        The cobras layout to use: none, hex, line, rails or full. 
        "none" means that the cobras central positions will be extracted 
        from a configuration file.
    
    Returns
    -------
    Object
        Complex numpy array with the cobras central positions.
    
    """
    if cobraLayout == "none":
        # Nothing to do. A configuration file will be used.
        centers = None
    elif cobraLayout == "hex":
        # The centers should follow an hexagon pattern.
        # (R3 collisions = 3.0e-4 w/ 100k targets)
        centers = np.zeros(7, dtype="complex")
        centers[1:] = COBRAS_SEPARATION * np.exp(np.arange(6) * 1j * np.pi / 3)
    elif cobraLayout == "line":
        # Line of 27 cobras. 
        # (R3 collisions = 4.5e-5 (fraction) w/ 1M targets)
        centers = COBRAS_SEPARATION * np.arange(27) + 1j
    elif cobraLayout == "rails":
        # One PFI sector with 14 cobra modules.
        # (R3 collisions = 9.9e-5 w/ 16k targets)
        centers = getFirstSectorCenters()
    elif cobraLayout == "full":
        # Full PFI bench (2394 cobras distributed in 3 rotated sectors). 
        # (R3 collisions = 1.7e-4 (fraction) w/ 1M targets)
        # Get the first sector
        firstSector = getFirstSectorCenters()
         
        # Create the cobras centers array
        cobrasPerSector = len(firstSector)
        centers = np.zeros(3 * cobrasPerSector, dtype="complex")
            
        # Add the first sector
        centers[:cobrasPerSector] = firstSector 
           
        # Add the second sector rotating the first sector 120 degrees
        centers[cobrasPerSector:-cobrasPerSector] = firstSector * np.exp(1j * 2 * np.pi / 3) 

        # Add the third sector rotating the first sector 240 degrees
        centers[-cobrasPerSector:] = firstSector * np.exp(1j * 4 * np.pi / 3)
    else:
        raise Exception(cobraLayout + " is not a valid cobra layout: none, hex, line, rails or full.")

    return centers


def getFirstSectorCenters():
    """Calculates the cobras central positions for the first PFI sector.
    
    Returns
    -------
    Object
        Complex numpy array with the cobras central positions.
    
    """
    # Create the cobras centers array
    cobrasPerModule = MODULE_FIRST_LINE_LENGTH + MODULE_SECOND_LINE_LENGTH
    centers = np.zeros(cobrasPerModule * MODULES_PER_SECTOR, dtype="complex")
    
    # Fill the first module
    firstModule = centers[:cobrasPerModule]
    firstModule[:MODULE_FIRST_LINE_LENGTH] = COBRAS_SEPARATION * np.arange(MODULE_FIRST_LINE_LENGTH)
    firstModule[MODULE_FIRST_LINE_LENGTH:] = COBRAS_SEPARATION * np.arange(MODULE_SECOND_LINE_LENGTH) + \
                                             COBRAS_SEPARATION * np.exp(1j * np.pi / 3)
    firstModule += COBRAS_SEPARATION * np.exp(1j * 2 * np.pi / 3)
    
    # Order the first module centers by the x coordinate
    orderedIndexes = firstModule.sort()
    firstModule[:] = firstModule[orderedIndexes]
    
    # Fill the rest of the modules
    modulesOffset = 2 * COBRAS_SEPARATION * np.exp(1j * 2 * np.pi / 3)
     
    for i in xrange(1, MODULES_PER_SECTOR):
        centers[i * cobrasPerModule:(i + 1) * cobrasPerModule] = firstModule + i * modulesOffset
    
    return centers


def plotCobrasCenters(centers):
    """Plot the cobras central positions.

    Parameters
    ----------
    centers: Object
        A numpy complex array with the cobras central positions. 
    
    """
    if centers is not None:
        plt.figure("Cobra centers", facecolor="white")
        plt.scatter(np.real(centers), np.imag(centers), s=2)
        plt.xlabel("x position")
        plt.ylabel("y position")
        plt.title("Cobra centers")
        plt.show(block=False)


if __name__ == "__main__":
    # Get the centers for the full PFI
    centers = getCobrasCenters("full")
    
    # Plot the centers
    plotCobrasCenters(centers)
    plt.show()
    