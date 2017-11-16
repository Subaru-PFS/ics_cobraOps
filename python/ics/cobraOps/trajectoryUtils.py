"""

Some utility methods related with the cobra trajectories.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

import ics.cobraOps.plotUtils as plotUtils


def animateCobraTrajecty(centralCobraIndex, trajectories, cobraColors=[0.0, 0.0, 1.0, 0.5], fileName=False):
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
    cobraCenters = trajectories.bench.cobras.centers
    rMax = trajectories.bench.cobras.rMax
    linkRadius = trajectories.bench.cobras.linkRadius
    
    # Get the bench precalculated nearest neighbors information 
    cobras = trajectories.bench.nearestNeighbors[0]
    nearbyCobras = trajectories.bench.nearestNeighbors[1]

    # Get the indices of the cobras that are near the central cobra
    nearCentralCobra = trajectories.bench.cobras.getCobraNeighbors(centralCobraIndex)
 
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
    trajectories.bench.cobras.addLinksToFigure(trajectories.fiberPositions[:, -1], colors=cobraColors, indices=np.logical_not(toAnimate))

    # Limit the relevant arrays to those cobras
    cobraCenters = cobraCenters[toAnimate]
    elbowPositions = trajectories.elbowPositions[toAnimate]
    fiberPositions = trajectories.fiberPositions[toAnimate]
    linkRadius = linkRadius[toAnimate]

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
        lineCollection = plotUtils.addLines(cobraCenters, elbowPositions[:, frame], edgecolor=cobraColors, linewidths=2)
        thickLineCollection = plotUtils.addThickLines(elbowPositions[:, frame], fiberPositions[:, frame], linkRadius, facecolors=cobraColors)  
    
        # Plot also their line trajectories
        trajectoryCollection = plotUtils.addTrajectories(np.vstack((elbowPositions[:, :frame + 1], fiberPositions[:, :frame + 1])), color="0.4", linewidth=1)

        # Log some animation information
        animationPercentage = int(100 * (frame + 1) / elbowPositions.shape[1])
        
        if animationPercentage % 20 == 0:
            print("Trajectory percentage: {}%".format(animationPercentage))

    # Add the animation to the current figure
    plotUtils.addAnimation(update, elbowPositions.shape[1], fileName=fileName)


def animateTrajectories(trajectories, cobraColors=[0.0, 0.0, 1.0, 0.5]):
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
    cobraCenters = trajectories.bench.cobras.centers
    linkRadius = trajectories.bench.cobras.linkRadius

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
        lineCollection = plotUtils.addLines(cobraCenters, trajectories.elbowPositions[:, frame], edgecolor=cobraColors, linewidths=2)
        thickLineCollection = plotUtils.addThickLines(trajectories.elbowPositions[:, frame], trajectories.fiberPositions[:, frame], linkRadius, facecolors=cobraColors)  
    
        # Plot also their line trajectories
        trajectoryCollection = plotUtils.addTrajectories(np.vstack((trajectories.elbowPositions[:, :frame + 1], trajectories.fiberPositions[:, :frame + 1])), color="0.4", linewidth=1)

        # Log some animation information
        animationPercentage = int(100 * (frame + 1) / trajectories.elbowPositions.shape[1])
        
        if animationPercentage % 20 == 0:
            print("Trajectory percentage: {}%".format(animationPercentage))

    # Add the animation to the current figure
    plotUtils.addAnimation(update, trajectories.elbowPositions.shape[1], fileName=None)
