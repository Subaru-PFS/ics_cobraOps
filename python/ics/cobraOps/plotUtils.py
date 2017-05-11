"""

Some utility methods to make plots using matplotlib.
  
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as collections


def createNewFigure(title, xLabel, yLabel, size=(7, 7), **kwargs):
    """Initializes a new matplotlib figure.
    
    Parameters
    ----------
    title: str
        The plot title.
    xLabel: str
        The label for the x axis.
    yLabel: str
        The label for the y axis.
    size: tuple, optional
        The figure size. Default is (7, 7).
    kwargs: plt.figure properties
        Any additional property that should be passed to the figure.
    
    """
    # Create the figure
    plt.figure(title, figsize=size, facecolor="white", tight_layout=True, **kwargs)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.show(block=False)

    # Fix the axes aspect ratio
    ax = plt.gca()
    ax.set_aspect("equal")


def setAxesLimits(xLim, yLim):
    """Sets the axes limits of an already initialized figure.
    
    Parameters
    ----------
    xLim: object
        A numpy array with the x axis limits.
    yLim: object
        A numpy array with the y axis limits.
    
    """
    ax = plt.gca()
    ax.set_xlim(xLim)
    ax.set_ylim(yLim)


def addPoints(points, **kwargs):
    """Adds an array of complex points to an already initialized figure.
    
    Parameters
    ----------
    points: object
        Complex numpy array with the points coordinates.
    kwargs: plt.scatter properties
        Any additional property that should be passed to the scatter function.
    
    """
    plt.scatter(np.real(points), np.imag(points), **kwargs)


def addCircles(centers, radii, **kwargs):
    """Adds a set of circles to an already initialized figure.

    Parameters
    ----------
    centers: object
        Complex numpy array with the circle central coordinates. 
    radii: object
        Numpy array with the circle radii. 
    kwargs: collections.PatchCollection properties
        Any additional property that should be passed to the patch collection.
         
    """
    # Create the circle collection
    circleList = [patches.Circle((c.real, c.imag), r) for c, r in zip(centers, radii)]
    circleCollection = collections.PatchCollection(circleList, **kwargs)

    # Plot the circles in the current figure
    plt.gca().add_collection(circleCollection)


def addLines(startPoints, endPoints, **kwargs):
    """Adds a set of lines to an already initialized figure.

    Parameters
    ----------
    startPoints: object
        Complex numpy array with the lines start coordinates. 
    endPoints: object
        Complex numpy array with the lines end coordinates. 
    kwargs: collections.LineCollection properties
        Any additional property that should be passed to the line collection.
         
    """
    # Create the line collection
    lineList = [[(p1.real, p1.imag), (p2.real, p2.imag)] for p1, p2 in zip(startPoints, endPoints)]
    lineCollection = collections.LineCollection(lineList, **kwargs)

    # Plot the line in the current figure
    plt.gca().add_collection(lineCollection)


def pauseExecution():
    """Pauses the general program execution to allow figure inspection.
    
    """
    plt.show()
