"""

Some utility methods to make plots using matplotlib.
  
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.collections as collections


def createNewFigure(title, xLabel, yLabel, size=(8, 8), **kwargs):
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
        The figure size. Default is (10, 10).
    kwargs: plt.figure properties
        Any additional property that should be passed to the figure.
    
    """
    # Create the figure
    plt.figure(figsize=size, facecolor="white", tight_layout=True, **kwargs)
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
    plt.scatter(points.real, points.imag, **kwargs)


def addCircles(centers, radii, **kwargs):
    """Adds a set of circles to an already initialized figure.

    Parameters
    ----------
    centers: object
        Complex numpy array with the circle central coordinates. 
    radii: object
        Numpy array with the circle radii. 
    kwargs: collections.EllipseCollection properties
        Any additional property that should be passed to the ellipse collection.
         
    """
    # Create the ellipse collection
    diameters = 2 * radii
    offsets = np.hstack((centers.real[:, np.newaxis], centers.imag[:, np.newaxis]))
    angles = np.zeros(len(centers))
    params = {"offsets":offsets, "units":"xy", "transOffset":plt.gca().transData}
    params.update(kwargs)
    ellipseCollection = collections.EllipseCollection(diameters, diameters, angles, **params)
 
    # Plot the ellipses in the current figure
    plt.gca().add_collection(ellipseCollection)


def addRings(centers, innerRadii, outerRadii, **kwargs):
    """Adds a set of rings to an already initialized figure.

    Parameters
    ----------
    centers: object
        Complex numpy array with the circle central coordinates. 
    innerRadii: object
        Numpy array with the ring inner radii. 
    outerRadii: object
        Numpy array with the ring outer radii. 
    kwargs: collections.EllipseCollection properties
        Any additional property that should be passed to the ellipse collection.
         
    """
    # Create the rings using 2 sets of circles
    addCircles(centers, outerRadii, **kwargs)
    addCircles(centers, innerRadii, facecolors="white")


def addRingsSlow(centers, innerRadii, outerRadii, **kwargs):
    """Adds a set of rings to an already initialized figure.

    Parameters
    ----------
    centers: object
        Complex numpy array with the circle central coordinates. 
    innerRadii: object
        Numpy array with the ring inner radii. 
    outerRadii: object
        Numpy array with the ring outer radii. 
    kwargs: collections.PatchCollection properties
        Any additional property that should be passed to the patch collection.
         
    """
    # Create the ring collection
    ringList = [patches.Wedge((c.real, c.imag), r1, 0, 360, r1 - r2) for c, r1, r2 in zip(centers, outerRadii, innerRadii)]
    ringCollection = collections.PatchCollection(ringList, **kwargs)

    # Plot the rings in the current figure
    plt.gca().add_collection(ringCollection)


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


def addThickLines(startPoints, endPoints, thicknesses, **kwargs):
    """Adds a set of thick lines to an already initialized figure.

    Parameters
    ----------
    startPoints: object
        Complex numpy array with the lines start coordinates. 
    endPoints: object
        Complex numpy array with the lines end coordinates. 
    thicknesses: object
        Numpy array with the lines thicknesses.
    kwargs: collections.LineCollection properties
        Any additional property that should be passed to the line collection.
         
    """
    # Create the thick line collection
    centers = (endPoints + startPoints) / 2
    difference = endPoints - startPoints
    widths = np.abs(difference)
    angles = np.angle(difference)
    thickLineList = [getThickLinePath(c, w, t, a) for c, w, t, a in zip(centers, widths, thicknesses, angles)]
    thickLineCollection = collections.PatchCollection(thickLineList, **kwargs)

    # Plot the thick lines in the current figure
    plt.gca().add_collection(thickLineCollection)


def getThickLinePath(center, width, thickness, angle):
    """Creates a matplotlib path representing a thick line.

    Parameters
    ----------
    center: complex
        Complex number with the line center coordinate. 
    width: float
        The line width.
    thicknesses: float
        The line thickness. The line height will be twice the thickness.
    angle: float
        The line rotation angle.
         
    """
    # Create a circle path with the appropriate radius centered at the origin
    circlePath = path.Path.circle((0, 0), thickness)
    
    # Extract the vertices and the spline codes
    circleVerts = circlePath.vertices
    circleCodes = circlePath.codes
    circlePoints = circlePath.vertices.shape[0]
    middleIndex = int(circlePoints / 2)

    # Create the thick line vertices array reusing the circle ones
    vertices = np.empty((circlePoints + 1, 2), dtype=circleVerts.dtype)
    vertices[:middleIndex] = circleVerts[:middleIndex] + [0.5 * width, 0.0]
    vertices[middleIndex] = [-0.5 * width, thickness]
    vertices[-middleIndex:] = circleVerts[-middleIndex:] - [0.5 * width, 0.0]

    # Rotate and translate the vertices to the given center and angle
    cos = np.cos(angle)
    sin = np.sin(angle)
    x = cos * vertices[:, 0] - sin * vertices[:, 1] + center.real
    vertices[:, 1] = sin * vertices[:, 0] + cos * vertices[:, 1] + center.imag
    vertices[:, 0] = x
    
    # Create the thick line codes array reusing the circle ones
    codes = np.empty(circlePoints + 1, dtype=circleCodes.dtype)
    codes[:middleIndex] = circleCodes[:middleIndex]
    codes[middleIndex] = path.Path.LINETO
    codes[-middleIndex:] = circleCodes[-middleIndex:]

    # Return the thick line path
    return patches.PathPatch(path.Path(vertices, codes))


def pauseExecution():
    """Pauses the general program execution to allow figure inspection.
    
    """
    plt.show()

