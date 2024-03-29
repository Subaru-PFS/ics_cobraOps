"""

BlackDotGroup class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from . import plotUtils
from .AttributePrinter import AttributePrinter
from .cobraConstants import BLACK_DOT_RADIUS
from .cobraConstants import BLACK_DOT_RELATIVE_POSITION


class BlackDotGroup(AttributePrinter):
    """Class describing the properties of a group of black dots.

    """

    def __init__(self, cobraCenters):
        """Constructs a new BlackDotGroup instance with default properties.

        Parameters
        ----------
        cobraCenters: object
            A complex numpy array with the cobras central positions.

        Returns
        -------
        object
            The BlackDotGroup instance.

        """
        # Set the number of black dots, their central positions and radius
        self.nBlackDots = len(cobraCenters)
        self.centers = cobraCenters + BLACK_DOT_RELATIVE_POSITION
        self.radius = np.full(self.nBlackDots, BLACK_DOT_RADIUS)

    def useCalibrationProduct(self, calibrationProduct):
        """Updates the black dot properties with the calibration product ones.

        Parameters
        ----------
        calibrationProduct: object
            The black dots calibration product containing the black dots
            properties.

        """
        # Copy the main information from the calibration product
        self.nBlackDots = calibrationProduct.nBlackDots
        self.centers = calibrationProduct.centers.copy()
        self.radius = calibrationProduct.radius.copy()

    def addToFigure(self, colors=np.array([0.0, 0.0, 0.0, 0.15]), indices=None):
        """Draws the black dots on top of an existing figure.

        Parameters
        ----------
        colors: object, optional
            The black dots colors. Default is very dark grey.
        indices: object, optional
            A numpy array with the black dots indices to use. If it is set to
            None, all the black dots will be used. Default is None.

        """
        # Extract some useful information
        centers = self.centers
        radius = self.radius

        # Select a subset of the black dots if necessary
        if indices is not None:
            centers = centers[indices]
            radius = radius[indices]

            if colors.ndim == 2:
                colors = colors[indices]

        # Draw the black dots
        plotUtils.addCircles(centers, radius, facecolors=colors)
