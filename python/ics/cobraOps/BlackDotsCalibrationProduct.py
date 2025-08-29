"""

BlackDotsCalibrationProduct class.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import csv
import numpy as np

from .AttributePrinter import AttributePrinter


class BlackDotsCalibrationProduct(AttributePrinter):
    """Class describing a black dots calibration product.

    """

    def __init__(self, fileName):
        """Constructs a new BlackDotsCalibrationProduct instance using the
        information contained in a csv file.

        Parameters
        ----------
        fileName: object
            The complete path to the csv file.

        Returns
        -------
        object
            The BlackDotsCalibrationProduct instance.

        """
        # Load the csv calibration file
        with open(fileName, newline="") as csvFile:
            blackDotsData = list(csv.reader(csvFile, delimiter=","))

        # Remove the first line because it contains the column names
        blackDotsData = blackDotsData[1:]

        # Get the total number of black dots
        self.nBlackDots = len(blackDotsData)

        # Create the calibration data arrays
        self.id = np.empty(self.nBlackDots, dtype="int")
        self.centers = np.empty(self.nBlackDots, dtype="complex")
        self.radius = np.empty(self.nBlackDots)

        # Fill the arrays
        for i, blackDotData in enumerate(blackDotsData):
            self.id[i] = int(blackDotData[0])
            self.centers[i] = complex(
                float(blackDotData[1]), float(blackDotData[2]))
            self.radius[i] = float(blackDotData[3])
