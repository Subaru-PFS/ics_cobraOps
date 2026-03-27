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

    def __init__(self, ids, centers, radius):
        """Constructs a new BlackDotsCalibrationProduct instance.

        Parameters
        ----------
        ids: object
            A numpy array with the black dots ids.
        centers: object
            A numpy complex array with the black dots centers.
        radius: object
            A numpy array with with the black dots radii in mm.

        Returns
        -------
        object
            The BlackDotsCalibrationProduct instance.

        """
        self.nBlackDots = len(ids)
        self.ids = ids.astype("int")
        self.centers = centers.astype("complex")
        self.radius = radius.astype("float")

    @classmethod
    def from_file(cls, fileName):
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
        nBlackDots = len(blackDotsData)

        # Create the calibration data arrays
        ids = np.empty(nBlackDots, dtype="int")
        centers = np.empty(nBlackDots, dtype="complex")
        radius = np.empty(nBlackDots)

        # Fill the arrays
        for i, blackDotData in enumerate(blackDotsData):
            ids[i] = int(blackDotData[0])
            centers[i] = complex(float(blackDotData[1]), float(blackDotData[2]))
            radius[i] = float(blackDotData[3])

        return cls(ids, centers, radius)

    @classmethod
    def from_pandas(cls, blackDotsData):
        """Constructs a new BlackDotsCalibrationProduct instance using the
        information contained in a pandas DataFrame.

        Parameters
        ----------
        blackDotsData: object
            A pandas DataFrame with the black dots data. It must contain the
            following columns: 'spotId', 'x', 'y' and 'r'.

        Returns
        -------
        object
            The BlackDotsCalibrationProduct instance.

        """
        ids = blackDotsData["spotId"].values
        centers = blackDotsData["x"].values + 1j * blackDotsData["y"].values
        radius = blackDotsData["r"].values

        return cls(ids, centers, radius)
