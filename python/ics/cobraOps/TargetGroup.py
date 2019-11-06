"""

TargetGroup class.

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
from .cobraConstants import (NULL_TARGET_INDEX,
                             NULL_TARGET_POSITION,
                             NULL_TARGET_ID)


class TargetGroup(AttributePrinter):
    """Class describing the properties of a group of targets.

    """

    def __init__(self, positions, ids=None, priorities=None):
        """Constructs a new target group instance.

        Parameters
        ----------
        positions: object
            A complex numpy array with the targets positions.
        ids: object, optional
            A unicode numpy array with the targets unique ids. If it is set to
            None, the targets will have consecutive ids starting from zero.
            Default is None.
        priorities: object, optional
            A numpy array with the targets priorities. If it is set to None,
            all the targets will have priority 1. Default is None.

        Returns
        -------
        object
            The target group instance.

        """
        # Save the number of targets and their positions
        self.nTargets = len(positions)
        self.positions = positions.copy()

        # Save the target ids or create some default ids
        if ids is not None:
            self.ids = ids.copy()
        else:
            self.ids = np.arange(self.nTargets).astype("<U10")

        # Save the target priorities or set them to 1
        if priorities is not None:
            self.priorities = priorities.copy()
        else:
            self.priorities = np.ones(self.nTargets)

        # Check which targets are not NULL
        self.notNull = self.ids != NULL_TARGET_ID

    @classmethod
    def fromFile(cls, fileName):
        """Constructs a new target group instance from an input file.

        Parameters
        ----------
        fileName: object
            The path to the file with the target group information.

        Returns
        -------
        object
            The target group instance.

        """
        # Read the input file and save the targets information in two lists
        positions = []
        ids = []

        with open(fileName, "r") as f:
            for line in f:
                # The data should be separated with comas
                data = line.split(",")
                positions.append(float(data[0]) + 1j * float(data[1]))

                # Check if the file contains the target ids
                if len(data) > 2:
                    ids.append(data[2].strip())

        # Transform the lists into numpy arrays
        positions = np.array(positions, dtype="complex")
        ids = np.array(ids, dtype="<U10")

        # Return a new target group instance
        if len(ids) == len(positions):
            return cls(positions, ids)
        else:
            return cls(positions)

    def saveToFile(self, fileName):
        """Saves the target group data to a file.

        Parameters
        ----------
        fileName: object
            The path to the output file where the target group information
            should be saved.

        """
        with open(fileName, "w") as f:
            for p, i in zip(self.positions, self.ids):
                f.write(", ".join((str(p.real), str(p.imag), i)) + "\n")

    def select(self, indices):
        """Selects a subset of the targets.

        Parameters
        ----------
        indices: object
            An numpy array with the indices of the targets to select.

        Returns
        -------
        object
            A new target group instance containing only the selected targets.

        """
        # Initialize the selected target positions and ids arrays
        selectedPositions = np.full(len(indices), NULL_TARGET_POSITION, dtype="complex")
        selectedIds = np.full(len(indices), NULL_TARGET_ID, dtype=self.ids.dtype)

        # Fill the arrays with the targets that are not NULL
        realIndices = indices != NULL_TARGET_INDEX
        selectedPositions[realIndices] = self.positions[indices[realIndices]]
        selectedIds[realIndices] = self.ids[indices[realIndices]]

        return TargetGroup(selectedPositions, selectedIds)

    def addToFigure(self, colors=np.array([0.4, 0.4, 0.4, 1.0]), indices=None):
        """Draws the targets on top of an existing figure.

        Parameters
        ----------
        colors: object, optional
            The target colors. Default is dark grey.
        indices: object, optional
            A numpy array with the target indices to use. If it is set to None,
            all the targets will be used. Default is None.

        """
        # Extract some useful information
        positions = self.positions
        notNull = self.notNull

        # Select a subset of the targets if necessary
        if indices is not None:
            positions = positions[indices]
            notNull = notNull[indices]

            if colors.ndim == 2:
                colors = colors[indices]

        # Remove the NULL targets
        positions = positions[notNull]

        if colors.ndim == 2:
            colors = colors[notNull]

        # Draw the targets
        plotUtils.addPoints(positions, s=2, facecolor=colors)
