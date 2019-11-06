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
                             NULL_TARGET_ID,
                             NULL_TARGET_PRIORITY)


class TargetGroup(AttributePrinter):
    """Class describing the properties of a group of targets.

    """

    def __init__(self, positions, ids=None, priorities=None):
        """Constructs a new TargetGroup instance.

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
            The TargetGroup instance.

        """
        # Save the number of targets and their positions
        self.nTargets = len(positions)
        self.positions = positions.copy()

        # Save the target ids or create some default ids
        if ids is not None:
            # Make sure the array length is correct
            if len(ids) != self.nTargets:
                raise ValueError("The ids array should have the same length "
                                 "as the positions array.")

            self.ids = ids.copy()
        else:
            self.ids = np.arange(self.nTargets).astype("<U10")

        # Save the target priorities or set them to 1
        if priorities is not None:
            # Make sure the array length is correct
            if len(priorities) != self.nTargets:
                raise ValueError("The priorities array should have the same "
                                 "length as the positions array.")

            self.priorities = priorities.copy()
        else:
            self.priorities = np.ones(self.nTargets)

        # Check which targets are not NULL
        self.notNull = self.ids != NULL_TARGET_ID

        # Use the default position and priority values for the NULL targets
        self.positions[~self.notNull] = NULL_TARGET_POSITION
        self.priorities[~self.notNull] = NULL_TARGET_PRIORITY

    @classmethod
    def fromFile(cls, fileName):
        """Constructs a new TargetGroup instance from an input file.

        The data for each target should be separated by commas:
        x, y, id, priority

        Parameters
        ----------
        fileName: object
            The path to the file with the target group information.

        Returns
        -------
        object
            The TargetGroup instance.

        """
        # Read the input file and save the targets information
        positions = []
        ids = []
        priorities = []

        with open(fileName, "r") as f:
            for line in f:
                # The data should be separated with commas
                data = line.split(",")
                positions.append(float(data[0]) + 1j * float(data[1]))

                # Check if the file contains the target ids
                if len(data) > 2:
                    ids.append(data[2].strip())

                # Check if the file contains the target priorities
                if len(data) > 3:
                    priorities.append(float(data[3]))

        # Transform the lists into numpy arrays
        positions = np.array(positions, dtype="complex")
        ids = np.array(ids, dtype="<U10") if len(ids) > 0 else None
        priorities = np.array(priorities) if len(priorities) > 0 else None

        # Return a new TargetGroup instance
        return cls(positions, ids, priorities)

    def saveToFile(self, fileName):
        """Saves the target group data to a file.

        Parameters
        ----------
        fileName: object
            The path to the output file where the target group information
            should be saved.

        """
        with open(fileName, "w+") as f:
            for pos, i, prio in zip(self.positions, self.ids, self.priorities):
                f.write("%s, %s, %s, %s\n" % (pos.real, pos.imag, i, prio))

    def select(self, indices):
        """Selects a subset of the targets.

        Parameters
        ----------
        indices: object
            An numpy array with the indices of the targets to select.

        Returns
        -------
        object
            A new TargetGroup instance containing only the selected targets.

        """
        # Initialize the selected target positions, ids and priorities arrays
        selectedPositions = np.full(
            len(indices), NULL_TARGET_POSITION, dtype=self.positions.dtype)
        selectedIds = np.full(
            len(indices), NULL_TARGET_ID, dtype=self.ids.dtype)
        selectedPriorities = np.full(
            len(indices), NULL_TARGET_PRIORITY, dtype=self.priorities.dtype)

        # Fill the arrays with the targets that are not NULL
        realTargets = indices != NULL_TARGET_INDEX
        realIndices = indices[realTargets]
        selectedPositions[realTargets] = self.positions[realIndices]
        selectedIds[realTargets] = self.ids[realIndices]
        selectedPriorities[realTargets] = self.priorities[realIndices]

        return TargetGroup(selectedPositions, selectedIds, selectedPriorities)

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
        priorities = self.priorities
        notNull = self.notNull

        # Select a subset of the targets if necessary
        if indices is not None:
            indices = indices[indices != NULL_TARGET_INDEX]
            positions = positions[indices]
            priorities = priorities[indices]
            notNull = notNull[indices]

            if colors.ndim == 2:
                colors = colors[indices]

        # Remove the NULL targets
        positions = positions[notNull]
        priorities = priorities[notNull]

        if colors.ndim == 2:
            colors = colors[notNull]

        # Draw the targets
        plotUtils.addPoints(positions, s=(2 * priorities), facecolor=colors)
