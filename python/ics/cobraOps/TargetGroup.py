"""

Target group class.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

import ics.cobraOps.plotUtils as plotUtils

from ics.cobraOps.AttributePrinter import AttributePrinter
from ics.cobraOps.cobraConstants import (NULL_TARGET_INDEX,
                                         NULL_TARGET_POSITION,
                                         NULL_TARGET_ID)


class TargetGroup(AttributePrinter):
    """
    
    Class describing the properties of a group of targets.
    
    """
    
    def __init__(self, positions, ids=None):
        """Constructs a new target group instance.
        
        Parameters
        ----------
        positions: object
            A complex numpy array with the targets positions.
        ids: object, optional
            A unicode numpy array with the targets unique ids. If it is set to
            None, the targets will have consecutive ids starting from zero.
            Default is None.
        
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
        
        # Check which targets are not NULL
        self.notNull = self.ids != NULL_TARGET_ID
    
    
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
    
    
    def addToFigure(self, colors=[0.4, 0.4, 0.4, 1.0], indices=None):
        """Draws the targets on top off an existing figure.
        
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
