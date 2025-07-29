"""

Collection of unit tests for the RandomTargetSelector class.

"""

import numpy as np

from ics.cobraOps.RandomTargetSelector import RandomTargetSelector
from ics.cobraOps.TargetGroup import TargetGroup


class TestRandomTargetSelector():
    """A collection of tests for the RandomTargetSelector class.

    """

    def test_orderAccessibleTargetsRandomly_method(self, bench, targets):
        # Create a RandomTargetSelector
        selector = RandomTargetSelector(bench, targets)

        # Calculate the accessible targets for each cobra
        selector.calculateAccessibleTargets()

        # Order the accessible targets randomly
        selector.orderAccessibleTargetsRandomly()

        # Check that the accessible targets are not ordered by distance
        orderedByDistance = True

        for i in range(bench.cobras.nCobras):
            # Get the accessible target indices and distances
            indices = selector.accessibleTargetIndices[i]
            distances = selector.accessibleTargetDistances[i]
            validTargets = indices != TargetGroup.NULL_TARGET_INDEX

            # Check if the distances are ordered
            nTargets = np.sum(validTargets)
            orderedByDistance = orderedByDistance and np.all(
                np.sort(distances[:nTargets]) == distances[:nTargets])

        assert not orderedByDistance

    def test_run_method(self, bench, targets):
        # Create a RandomTargetSelector
        selector = RandomTargetSelector(bench, targets)

        # Run the complete target selection
        selector.run()

        # Get the selected targets
        selectedTargets = selector.getSelectedTargets()

        # Check that the result makes sense
        assert selectedTargets.nTargets == bench.cobras.nCobras
        assert np.sum(selectedTargets.notNull) > 0
