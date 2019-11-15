"""

Collection of unit tests for the DistanceTargetSelector class.

"""

import numpy as np

from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector


class TestDistanceTargetSelector():
    """A collection of tests for the DistanceTargetSelector class.

    """

    def test_run_method(self, bench, targets):
        # Create a DistanceTargetSelector
        selector = DistanceTargetSelector(bench, targets)

        # Run the complete target selection
        selector.run()

        # Get the selected targets
        selectedTargets = selector.getSelectedTargets()

        # Check that the result makes sense
        assert selectedTargets.nTargets == bench.cobras.nCobras
        assert np.sum(selectedTargets.notNull) > 0

    def test_compare_with_random_selection(self, bench, targets):
        # Create a DistanceTargetSelector
        selector = DistanceTargetSelector(bench, targets)

        # Run the complete target selection
        selector.run()

        # Get the selected targets
        selectedTargets = selector.getSelectedTargets()

        # Create a RandomTargetSelector
        randomSelector = RandomTargetSelector(bench, targets)

        # Run the complete target selection
        randomSelector.run()

        # Get the selected random targets
        randomSelectedTargets = randomSelector.getSelectedTargets()

        # Check that the average of the target distances is larger for the
        # random selector method
        cobraCenters = bench.cobras.centers
        distances = np.abs(cobraCenters - selectedTargets.positions)
        distances = distances[selectedTargets.notNull]
        randomDistances = np.abs(cobraCenters - randomSelectedTargets.positions)
        randomDistances = randomDistances[randomSelectedTargets.notNull]
        assert np.mean(distances) < np.mean(randomDistances)
