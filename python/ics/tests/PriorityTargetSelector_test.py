"""

Collection of unit tests for the PriorityTargetSelector class.

"""

import numpy as np

from ics.cobraOps.PriorityTargetSelector import PriorityTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector
from ics.cobraOps.cobraConstants import (NULL_TARGET_INDEX,
                                         NULL_TARGET_PRIORITY)


class TestPriorityTargetSelector():
    """A collection of tests for the PriorityTargetSelector class.

    """

    def test_orderAccessibleTargetsByPriority_method(self, bench, targets):
        # Create a PriorityTargetSelector
        selector = PriorityTargetSelector(bench, targets)

        # Calculate the accessible targets for each cobra
        selector.calculateAccessibleTargets()

        # Order the accessible targets by priority
        selector.orderAccessibleTargetsByPriority()

        # Check that the accessible targets are not ordered by distance, but
        # are ordered by priority
        orderedByDistance = True
        orderedByPriority = True

        for i in range(bench.cobras.nCobras):
            # Get the accessible target indices, distances and priorities
            indices = selector.accessibleTargetIndices[i]
            distances = selector.accessibleTargetDistances[i]
            priorities = selector.accessibleTargetPriorities[i]
            validTargets = indices != NULL_TARGET_INDEX

            # Check if the distances and the priorities are ordered
            nTargets = np.sum(validTargets)
            orderedByDistance = orderedByDistance and np.all(
                np.sort(distances[:nTargets]) == distances[:nTargets])
            orderedByPriority = orderedByPriority and np.all(
                np.sort(priorities[:nTargets])[::-1] == priorities[:nTargets])

        assert not orderedByDistance
        assert orderedByPriority

    def test_run_method(self, bench, targets):
        # Create a PriorityTargetSelector
        selector = PriorityTargetSelector(bench, targets)

        # Run the complete target selection
        selector.run()

        # Get the selected targets
        selectedTargets = selector.getSelectedTargets()

        # Check that the result makes sense
        assert selectedTargets.nTargets == bench.cobras.nCobras
        assert np.sum(selectedTargets.notNull) > 0

    def test_compare_with_random_selection(self, bench, targets):
        # Add some random priorities
        targets.priorities = 1 + 9 * np.random.random(targets.nTargets)
        targets.priorities[~targets.notNull] = NULL_TARGET_PRIORITY

        # Create a PriorityTargetSelector
        selector = PriorityTargetSelector(bench, targets)

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

        # Check that the average of the target priorities is smaller for the
        # random selector method
        priorities = selectedTargets.priorities[selectedTargets.notNull]
        randomPriorities = randomSelectedTargets.priorities[
            randomSelectedTargets.notNull]
        assert np.mean(priorities) > np.mean(randomPriorities)
