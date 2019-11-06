"""

Collection of unit tests for the TargetGroup class.

"""

import numpy as np

from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps.cobraConstants import (NULL_TARGET_ID,
                                         NULL_TARGET_POSITION,
                                         NULL_TARGET_PRIORITY)


class TestTargetGroup():
    """A collection of tests for the TargetGroup class.

    """

    def test_constructor(self, targetPositions, targetIds, targetPriorities):
        # Create a new TargetGroup instance using only the target positions
        targets = TargetGroup(targetPositions)

        # Check that the number of targets and the positions are correct
        assert targets.nTargets == len(targetPositions)
        assert np.all(targets.positions == targetPositions)

        # Check that the default ids and priorities are correct
        assert np.all(targets.ids == ["0", "1", "2", "3"])
        assert np.all(targets.priorities == [1, 1, 1, 1])

        # Create a new TargetGroup instance with all the inputs
        targets = TargetGroup(targetPositions, targetIds, targetPriorities)

        # Check that the ids, positions and priorities are correct
        notNull = targetIds != NULL_TARGET_ID
        assert np.all(targets.ids[notNull] == targetIds[notNull])
        assert np.all(targets.positions[notNull] == targetPositions[notNull])
        assert np.all(targets.priorities[notNull] == targetPriorities[notNull])
        assert np.all(targets.ids[~notNull] == NULL_TARGET_ID)
        assert np.all(targets.positions[~notNull] == NULL_TARGET_POSITION)
        assert np.all(targets.priorities[~notNull] == NULL_TARGET_PRIORITY)
        print(targets)
