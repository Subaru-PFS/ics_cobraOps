"""

Collection of unit tests for the DistanceTargetSelector class.

"""

from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector


class TestDistanceTargetSelector():
    """A collection of tests for the DistanceTargetSelector class.

    """

    def test_constructor(self, bench, targets):
        DistanceTargetSelector(bench, targets)
