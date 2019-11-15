"""

Collection of unit tests for the RandomTargetSelector class.

"""

from ics.cobraOps.RandomTargetSelector import RandomTargetSelector


class TestRandomTargetSelector():
    """A collection of tests for the RandomTargetSelector class.

    """

    def test_constructor(self, bench, targets):
        RandomTargetSelector(bench, targets)
