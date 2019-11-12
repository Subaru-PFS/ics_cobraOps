"""

Collection of unit tests for the TargetSelector abstract class.

"""

import pytest

from ics.cobraOps.TargetSelector import TargetSelector
from ics.cobraOps import targetUtils


class TargetSelectorSubclass(TargetSelector):
    """Dummy TargetSelector subclass used only for tests.

    """

    def run(self):
        return

    def selectTargets(self):
        return


class TestTargetSelector():
    """A collection of tests for the TargetSelector abstract class.

    """

    def test_constructor_exception(self, bench, targets):
        # Check that the constructor raises a exception because we are trying
        # to instantiate an abstract class
        with pytest.raises(TypeError):
            TargetSelector(bench, targets)

    def test_subclass_constructor(self, bench, targets):
        # Check that we don't get an exception if we subclass the abstract
        # class and implement the abstract methods
        TargetSelectorSubclass(bench, targets)

    def test_calculateAccessibleTargets_method(self, bench):
        targets = targetUtils.generateRandomTargets(20, bench)
        selector = TargetSelectorSubclass(bench, targets)
        selector.calculateAccessibleTargets()
