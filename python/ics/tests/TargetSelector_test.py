"""

Collection of unit tests for the TargetSelector abstract class.

"""

import pytest
import numpy as np

from ics.cobraOps.TargetSelector import TargetSelector
from ics.cobraOps.Bench import Bench
from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps import targetUtils
from ics.cobraOps.cobraConstants import NULL_TARGET_INDEX


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

    def test_constructKDTree_method(self):
        # Create a basic bench with 2 cobras
        cobraCenters = np.array([0 + 0j, 0 + 5j])
        bench = Bench(cobraCenters)

        # Create some targets
        targetPositions = np.array([0, cobraCenters.mean(), 0 - 3j, 0 + 6j])
        targets = TargetGroup(targetPositions)

        # Create a dummy TargetSelector
        selector = TargetSelectorSubclass(bench, targets)

        # Construct the KD tree
        selector.constructKDTree()

        # Check that the KD tree is available
        assert selector.kdTree is not None

        # Check that it contains the expected data
        assert np.all(selector.kdTree.data == np.column_stack(
            (targetPositions.real, targetPositions.imag)))

    def test_getTargetsInsidePatrolArea_method(self):
        # Create a basic bench with 2 cobras
        cobraCenters = np.array([0 + 0j, 0 + 5j])
        bench = Bench(cobraCenters)

        # Create some targets
        targetPositions = np.array([0, cobraCenters.mean(), 0 - 3j, 0 + 6j])
        targets = TargetGroup(targetPositions)

        # Create a dummy TargetSelector
        selector = TargetSelectorSubclass(bench, targets)

        # Get the targets that fall inside the first cobra
        indices, positions, distances = selector.getTargetsInsidePatrolArea(0)

        # Check that we get what we expected
        expectedIndices = np.array([1, 2])
        expectedPositions = targetPositions[expectedIndices]
        expectedDistances = np.abs(cobraCenters[0] - expectedPositions)
        assert np.all(indices == expectedIndices)
        assert np.all(positions == expectedPositions)
        assert np.all(np.abs(distances - expectedDistances) < 1e-10)

        # Get the targets that fall inside the second cobra
        indices, positions, distances = selector.getTargetsInsidePatrolArea(1)

        # Check that we get what we expected
        expectedIndices = np.array([3, 1])
        expectedPositions = targetPositions[expectedIndices]
        expectedDistances = np.abs(cobraCenters[1] - expectedPositions)
        assert np.all(indices == expectedIndices)
        assert np.all(positions == expectedPositions)
        assert np.all(np.abs(distances - expectedDistances) < 1e-10)

    def test_calculateAccessibleTargets_method(self, bench):
        # Create a random set of targets
        targets = targetUtils.generateRandomTargets(0.1, bench)
        selector = TargetSelectorSubclass(bench, targets)
        selector.calculateAccessibleTargets()

        # Check that the targets fall inside the cobras patrol areas
        cobraCenters = bench.cobras.centers
        targetIndices = selector.accessibleTargetIndices
        targetPositions = targets.positions[targetIndices]
        validTargets = targetIndices != NULL_TARGET_INDEX
        distances = np.abs(cobraCenters[:, np.newaxis] - targetPositions)
        assert np.all(
            (distances < bench.cobras.rMax[:, np.newaxis])[validTargets])
        
