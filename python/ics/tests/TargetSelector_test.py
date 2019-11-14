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


@pytest.fixture(scope="function", params=[(True, 0.1, True),
                                          (True, 0.1, False),
                                          (True, 10, True),
                                          (True, 10, False),
                                          (False, 1, True),
                                          (False, 1, False)])
def testParameters(request):
    # Extract the test parameters
    randomTargetDistribution, targetDensity, useKDTree = request.param

    return randomTargetDistribution, targetDensity, useKDTree


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

        # Create some random targets
        targetPositions = 8 * np.random.random(50) + 8j * np.random.random(50)
        targets = TargetGroup(targetPositions)

        # Create a dummy TargetSelector
        selector = TargetSelectorSubclass(bench, targets)

        # Check that the KD tree is not available
        assert selector.kdTree is None

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

        # Construct the KD tree
        selector.constructKDTree()

        # Get again the targets that fall inside the second cobra
        indices, positions, distances = selector.getTargetsInsidePatrolArea(1)

        # Check that we get what we expected
        expectedIndices = np.array([3, 1])
        expectedPositions = targetPositions[expectedIndices]
        expectedDistances = np.abs(cobraCenters[1] - expectedPositions)
        assert np.all(indices == expectedIndices)
        assert np.all(positions == expectedPositions)
        assert np.all(np.abs(distances - expectedDistances) < 1e-10)

    def test_calculateAccessibleTargets_method(self, bench, testParameters):
        # Extract the test parameters
        randomTargetDistribution, targetDensity, useKDTree = testParameters

        # Create the targets
        if randomTargetDistribution:
            targets = targetUtils.generateRandomTargets(targetDensity, bench)
        else:
            targets = targetUtils.generateOneTargetPerCobra(bench)

        # Create a dummy TargetSelector
        selector = TargetSelectorSubclass(bench, targets)

        # Construct the KD tree if requested
        if useKDTree:
            selector.constructKDTree()

        # Calculate the accessible targets for each cobra
        selector.calculateAccessibleTargets()

        # Check that the targets fall inside the cobras patrol areas
        cobraCenters = bench.cobras.centers
        indices = selector.accessibleTargetIndices
        positions = targets.positions[indices]
        distances = np.abs(cobraCenters[:, np.newaxis] - positions)
        validTargets = indices != NULL_TARGET_INDEX
        assert np.all(
            (distances < bench.cobras.rMax[:, np.newaxis])[validTargets])
        assert np.all(
            (distances > bench.cobras.rMin[:, np.newaxis])[validTargets])

        # Check that the accessible targets are ordered by distance
        for i in range(bench.cobras.nCobras):
            indices = selector.accessibleTargetIndices[i]
            distances = selector.accessibleTargetDistances[i]
            validTargets = indices != NULL_TARGET_INDEX
            nTargets = np.sum(validTargets)
            assert np.all(
                distances[:nTargets] == np.sort(distances[validTargets]))
