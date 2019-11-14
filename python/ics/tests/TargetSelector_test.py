"""

Collection of unit tests for the TargetSelector abstract class.

"""

import pytest
import numpy as np

from ics.cobraOps.TargetSelector import TargetSelector
from ics.cobraOps.Bench import Bench
from ics.cobraOps.TargetGroup import TargetGroup
from ics.cobraOps import targetUtils
from ics.cobraOps.cobraConstants import (NULL_TARGET_INDEX,
                                         NULL_TARGET_POSITION,
                                         NULL_TARGET_ID,
                                         NULL_TARGET_PRIORITY)


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
        cobraCenters = np.array([0, 5], dtype=np.complex)
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
        cobraCenters = np.array([0, 5], dtype=np.complex)
        bench = Bench(cobraCenters)

        # Create some targets
        targetPositions = np.array([0, cobraCenters.mean(), -3, 6])
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

    def test_solveEndPointCollisions_method(self):
        # Create a basic bench with 2 cobras
        cobraCenters = np.array([0, 5], dtype=np.complex)
        bench = Bench(cobraCenters)

        # Create some targets
        rMin = bench.cobras.rMin
        targetPositions = np.array([cobraCenters[0] + 1.2 * rMin[0],
                                    cobraCenters.mean() - 0.1,
                                    cobraCenters.mean(),
                                    cobraCenters[1] + 1.2 * rMin[1]])
        targets = TargetGroup(targetPositions)

        # Create a dummy TargetSelector
        selector = TargetSelectorSubclass(bench, targets)

        # Calculate the accessible targets for each cobra
        selector.calculateAccessibleTargets()

        # Assign the targets in a way that we have an end-point collision
        selector.assignedTargetIndices = np.array([1, 2])

        # Solve the cobra end-point collisions
        selector.solveEndPointCollisions()

        # Check that we get the correct assignment
        assert np.all(selector.assignedTargetIndices == [0, 3])

        # Assign the targets in a way that we don't have an end-point collision
        selector.assignedTargetIndices = np.array([0, 2])

        # Solve the cobra end-point collisions
        selector.solveEndPointCollisions()

        # Check that we get the same assignment
        assert np.all(selector.assignedTargetIndices == [0, 2])

    def test_getSelectedTargets_method(self):
        # Create a basic bench with 2 cobras
        cobraCenters = np.array([0, 5], dtype=np.complex)
        bench = Bench(cobraCenters)

        # Create some targets
        targetPositions = np.array([0, 3, NULL_TARGET_POSITION, 6])
        targets = TargetGroup(targetPositions)

        # Create a dummy TargetSelector
        selector = TargetSelectorSubclass(bench, targets)

        # Assign the cobras to valid targets
        indices = np.array([3, 0])
        selector.assignedTargetIndices = indices

        # Get the selected targets
        selectedTargers = selector.getSelectedTargets()

        # Check that we get the correct targets
        assert np.all(selectedTargers.positions == targetPositions[indices])
        assert np.all(selectedTargers.ids == indices.astype(np.str))
        assert np.all(selectedTargers.priorities == 1)

        # Assign one of the cobras to the NULL target
        indices = np.array([1, 2])
        selector.assignedTargetIndices = indices

        # Get the selected targets
        selectedTargers = selector.getSelectedTargets()

        # Check that we get the NULL target
        assert selectedTargers.positions[1] == NULL_TARGET_POSITION
        assert selectedTargers.ids[1] == NULL_TARGET_ID
        assert selectedTargers.priorities[1] == NULL_TARGET_PRIORITY
