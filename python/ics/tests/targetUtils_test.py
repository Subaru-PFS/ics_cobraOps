"""

Collection of unit tests for the targetUtils module.

"""

import numpy as np

from ics.cobraOps import targetUtils
from ics.cobraOps.cobraConstants import NULL_TARGET_ID


class TestTargetUtils():
    """A collection of tests for the targetUtils module.

    """

    def test_generateOneTargetPerCobra_method(self, bench):
        # Generate the targets
        targets = targetUtils.generateOneTargetPerCobra(bench)

        # Check that there are as many targets as cobras in the bench
        assert targets.nTargets == bench.cobras.nCobras

        # Check that the target separation is correct
        distances = np.abs(targets.positions - bench.cobras.centers)
        rMin = bench.cobras.rMin
        rMax = bench.cobras.rMax
        epsilon = 1e-6
        assert np.all(np.logical_and(distances + epsilon > rMin,
                                     distances - epsilon < rMax))

        # Generate a new set of targets with some maximum distance limit
        maximumDistance = 5
        targets = targetUtils.generateOneTargetPerCobra(bench, maximumDistance)

        # Check that the target separation is correct
        distances = np.abs(targets.positions - bench.cobras.centers)
        rMax = rMax.copy()
        rMax[rMax > maximumDistance] = maximumDistance
        rMax[rMax < rMin] = rMin[rMax < rMin]
        assert np.all(np.logical_and(distances + epsilon > rMin,
                                     distances - epsilon < rMax))

        # Generate a new set of targets with zero maximum distance
        maximumDistance = 0
        targets = targetUtils.generateOneTargetPerCobra(bench, maximumDistance)

        # Check that all the targets are NULL targets
        assert np.all(targets.ids == NULL_TARGET_ID)

    def test_generateRandomTargets_method(self, bench):
        # Generate the targets
        density = 5
        targets = targetUtils.generateRandomTargets(density, bench)

        # Check that the targets fall within the bench
        distances = np.abs(targets.positions - bench.center)
        assert np.all(distances <= bench.radius)

        # Check that we get the expected density at different test radii
        testRadii = bench.radius * np.linspace(0.2, 1, 9)
        testTargets = np.sum(distances < testRadii[:, np.newaxis], axis=1)
        testDensities = testTargets * (
            np.mean(bench.cobras.rMax ** 2) / testRadii ** 2)
        assert np.all(np.abs(testDensities - density) < 0.2 * density)
