"""

Collection of fixtures for the unit tests.

"""

import pytest
import numpy as np

from ics.cobraOps.Bench import Bench
from ics.cobraOps import targetUtils


@pytest.fixture(scope="function")
def targetPositions():
    return np.array([-1j, 1 + 1j, 1 - 1j, 2j])


@pytest.fixture(scope="function")
def targetIds():
    return np.array(["source00+11", "source_b", "+c", "-d"])


@pytest.fixture(scope="function")
def targetPriorities():
    return np.array([1, 1, 2, 8])


@pytest.fixture(scope="function")
def bench():
    return Bench(layout="full")


@pytest.fixture(scope="function")
def targets(bench):
    return targetUtils.generateRandomTargets(2, bench)
