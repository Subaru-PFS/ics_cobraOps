"""

Collection of fixtures for the unit tests.

"""

import pytest
import numpy as np

from ics.cobraOps.cobraConstants import NULL_TARGET_ID


@pytest.fixture(scope="function")
def targetPositions():
    return np.array([-1j, 1 + 1j, 1 - 1j, 2j])


@pytest.fixture(scope="function")
def targetIds():
    return np.array(["a", NULL_TARGET_ID, "c", "d"])


@pytest.fixture(scope="function")
def targetPriorities():
    return np.array([1, 1, 2, 3])
