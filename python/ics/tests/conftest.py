"""

Collection of fixtures for the unit tests.

"""

import os
import pytest
import numpy as np

from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach

from ics.cobraOps.Bench import Bench
from ics.cobraOps import targetUtils


@pytest.fixture(scope="function")
def finalPositions():
    return np.array([-1j, 1 + 1j, 1 - 1j, 2j])


@pytest.fixture(scope="function")
def targetIds():
    return np.array(["source00+11", "source_b", "+c", "-d"])


@pytest.fixture(scope="function")
def targetPriorities():
    return np.array([1, 1, 2, 8])


@pytest.fixture(scope="function")
def bench():
    os.environ["PFS_INSTDATA_DIR"] = "/home/jgracia/github/pfs_instdata"
    cobraCoach = CobraCoach(
        loadModel=True, trajectoryMode=True, rootDir="/home/jgracia/testPFI/")
    calibrationFileName = os.path.join(
       os.environ["PFS_INSTDATA_DIR"],"data/pfi/dot", "black_dots_mm.csv")
    blackDotsCalibrationProduct = BlackDotsCalibrationProduct(calibrationFileName)

    return Bench(cobraCoach, blackDotsCalibrationProduct)


@pytest.fixture(scope="function")
def targets(bench):
    return targetUtils.generateRandomTargets(2, bench)
