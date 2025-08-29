"""

Some utility methods related with PFS science targets.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

"""

import numpy as np

from .TargetGroup import TargetGroup


def generateOneTargetPerCobra(bench, maximumDistance=np.inf):
    """Generates a single target per cobra using a uniform radial distribution.

    Parameters
    ----------
    bench: object
        The PFI bench instance.
    maximumDistance: float, optional
        The maximum radial distance between the targets and the cobra centers.
        Default is no limit.

    Returns
    -------
    Object
        A TargetGroup instance.

    """
    # Extract some useful information
    nCobras = bench.cobras.nCobras
    cobraCenters = bench.cobras.centers
    rMin = bench.cobras.rMin
    rMax = bench.cobras.rMax

    # Calculate the maximum target distance allowed for each cobra
    rMax = rMax.copy()
    rMax[rMax > maximumDistance] = maximumDistance
    rMax[rMax < rMin] = rMin[rMax < rMin]

    # Calculate the target positions
    ang = 2 * np.pi * np.random.random(nCobras)
    radius = np.sqrt(
        (rMax ** 2 - rMin ** 2) * np.random.random(nCobras) + rMin ** 2)
    finalPositions = cobraCenters + radius * np.exp(1j * ang)

    # Set to NULL those targets where the maximum distance is smaller than rMin
    finalPositions[maximumDistance < rMin] = TargetGroup.NULL_TARGET_POSITION

    return  TargetGroup(finalPositions)


def generateRandomTargets(density, bench):
    """Generates a set of targets uniformly distributed over the bench field of
    view.

    Parameters
    ----------
    density: float
        The average number of targets per cobra patrol area.
    bench: object
        The PFI bench instance.

    Returns
    -------
    Object
        A target group instance.

    """
    # Calculate the total number of targets based on the bench properties
    averagePatrolArea = np.mean(np.pi * bench.cobras.rMax ** 2)
    benchArea = np.pi * bench.radius ** 2
    nTargets = int(np.ceil(density * benchArea / averagePatrolArea))

    # Calculate the uniformly distributed target positions
    ang = 2 * np.pi * np.random.random(nTargets)
    radius = bench.radius * np.sqrt(np.random.random(nTargets))
    finalPositions = bench.center + radius * np.exp(1j * ang)

    return TargetGroup(finalPositions)
