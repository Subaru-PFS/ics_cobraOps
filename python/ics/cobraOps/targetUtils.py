"""

Some utility methods related with the science targets.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np

from ics.cobraOps.TargetGroup import TargetGroup


def generateOneTargetPerCobra(bench, maximumDistance=np.Inf):
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
        A target group instance.
    
    """
    # Extract some useful information
    nCobras = bench.cobras.nCobras
    cobraCenters = bench.cobra.centers
    rMin = bench.cobras.rMin
    rMax = bench.cobras.rMax
    
    # Calculate the maximum target distance allowed for each cobra
    rMax = rMax.copy()
    rMax[rMax > maximumDistance] = maximumDistance
    rMax[rMax < rMin] = rMin[rMax < rMin]
    
    # Calculate the target positions
    ang = 2 * np.pi * np.random.random(nCobras)
    radius = np.sqrt((rMax ** 2 - rMin ** 2) * np.random.random(nCobras) + rMin ** 2)
    targetPositions = cobraCenters + radius * np.exp(1j * ang)
    
    return  TargetGroup(targetPositions)


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
    medianPatrolAreaRadius = np.median(bench.cobras.rMax)
    nTargets = int(np.ceil(density * (bench.radius / medianPatrolAreaRadius) ** 2))
    
    # Calculate the uniformly distributed target positions
    ang = 2 * np.pi * np.random.random(nTargets)
    radius = bench.radius * np.sqrt(np.random.random(nTargets))
    targetPositions = bench.center + radius * np.exp(1j * ang)
    
    return TargetGroup(targetPositions)
