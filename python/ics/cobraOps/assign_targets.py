
import numpy as np

from getCentersRails import getCentersRails
from defineBenchGeometry import defineBenchGeometry
from generateTargets import generateTargets


def assign_targets(tgt, bench=None, makefigs=False):
    """Assign targets to cobras.

    Parameters
    ----------
    tgt: object
        Target list (complex) or density (real scalar).
    bench: object
        The bench to use. If None, a full bench will be generated.
        Default is None.
    makefigs: bool
        If true, some figures will be made. Default is False.
    
    Returns
    -------
    Object
        ...
            
    """
    # Define the bench if necessary
    if bench is None:
        # Get the first sector
        firstSector = getCentersRails(14)
            
        # Create the centers array
        cobrasPerSector = len(firstSector)
        centers = np.zeros(3 * cobrasPerSector, dtype="complex")
            
        # Add the first sector
        centers[:cobrasPerSector] = firstSector 
            
        # Add the second sector rotating the first sector 120 degrees
        centers[cobrasPerSector:2 * cobrasPerSector] = firstSector * np.exp(1j * 2 * np.pi / 3) 

        # Add the third sector rotating the first sector 240 degrees
        centers[2 * cobrasPerSector:] = firstSector * np.exp(1j * 4 * np.pi / 3)

        bench = defineBenchGeometry(centers, True, True)
        bench["alpha"] = 0

    ncobras = len(bench["center"])

    if (isinstance(tgt, int) or isinstance(tgt, float)):
        tgtdensity = tgt
        tgt = generateTargets(tgtdensity, bench)
   
    ntgt = len(tgt)

    # Check min/max radius compliance
    LL = np.zeros((ncobras, ntgt))
    
    for i in xrange(ncobras):
        LL[i] = np.abs(tgt - bench["center"][i]) 
    
    LL = LL < np.max(bench["rMax"])  # Who is in patrol area

    (cc, tt) = np.where(LL)
   
    XY = np.zeros(LL.shape, dtype="complex")  # initialize XY with zeros

    for i in xrange(len(cc)):  # for every cobra check assigment 
        xy = tgt[tt[i]] - bench["center"][cc[i]]  # get the local coordinate
        dst = np.abs(xy)
        
        if (dst > bench["rMin"][cc[i]] and dst < bench["rMax"][cc[i]]):
            XY[cc[i], tt[i]] = xy
    
    output = {}
    output["tgt"] = tgt

    return output
