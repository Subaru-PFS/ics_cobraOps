
import numpy as np


def generateTargets(density, bench, makefigs=True):
    """Generate a uniform density set of targets over the field of view of the bench.

    Parameters
    ----------
    density: int
        The # targets per patrol area.
    bench: object
        The bench to use. If None, a full bench will be generated.
        Default is None.
    makefigs: bool
        If true, some figures will be made. Default is True.
    
    Returns
    -------
    Object
        The target list.
            
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

    # uniformly generate targets over field of view.
    ntgt = int(np.ceil(density * (bench["field"]["R"] / np.median(bench["rMax"])) ** 2))
    print("Generating " + str(ntgt) + " targets for " + str(ncobras) + " positioners")

    THT = np.random.random(ntgt) * 2 * np.pi
    RR = np.sqrt(np.random.random(ntgt)) * bench["field"]["R"]
    tgt = RR * np.exp(1j * THT)  # generated targets, centered on cobra CM
    tgt += bench["field"]["cm"]

    return tgt
