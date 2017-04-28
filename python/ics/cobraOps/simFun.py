
import cmath

def getCentersRails(numOfRails):
    # Create the first line
    centers = [complex(8 * i, 0) for i in xrange(0, 29)]
    
    # Add a second line
    centers.extend([i + 8 * cmath.exp(1j * cmath.pi / 3) for i in centers[:-1]])
    
    # Shift the two lines
    centers = [i + 8 * cmath.exp(1j * 2 * cmath.pi / 3) for i in centers]
    
    # Order them by the x coordinate
    centers = sorted(centers, key=lambda x: x.real)
    
    # Create all the rails iteratively
    rails = centers[:]
    
    for i in xrange(1, numOfRails):
        rails.extend([c + i * 16 * cmath.exp(1j * 2 * cmath.pi / 3) for c in centers]);
    
    return rails


def defineBenchGeometry():
    return 


def simFun(numtrg=1, cobraLayout="none", useRealMaps=True, useRealLinks=True, varargin=None):
    """Collision simulator for PFS.

    Parameters
    ----------
    numtrg: int
        Number of targets to generate. Default is 1.
    cobraLayout: str
        The cobras layout to use: "none", "hex", "line", "rails" or "full". "none" means
        that the layout should be extracted from a bench configuration file. Default is "none".
    useRealMaps: bool
        If true, use xml data for maps.
    useRealLinks: bool
        If true, use xml data for geometry.
    varargin: dict
        list of optional parameters
    
    Returns
    -------
    Object
        Complex object with the following elements:
            targets: Nx1 complex array of target positions
            Traj: trajectory structure from realizeTrajectory2
            Coll: collision structure from detectCollisionsSparse
            bench: bench definition from defineBenchGeometry
            minDist: 6Nx1 array of minimum distances between nearest neighbors
            caats: "collide at any time" list of nearest neighbors
            IR1_colliders: not sure
            IR1: something to do with interference replan 1
    
    """
    # General application control dictionary
    toggle = {}
    toggle["info"] = "Logical switches"
    toggle["showMoves"] = False
    toggle["showFigures"] = True
    toggle["SkipTargetReplan"] = True
    toggle["verbosity"] = 0

    # Check if there is some special configuration
    if varargin is not None:
        if "alpha" in varargin:
            alpha = varargin["alpha"]
        if "showMoves" in varargin:
            toggle["showMoves"] = varargin["showMoves"]
        if "SkipTargetReplan" in varargin:
            toggle["SkipTargetReplan"] = varargin["SkipTargetReplan"]
        if "verbosity" in varargin:
            toggle["verbosity"] = varargin["verbosity"]

    # Initialize the performance metrics
    PM = {}
    PM["info"] = "Performance metrics"
    PM["total_primary_collisions"] = 0
    PM["total_collisions"] = 0
    PM["total_targets"] = 0
    
    # Define the bench
    if varargin is not None and "UseThisBench" in varargin:
        bench = varargin["UseThisBench"]
    else:
        if cobraLayout == "none":
            # Use the configuration file geometry
            centers = []
            useRealMaps = True
        elif cobraLayout == "hex":
            # The centers will follow a hexagon pattern
            # (R3 collisions = 3.0e-4 w/ 100k targets)
            centers = [complex(0, 0)]
            centers.extend([8 * cmath.exp(i * 1j * cmath.pi / 3) for i in xrange(0, 6)])            
        elif cobraLayout == "line":
            # Line of cobras 
            # (R3 collisions = 4.5e-5 (fraction) w/ 1M targets)
            centers = [complex(8 * i , 1) for i in xrange(0, 27)]
        elif cobraLayout == "rails":
            # N hard coded rails 
            # (R3 collisions = 9.9e-5 w/ 16k targets)
            centers = getCentersRails(14)
        elif cobraLayout == "full":
            # Full PFI bench 
            # (R3 collisions = 1.7e-4 (fraction) w/ 1M targets)
            tempCenters = getCentersRails(14)
            centers = tempCenters[:]
            centers.extend([i * cmath.exp(1j * 2 * cmath.pi / 3) for i in tempCenters])
            centers.extend([i * cmath.exp(1j * 4 * cmath.pi / 3) for i in tempCenters])
        else:
            print("Valid cobra layouts: 'none', 'hex', 'line', 'rails' or 'full'")
            print("Exiting...")
            return

        # Define a structure with the bench geometry
        bench = defineBenchGeometry(centers, useRealMaps, useRealLinks)

    # TBD...

    