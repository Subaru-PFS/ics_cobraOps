
import numpy as np
import matplotlib.pyplot as plt

from getCentersRails import getCentersRails
from defineBenchGeometry import defineBenchGeometry
from assign_targets import assign_targets


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
        List of optional parameters.
    
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
        for key in toggle:
            if key in varargin:
                toggle[key] = varargin[key]

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
            centers = None
            useRealMaps = True
        elif cobraLayout == "hex":
            # The centers will follow a hexagon pattern
            # (R3 collisions = 3.0e-4 w/ 100k targets)
            centers = np.zeros(7, dtype="complex")
            centers[1:] = 8 * np.exp(np.arange(6) * 1j * np.pi / 3)
        elif cobraLayout == "line":
            # Line of cobras 
            # (R3 collisions = 4.5e-5 (fraction) w/ 1M targets)
            centers = 8 * np.arange(27) + 1j
        elif cobraLayout == "rails":
            # N hard coded rails 
            # (R3 collisions = 9.9e-5 w/ 16k targets)
            centers = getCentersRails(14)
        elif cobraLayout == "full":
            # Full PFI bench 
            # (R3 collisions = 1.7e-4 (fraction) w/ 1M targets)
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
        else:
            print("Valid cobra layouts: 'none', 'hex', 'line', 'rails' or 'full'")
            print("Exiting...")
            return

        # Plot the centers
        if toggle["showFigures"] and centers is not None:
            plt.scatter(np.real(centers), np.imag(centers), s=2)
            plt.xlabel("x position")
            plt.ylabel("y position")
            plt.title("Cobra centers (" + cobraLayout + ")")
            plt.show()

        # Define a structure with the bench geometry
        bench = defineBenchGeometry(centers, useRealMaps, useRealLinks)
    
    #----------------------------
    #    bench defined      
    #----------------------------

    # reassign bench alpha/beta here (if desired)
    if varargin is not None and "alpha" in varargin:
        alpha = varargin["alpha"]
        bench["alpha"] = alpha 
    else:
        alpha = bench["alpha"]

    numPos = len(bench["center"])

    if isinstance(numtrg, int):
        TGT_GEN_STRATEGY = "field"  # uniform over... ('field|'patrol')
    else:
        TGT_GEN_STRATEGY = "targetlist"

    assignments = {}
    
    if TGT_GEN_STRATEGY == "targetlist":
        numFields = 1
        targets = numtrg.copy()
        assignments["tgt"] = targets  # for partial compatibility with 'field' case.
    elif TGT_GEN_STRATEGY == "field":
        # Use numtrg as density
        numFields = 1
        targets = np.zeros((numFields, numPos), dtype="complex")
        
        for i in xrange(0, numFields):
            assignments = assign_targets(numtrg, bench)
            #targets[i] = assignments["tgt"]
        
        PM["R2_percentColl"] = None
    elif TGT_GEN_STRATEGY == "patrol":
        # TODO See matlab code
        print("Impossible path!!!")
    else:
        print("Target strategy is uniform over 'field' or uniform over 'patrol'")
        return

    '''
    if toggle["showFigures"] and centers is not None:
        plt.scatter(np.real(targets), np.imag(targets), s=2)
        plt.xlabel("x position")
        plt.ylabel("y position")
        plt.title("Tarject centers (" + TGT_GEN_STRATEGY + ")")
        plt.show()
    '''
    
    #------------------------------------------------------------------
    #    Targets defined, no end-point physical interferences (Rule 2)      
    #------------------------------------------------------------------

    # Decide on which hard stop to use.
    # There is a ~25 deg overlap in the tht patrol range.  to prevent having to go after targets just
    # on the wrong side of the HS, set the allowable range from DeltaTheta/2 to MaxThtTravel (~Pi).
    #
    # $$$ targetsTP = XY2TP(bsxfun(@minus, targets, bench.center), bench.L1, bench.L2);
    # delta angles out of same (0) and opposite (1) set ups.
    # $$$ dtht0 = mod(targetsTP.tht - bench.tht0, 2*pi)
    # $$$ dtht1 = mod(targetsTP.tht - bench.tht1, 2*pi) - 2*pi
    # $$$ tht_overlap = mod(bench.tht1 - bench.tht0,2*pi)
    # $$$ max_tht_trvl = pi + tht_overlap
    # in deciding which hard stop to use, we allow some overlap both at small angles (1/3 of small
    # overlap angle) and large angles (pi + overlap)
    # what follows are logical array cuts on the positioners.
    # $$$ ss_allowed =  dtht0 > tht_overlap/3 &  dtht0 < max_tht_trvl
    # $$$ os_allowed = -dtht1 > tht_overlap/3 & -dtht1 < max_tht_trvl 
    # $$$ only_ss = ss_allowed & ~os_allowed
    # $$$ only_os = os_allowed & ~ss_allowed
    # $$$ s_or_os = ss_allowed &  os_allowed
    # these take the shorter path, exept when it's too close to zero.
    # $$$ use_ss  = s_or_os & (dtht0 < -dtht1) | only_ss
    # $$$ use_os  = s_or_os & (dtht0 > -dtht1) | only_os

    # TBD...
    return bench
    
aa = simFun(cobraLayout="full")
print(len(aa["center"]))
print(len(aa["L1"]))
