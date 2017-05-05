
import numpy as np
import matplotlib.pyplot as plt

import cobraUtils as cobraUtils
import benchUtils as benchUtils
from assign_targets import assign_targets

def simFun(numtrg=1, cobraLayout="none", useRealMaps=True, useRealLinks=True, varargin=None):
    """Cobras collision simulator for PFS.

    Parameters
    ----------
    numtrg: Object
        Number of targets to generate. It could be a number describing the target density 
        to simulate or a numpy array with the target positions. Default is 1.
    cobraLayout: str
        The cobras layout to use: "none", "hex", "line", "rails" or "full". "none" means
        that the layout should be extracted from a bench configuration file. Default is "none".
    useRealMaps: bool
        If true, use xml data for maps.
    useRealLinks: bool
        If true, use xml data for geometry.
    varargin: dict
        List of optional configuration parameters.
    
    Returns
    -------
    Object
        Complex object with the following elements:
        - targets: Nx1 complex array of target positions.
        - Traj: trajectory structure from realizeTrajectory2.
        - Coll: collision structure from detectCollisionsSparse.
        - bench: bench definition from defineBenchGeometry.
        - minDist: 6Nx1 array of minimum distances between nearest neighbors.
        - caats: "collide at any time" list of nearest neighbors.
        - IR1_colliders: not sure.
        - IR1: something to do with interference replan 1.
    
    """
    # Define the application control parameters defaults
    toggle = {}
    toggle["info"] = "Logical switches"
    toggle["showMoves"] = False
    toggle["showFigures"] = True
    toggle["SkipTargetReplan"] = True
    toggle["verbosity"] = 0

    # Check if the user provided some special configuration
    if varargin is not None:
        for key in toggle:
            if key in varargin:
                toggle[key] = varargin[key]
    
    # Define the bench
    if varargin is not None and "UseThisBench" in varargin:
        bench = varargin["UseThisBench"]
    else:
        # Get the cobras central positions for the requested layout
        centers = cobraUtils.getCobrasCenters(cobraLayout)

        # Plot the centers if necessary
        if toggle["showFigures"]:
            cobraUtils.plotCobrasCenters(centers)

        # Define the bench geometry
        bench = benchUtils.defineBenchGeometry(centers, useRealMaps, useRealLinks)

    # Reassign the bench alpha value if desired
    if varargin is not None and "alpha" in varargin:
        bench["alpha"] = varargin["alpha"] 

    alpha = bench["alpha"]

    numPos = len(bench["center"])

    # Initialize the performance metrics
    PM = {}
    PM["info"] = "Performance metrics"
    PM["total_targets"] = 0
    PM["total_collisions"] = 0
    PM["total_primary_collisions"] = 0

    if isinstance(numtrg, int) or isinstance(numtrg, float):
        # uniform over... ('field|'patrol')
        TGT_GEN_STRATEGY = "field" 
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
            # targets[i] = assignments["tgt"]
        
        PM["R2_percentColl"] = None
    elif TGT_GEN_STRATEGY == "patrol":
        # TODO See matlab code
        raise Exception("Impossible path: TGT_GEN_STRATEGY = 'patrol'")

    """
    if toggle["showFigures"] and centers is not None:
        plt.figure("Tarject centers", facecolor="white")
        plt.scatter(np.real(targets), np.imag(targets), s=2)
        plt.xlabel("x position")
        plt.ylabel("y position")
        plt.title("Tarject centers (" + TGT_GEN_STRATEGY + ")")
        plt.show(block=False)
    """
    
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

    plt.show()

    # TBD...
    return bench
    
#aa = simFun(cobraLayout="rails")
aa = benchUtils.getBenchCalibrationData("updatedMotorMapsFromThisRun2.xml")
print(len(aa["centers"]))
print(len(aa["L1"]))
print(aa)
