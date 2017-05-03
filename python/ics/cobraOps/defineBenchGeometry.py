
import numpy as np

from loadCfgXml import loadCfgXml


def defineBenchGeometry(centers, useRealMaps, useRealLinks):
    """Defines the bench geometry.
    
    thetaDIR is > 0 for positive move out of home (opposite of phi)    
    thetaDIR is < 0 for negative move out of home (same as phi)
    
    output is a structure with
        center (Mx1)
        L1     (Mx1)
        L2     (Mx1)
        phiIn  (Mx1)
        phiOut (Mx1)
        tht0   (Mx1)
        NNmap  (MxM), logical
        rf
        distCobras
        minDist
        pids 
        mids 
        S1Nm 
        S1Pm, S2Pm, S2Nm, binWidth);
        
    """
    KeepOutAngle = 0.1  # radians (~ 5.5 deg)
    alpha = 0.07
    beta = 0.50
    pids = np.array([], dtype="int")
    mids = np.array([], dtype="int")
    mapAssignment = None

    # This is an epsilon to make sure that values near zero in theta are
    # intepreted on the positive (or negative) side of the cut,
    # respectively. Make it negative for same same direction moveouts (positive
    # for positive hardstop).
    #
    # tht1 implemented below is an opposite-sense hard stop.  when
    # that's working, theta direction will be specified on a per cobra
    # basis, not on a full instrument basis, as is implemented here.
    thteps = 2e-10  # vesitigial

    if useRealMaps or useRealLinks:
        # Load the configuration file
        cobraConfig = loadCfgXml("updatedMotorMapsFromThisRun2.xml")
        
        # Get the arm data containers
        dataContainers = cobraConfig.findall("ARM_DATA_CONTAINER")
        nDataContainers = len(dataContainers)

        # Get the number of steps from the first data container        
        slowCalibrationTable = dataContainers[0].find("SLOW_CALIBRATION_TABLE")
        
        if slowCalibrationTable.find("Joint1_fwd_stepsizes") is not None:
            nSteps = int(slowCalibrationTable.find("Joint1_fwd_stepsizes").text.split(",")[0])
        else:
            nSteps = 0

        # Save some of the configuration information
        pids = np.zeros(nDataContainers)
        mids = np.zeros(nDataContainers)
        cx = np.zeros(nDataContainers)
        cy = np.zeros(nDataContainers)
        tht0 = np.zeros(nDataContainers)  # hard stop angle for same sense move-out
        tht1 = np.zeros(nDataContainers)  # hard stop angle for opposite sense move-out
        phiIn = np.zeros(nDataContainers)
        phiOut = np.zeros(nDataContainers)
        L1 = np.zeros(nDataContainers)
        L2 = np.zeros(nDataContainers)
        S1Pm = np.zeros((nDataContainers, nSteps))
        S2Pm = np.zeros((nDataContainers, nSteps))
        S1Nm = np.zeros((nDataContainers, nSteps))
        S2Nm = np.zeros((nDataContainers, nSteps))

        for i in xrange(nDataContainers):
            dataContainer = dataContainers[i]
            
            header = dataContainer.find("DATA_HEADER")
            pids[i] = int(header.find("Positioner_Id").text)
            mids[i] = int(header.find("Module_Id").text)
            
            kinematics = dataContainer.find("KINEMATICS")
            cx[i] = float(kinematics.find("Global_base_pos_x").text)
            cy[i] = float(kinematics.find("Global_base_pos_y").text)
            tht0[i] = float(kinematics.find("CCW_Global_base_ori_z").text) * np.pi / 180
            tht1[i] = float(kinematics.find("CW_Global_base_ori_z").text) * np.pi / 180
            phiIn[i] = float(kinematics.find("Joint2_CCW_limit_angle").text) * np.pi / 180 - np.pi
            phiOut[i] = float(kinematics.find("Joint2_CW_limit_angle").text) * np.pi / 180 - np.pi
                        
            if kinematics.find("Link1_Link_Length") is not None:
                L1[i] = float(kinematics.find("Link1_Link_Length").text)
            else:
                L1[i] = 25.0

            if kinematics.find("Link2_Link_Length") is not None:
                L2[i] = float(kinematics.find("Link2_Link_Length").text)
            else:
                L2[i] = 25.0
            
            slowCalibrationTable = dataContainer.find("SLOW_CALIBRATION_TABLE")

            if slowCalibrationTable.find("Joint1_fwd_stepsizes") is not None:
                # Regions are expected to be standard bin width of 3.6deg
                # phi out of HS in positive direction
                # tht out of HS in negative direction
                Joint1Fwd = slowCalibrationTable.find("Joint1_fwd_stepsizes").text.split(",")[2:-1]
                Joint2Fwd = slowCalibrationTable.find("Joint2_fwd_stepsizes").text.split(",")[2:-1]
                Joint1Rev = slowCalibrationTable.find("Joint1_rev_stepsizes").text.split(",")[2:-1]
                Joint2Rev = slowCalibrationTable.find("Joint2_rev_stepsizes").text.split(",")[2:-1]
                S1Pm[i] = [3.6 / float(val) for val in Joint1Fwd]
                S2Pm[i] = [3.6 / float(val) for val in Joint2Fwd]
                S1Nm[i] = [3.6 / float(val) for val in Joint1Rev]
                S2Nm[i] = [3.6 / float(val) for val in Joint2Rev]

        # Calculated map ranges (not physical range of motion).  other
        # option is to read it directly from regions.
        map_range = {}
        map_range["tht"] = np.array([0, nSteps * 3.6 * np.pi / 180])
        map_range["phi"] = np.array([0, nSteps * 3.6 * np.pi / 180]) - np.pi

        phiOut -= KeepOutAngle
        phiIn = np.maximum(phiIn, -np.pi)  # HACK SOLUTION TO BAD PHI HOMES [2016-07-15]
        phiIn += KeepOutAngle
        center = cx + cy * 1j

        # Physical parameters
        rf = 1000.0 / 90  # fiber holder radius
        distCobras = 8.0 * 1000.0 / 90  # distance between cobras for custom geometries.
        minDist = 2 * rf  # minimum distance between fiber and anything
        
        # variables defined in this branch that can propagate to different centers:
        # L1
        # L2
        # phiIn
        # phiOut
        # S[12][FR]m
        configData = {}
        configData["pids"] = pids
        configData["mids"] = mids
        configData["phiIn"] = phiIn
        configData["phiOut"] = phiOut
        configData["L1"] = L1
        configData["L2"] = L2
        configData["S1Pm"] = S1Pm
        configData["S2Pm"] = S2Pm
        configData["S1Nm"] = S1Nm
        configData["S2Nm"] = S2Nm
        configData["distCobras"] = distCobras

    if centers is not None:
        # Cobra center locations
        center = centers
        csize = len(center)
        ONE = np.ones(csize)
    
        # Physical parameters
        rf = 1.0  # fiber holder radius
        distCobras = 8.0  # distance between cobras for custom geometries.
        minDist = 2 * rf  # minimum distance between fiber and anything

        thtrange = 385 * np.pi / 180  # range of motion of the theta stage.

        tht0 = 2 * np.pi * np.random.random(csize)  # same sense hard stop
        tht1 = (tht0 + thtrange) % (2 * np.pi)  # opp sense hard stop

        if useRealMaps:
            nDataContainers = configData["S1Pm"].shape[0]
            nSteps = configData["S1Pm"].shape[1]
            S1Pm = np.zeros((csize, nSteps))
            S2Pm = np.zeros((csize, nSteps))
            S1Nm = np.zeros((csize, nSteps))
            S2Nm = np.zeros((csize, nSteps))
            
            mapAssignment = (nDataContainers * np.random.random((4, csize))).astype("int")
            
            S1Pm = configData["S1Pm"][mapAssignment[0]]
            S2Pm = configData["S2Pm"][mapAssignment[1]]
            S1Nm = configData["S1Nm"][mapAssignment[2]]
            S2Nm = configData["S2Nm"][mapAssignment[3]]
        else:
            S1Pm = np.array([])
            S2Pm = np.array([])
            S1Nm = np.array([])
            S2Nm = np.array([])
        
        if useRealLinks:
            # take L1,L2, phiIn/Out from CobraConfig.
            pix2mm = distCobras / configData["distCobras"]
            nDataContainers = configData["L1"].shape[0]

            indexes = (nDataContainers * np.random.random((4, csize))).astype("int")
            
            phiIn = configData["phiIn"][indexes[0]] 
            phiOut = configData["phiOut"][indexes[1]] 
            L1 = configData["L1"][indexes[2]] * pix2mm              
            # ## There was a bug here it was using L1
            L2 = configData["L2"][indexes[3]] * pix2mm              
        else: 
            DeltaPhi = np.pi - 2 * KeepOutAngle  # throw of the phi arm
            phiOut = -KeepOutAngle * (1.0 + 0.2 * np.random.random(csize))
            phiIn = phiOut - DeltaPhi

            linkLength = 2.375
            L1 = linkLength * np.ones(csize)
            L2 = linkLength * np.ones(csize)

    rMin = np.abs(L1 + L2 * np.exp(1j * phiIn))
    rMax = np.abs(L1 + L2 * np.exp(1j * phiOut))
    
    binWidth = 2 * np.pi / 100;  # this should be defined higher up in the code.[PHM 20160901
    
    # parameters for populating the annular patrol areas uniformly
    dA = 1.0 / ((rMax / rMin) ** 2 - 1.0)  # fractional keepout area
    rRange = np.sqrt(rMax ** 2 - rMin ** 2)

    # by default, phi home is phi in.  existence of mat files *may* change phi home.
    phiHome = phiIn.copy()

    # for the realities of the test bench, having phi too far in is
    # inconvenient.  Use a larger value so that the theta arm does not
    # go crazy.
    phiHome += 0.5

    # read in phiHome from elsewhere
    """
    # TODO!!!
    matfiles = ls2cell('*.mat')
    
    if len(matfiles) > 0:
        print('Warning -- using mat file positions for home in defineBenchGeometry')
        data = loadTestData('', CobraConfig, 1) 
        
        for id = 1:numel(data):
            if ~isempty(data(id).xst):
                kk = find(pids == id)
                phiHome(kk) = mean(data(id).s2.startAngle(1,:)) * pi / 180 - pi
    """
    
    # home positions (0 = ss, 1 = os)
    home0 = center + L1 * np.exp(1j * tht0) + L2 * np.exp(1j * (tht0 + phiHome))
    home1 = center + L1 * np.exp(1j * tht1) + L2 * np.exp(1j * (tht1 + phiHome))

    # theta overlap
    tht_overlap = ((tht1 - tht0 + np.pi) % (2 * np.pi)) - np.pi;

    # generate the nearest neighbor map/structure
    csize = len(center)
    epsilon = 1e-9  # small number to take care of roundoff errors.
    CCdist = np.zeros((csize, csize))
    
    for i in xrange(csize):
        CCdist[i] = np.abs(center - center[i])
        
    distCent = 2 * np.median(L1 + L2) * (8.0 / 9.5)
    nnMap = np.logical_and(CCdist < 1.5 * distCent , CCdist > epsilon)
    (row, col) = np.where(nnMap)
    # inUpperTriangle = row < col
    NN = {}
    NN["row"] = row  # (inUpperTriangle)
    NN["col"] = col  # (inUpperTriangle)
    NN["xy"] = (center[row] + center[col]) / 2

    # Calculate some parameters for the field geometry
    field = {}
    field["cm"] = np.mean(center)
    field["R"] = np.max(np.abs(center - field["cm"]) + rMax)
    # Fit a line to the centers and calculate the slope
    x = np.real(center)
    y = np.imag(center)
    slope = (csize * (x * y).sum() - x.sum() * y.sum()) / (csize * (x * x).sum() - x.sum() ** 2);
    field["ang"] = np.arctan([slope])

    output = {}
    output["center"] = center
    output["L1"] = L1
    output["L2"] = L2
    output["phiIn"] = phiIn
    output["phiOut"] = phiOut
    output["tht0"] = tht0
    output["tht1"] = tht1
    output["rMin"] = rMin
    output["rMax"] = rMax
    output["dA"] = dA
    output["rRange"] = rRange
    output["home0"] = home0
    output["home1"] = home1
    output["nnMap"] = nnMap
    output["NN"] = NN
    output["rf"] = rf
    output["distCobras"] = distCobras
    output["minDist"] = minDist
    output["S1Nm"] = S1Nm
    output["S1Pm"] = S1Pm
    output["S2Pm"] = S2Pm
    output["S2Nm"] = S2Nm
    output["map_range"] = map_range
    output["binWidth"] = binWidth
    output["alpha"] = alpha
    output["beta"] = beta
    output["pids"] = pids
    output["mids"] = mids
    output["thteps"] = thteps
    output["tht_overlap"] = tht_overlap
    output["field"] = field
    output["map"] = mapAssignment

    return output
