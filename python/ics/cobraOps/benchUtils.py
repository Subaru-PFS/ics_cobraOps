"""

Some utility methods related with the PFS bench.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T
  
"""

import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET


KEEP_OUT_ANGLE = 0.1
"""Safety angle value in radians (~ 5.7 deg)."""

ALPHA = 0.07
"""The motor noise normalization factor."""

BETA = 0.50
"""The motor noise exponent."""

LINK_LENGTH = 2.375
"""The theoretical link length in mm."""
 
COBRAS_SEPARATION = 8.0
"""The theoretical separation in mm between two consecutive cobras."""
           
THT_EPS = 2e-10
"""This is an epsilon to make sure that values near zero in theta are
intepreted on the positive (or negative) side of the cut, respectively. 
Make it negative for same same direction moveouts (positive for 
positive hardstop).

tht1 implemented below is an opposite-sense hard stop.  when that's 
working, theta direction will be specified on a per cobra basis, not 
on a full instrument basis, as is implemented here."""


def getBenchCalibrationData(fileName):
    """Returns some bench calibration data stored in a XML file.

    Parameters
    ----------
    fileName: str
        The path to the XML calibration file.
    
    Returns
    -------
    Object
        The bench calibration data, containing the following elements:
        - mids: The module ids (n).
        - pids: The positioner ids (n).
        - centers: The central positions in pixel units (n).
        - tht0: The hard stop angle for same sense move-out in radians (n).
        - tht1: The hard stop angle for opposite sense move-out in radians (n).
        - phiIn: ... in radians (n).
        - phiOut: ... in radians (n).
        - pixelScale: The pixel size in millimeters (n).
        - L1: The link1 lengths in pixel units (n).
        - L2: The link2 lengths in pixel units (n).
        - angularStep: The angular steps size in degrees used in the motor measurements (n).
        - mapRangeTht: The theta angle map range in radias (nx2).
        - mapRangePhi: The phi angle map range in radias (nx2).
        - S1Pm: The joint1 forward motor steps required to move an angular step (nxm).
        - S1Nm: The joint1 reverse motor steps required to move an angular step (nxm).
        - S2Pm: The joint2 forward motor steps required to move an angular step (nxm).
        - S2Nm: The joint2 reverse motor steps required to move an angular step (nxm).        
    
    """
    # Load the bench calibration file
    calibrationFileRootElement = ET.parse(fileName).getroot()
        
    # Get all the arm data container elements
    dataContainers = calibrationFileRootElement.findall("ARM_DATA_CONTAINER")
    nDataContainers = len(dataContainers)

    # Get the number of motor steps from the first data container
    slowCalTable = dataContainers[0].find("SLOW_CALIBRATION_TABLE")
        
    if slowCalTable is not None and slowCalTable.find("Joint1_fwd_stepsizes") is not None:
        # The number of motor steps is saved in the first value of the arrays
        nSteps = int(slowCalTable.find("Joint1_fwd_stepsizes").text.split(",")[0])
    else:
        nSteps = 0
    
    # Save some of the calibration information
    mids = np.zeros(nDataContainers, dtype="int")
    pids = np.zeros(nDataContainers, dtype="int")
    centers = np.zeros(nDataContainers, dtype="complex")
    tht0 = np.zeros(nDataContainers)
    tht1 = np.zeros(nDataContainers)
    phiIn = np.zeros(nDataContainers)
    phiOut = np.zeros(nDataContainers)
    pixelScale = np.zeros(nDataContainers)
    L1 = np.zeros(nDataContainers)
    L2 = np.zeros(nDataContainers)
    angularStep = np.zeros(nDataContainers)
    mapRangeTht = np.zeros((nDataContainers, 2))
    mapRangePhi = np.zeros((nDataContainers, 2))
    S1Pm = np.zeros((nDataContainers, nSteps))
    S1Nm = np.zeros((nDataContainers, nSteps))
    S2Pm = np.zeros((nDataContainers, nSteps))
    S2Nm = np.zeros((nDataContainers, nSteps))

    for i in xrange(nDataContainers):
        # Save some of the data header information
        header = dataContainers[i].find("DATA_HEADER")
        mids[i] = int(header.find("Module_Id").text)
        pids[i] = int(header.find("Positioner_Id").text)

        # Save some of the kinematics information
        kinematics = dataContainers[i].find("KINEMATICS")
        centers[i] = float(kinematics.find("Global_base_pos_x").text) + float(kinematics.find("Global_base_pos_y").text) * 1j
        tht0[i] = np.deg2rad(float(kinematics.find("CCW_Global_base_ori_z").text))
        tht1[i] = np.deg2rad(float(kinematics.find("CW_Global_base_ori_z").text))
        phiIn[i] = np.deg2rad(float(kinematics.find("Joint2_CCW_limit_angle").text)) - np.pi
        phiOut[i] = np.deg2rad(float(kinematics.find("Joint2_CW_limit_angle").text)) - np.pi
        pixelScale[i] = float(kinematics.find("Pixel_scale").text) / 1000

        if kinematics.find("Link1_Link_Length") is not None:
            L1[i] = float(kinematics.find("Link1_Link_Length").text)
        else:
            L1[i] = 25.0

        if kinematics.find("Link2_Link_Length") is not None:
            L2[i] = float(kinematics.find("Link2_Link_Length").text)
        else:
            L2[i] = 25.0

        # Save some of the slow calibration table information
        slowCalTable = dataContainers[i].find("SLOW_CALIBRATION_TABLE")

        if slowCalTable is not None and slowCalTable.find("Joint1_fwd_stepsizes") is not None:
            # Get the angular step used in the measurements
            angularPositions = slowCalTable.find("Joint1_fwd_regions").text.split(",")[2:-1]
            angularStep[i] = float(angularPositions[1]) - float(angularPositions[0])

            # Calculate the map ranges
            mapRangeTht[i] = np.deg2rad([0, nSteps * angularStep[i]])
            mapRangePhi[i] = np.deg2rad([0, nSteps * angularStep[i]]) - np.pi

            # Get the cobra motors speeds in degrees per step
            Joint1Fwd = slowCalTable.find("Joint1_fwd_stepsizes").text.split(",")[2:-1]
            Joint1Rev = slowCalTable.find("Joint1_rev_stepsizes").text.split(",")[2:-1]
            Joint2Fwd = slowCalTable.find("Joint2_fwd_stepsizes").text.split(",")[2:-1]
            Joint2Rev = slowCalTable.find("Joint2_rev_stepsizes").text.split(",")[2:-1]

            # Calculate the steps required to move that angular step
            S1Pm[i] = angularStep[i] / np.array(list(map(float, Joint1Fwd)))
            S1Nm[i] = angularStep[i] / np.array(list(map(float, Joint1Rev)))
            S2Pm[i] = angularStep[i] / np.array(list(map(float, Joint2Fwd)))
            S2Nm[i] = angularStep[i] / np.array(list(map(float, Joint2Rev)))

    # HACK SOLUTION TO BAD PHI HOMES [2016-07-15] 
    phiIn = np.maximum(-np.pi, phiIn)

    # Add some safety range to the phi limits
    phiIn = phiIn + KEEP_OUT_ANGLE
    phiOut = phiOut - KEEP_OUT_ANGLE

    # Return the calibration data as a python dictionary
    calibrationData = {}
    calibrationData["mids"] = mids
    calibrationData["pids"] = pids
    calibrationData["centers"] = centers
    calibrationData["tht0"] = tht0
    calibrationData["tht1"] = tht1
    calibrationData["phiIn"] = phiIn
    calibrationData["phiOut"] = phiOut
    calibrationData["pixelScale"] = pixelScale
    calibrationData["L1"] = L1
    calibrationData["L2"] = L2
    calibrationData["angularStep"] = angularStep
    calibrationData["mapRangeTht"] = mapRangeTht
    calibrationData["mapRangePhi"] = mapRangePhi
    calibrationData["S1Pm"] = S1Pm
    calibrationData["S1Nm"] = S1Nm
    calibrationData["S2Pm"] = S2Pm
    calibrationData["S2Nm"] = S2Nm

    return calibrationData


def defineBenchGeometry(centers, useRealMaps, useRealLinks):
    """Defines the PFI bench geometry.
    
    thetaDIR is > 0 for positive move out of home (opposite of phi)    
    thetaDIR is < 0 for negative move out of home (same as phi)
    
    Parameters
    ----------
    centers: Object
        A numpy array with the cobras central positions. If none, the cobras
        positions from a configuration file will be used.
    useRealMaps: bool
        If true, use the configuration file data for the maps.
    useRealLinks: bool
        If true, use the configuration file data for the geometry.

    Returns
    -------
    Object
        The bench object, containing the following elements:
        - center (Mx1)
        - L1     (Mx1)
        - L2     (Mx1)
        - phiIn  (Mx1)
        - phiOut (Mx1)
        - tht0   (Mx1)
        - NNmap  (MxM), logical
        - rf
        - distCobras
        - minDist
        - pids 
        - mids 
        - S1Nm 
        - S1Pm, S2Pm, S2Nm, binWidth
        
    """
    # Check if we should read the bench calibration data from an XML file
    if useRealMaps or useRealLinks or centers is None:
        calibrationData = getBenchCalibrationData("updatedMotorMapsFromThisRun2.xml")
    else:
        calibrationData = None        

    # Check if we should use the provided cobra center locations
    if centers is not None:
        # Get the total number of cobras 
        csize = len(centers)
        
        # Some relative units
        rf = 1 
        
        # Simulate some theta ranges for each cobra
        thtRange = np.deg2rad(385)
        tht0 = 2 * np.pi * np.random.random(csize)
        tht1 = (tht0 + thtRange) % (2 * np.pi)

        # Check if we should use realistic motor maps
        if useRealMaps:
            # Randomize the data container indices
            nDataContainers = calibrationData["S1Pm"].shape[0]
            indices = np.random.randint(nDataContainers, size=(4, csize))
            
            # Assign a random motor map to each cobra
            S1Pm = calibrationData["S1Pm"][indices[0]]
            S2Pm = calibrationData["S2Pm"][indices[1]]
            S1Nm = calibrationData["S1Nm"][indices[2]]
            S2Nm = calibrationData["S2Nm"][indices[3]]
        else:
            S1Pm = None
            S2Pm = None
            S1Nm = None
            S2Nm = None
        
        # Check if we should use realistic link properties
        if useRealLinks:
            # Randomize the container indices
            nDataContainers = len(calibrationData["L1"])
            indices = np.random.randint(nDataContainers, size=(4, csize))

            # Assign random link properties to each cobra
            phiIn = calibrationData["phiIn"][indices[0]] 
            phiOut = calibrationData["phiOut"][indices[1]] 
            L1 = calibrationData["L1"][indices[2]] * calibrationData["pixelScale"][indices[2]]
            L2 = calibrationData["L2"][indices[3]] * calibrationData["pixelScale"][indices[3]]
        else: 
            # Create some random values based on some realistic values
            deltaPhi = np.pi - 2 * KEEP_OUT_ANGLE
            phiOut = -KEEP_OUT_ANGLE * (1.0 + 0.2 * np.random.random(csize))
            phiIn = phiOut - deltaPhi
            L1 = np.full(csize, LINK_LENGTH)
            L2 = np.full(csize, LINK_LENGTH)
    else:
        # Use the calibration data
        centers = calibrationData["centers"]
        tht0 = calibrationData["tht0"]
        tht1 = calibrationData["tht1"]
        phiIn = calibrationData["phiIn"]
        phiOut = calibrationData["phiOut"]
        L1 = calibrationData["L1"]
        L2 = calibrationData["L2"]       
        S1Pm = calibrationData["S1Pm"]
        S1Nm = calibrationData["S1Nm"]
        S2Pm = calibrationData["S2Pm"]
        S2Nm = calibrationData["S2Nm"]
        
        # Some relative units
        rf = 1.0 / calibrationData["pixelScale"][0]
    
    # Calculate the maximum and minimum radius that the cobras can reach
    rMin = np.abs(L1 + L2 * np.exp(1j * phiIn))
    rMax = np.abs(L1 + L2 * np.exp(1j * phiOut))
    
    # Parameters for populating the annular patrol areas uniformly
    dA = 1.0 / ((rMax / rMin) ** 2 - 1.0)  # fractional keepout area
    rRange = np.sqrt(rMax ** 2 - rMin ** 2)

    # By default phi home is phi in.  
    # Existence of mat files *may* change phi home.
    phiHome = phiIn.copy()

    # For the realities of the test bench, having phi too far in is
    # inconvenient.  Use a larger value so that the theta arm does not
    # go crazy.
    phiHome += 0.5

    # Read in phiHome from elsewhere
    # TBD (See matlab code)
        
    # Home positions (0 = ss, 1 = os)
    home0 = centers + L1 * np.exp(1j * tht0) + L2 * np.exp(1j * (tht0 + phiHome))
    home1 = centers + L1 * np.exp(1j * tht1) + L2 * np.exp(1j * (tht1 + phiHome))

    # Calculate the theta overlap
    thtOverlap = ((tht1 - tht0 + np.pi) % (2 * np.pi)) - np.pi

    # Generate the nearest neighbor map/structure
    csize = len(centers)
    CCdist = np.zeros((csize, csize))

    for i in xrange(csize):
        CCdist[i] = np.abs(centers - centers[i])
    
    distCenters = COBRAS_SEPARATION * (2 * np.median(L1 + L2)) / (2 * LINK_LENGTH)
    nnMap = np.logical_and(CCdist > 1e-9, CCdist < 1.5 * distCenters)
    (row, col) = np.where(nnMap)
    
    NN = {}
    NN["row"] = row
    NN["col"] = col
    NN["xy"] = (centers[row] + centers[col]) / 2

    # Calculate some parameters for the field geometry
    field = {}
    field["cm"] = np.mean(centers)
    field["R"] = np.max(np.abs(centers - field["cm"]) + rMax)
    
    # Fit a line to the centers and calculate the slope
    x = np.real(centers)
    y = np.imag(centers)
    slope = (csize * (x * y).sum() - x.sum() * y.sum()) / (csize * (x * x).sum() - x.sum() ** 2);
    field["ang"] = np.arctan([slope])

    # Save the bench results
    bench = {}
    bench["center"] = centers
    bench["tht0"] = tht0
    bench["tht1"] = tht1
    bench["phiIn"] = phiIn
    bench["phiOut"] = phiOut
    bench["L1"] = L1
    bench["L2"] = L2
    bench["S1Pm"] = S1Pm
    bench["S1Nm"] = S1Nm
    bench["S2Pm"] = S2Pm
    bench["S2Nm"] = S2Nm
    bench["rf"] = rf
    bench["distCobras"] = COBRAS_SEPARATION * rf
    bench["minDist"] = 2 * rf
    bench["rMin"] = rMin
    bench["rMax"] = rMax
    bench["dA"] = dA
    bench["rRange"] = rRange
    bench["home0"] = home0
    bench["home1"] = home1
    bench["thtOverlap"] = thtOverlap
    bench["nnMap"] = nnMap
    bench["NN"] = NN
    bench["field"] = field
    bench["binWidth"] = 2 * np.pi / 100
    bench["alpha"] = ALPHA
    bench["beta"] = BETA
    bench["thteps"] = THT_EPS

    if calibrationData is not None:
        bench["mapRangeTht"] = calibrationData["mapRangeTht"]
        bench["mapRangePhi"] = calibrationData["mapRangePhi"]
        bench["pids"] = calibrationData["pids"]
        bench["mids"] = calibrationData["mids"]
    else:
        bench["mapRangeTht"] = None
        bench["mapRangePhi"] = None
        bench["pids"] = None
        bench["mids"] = None
        
    return bench


if __name__ == "__main__":
    # Get the bench
    bench = defineBenchGeometry(None, 1, 1)
    
    # Print the data in the console
    for key in bench:
        print("Bench data: " + key)
        print(bench[key])
    
    print(bench["home0"])
