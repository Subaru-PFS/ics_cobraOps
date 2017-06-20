"""

Some utility methods related with the PFS bench.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T
  
"""

import numpy as np
import xml.etree.ElementTree as ET

import plotUtils as plotUtils


KEEP_OUT_ANGLE = 0.1
"""Safety angle value in radians (~ 5.7 deg)."""

ALPHA = 0.07
"""The motor noise normalization factor."""

BETA = 0.50
"""The motor noise exponent."""

LINK_LENGTH = 2.375
"""The theoretical link length in mm."""
 
COBRAS_SEPARATION = 8.0
"""The theoretical separation between two consecutive cobras in mm."""
           
THT_EPS = 2e-10
"""This is an epsilon to make sure that values near zero in theta are
interpreted on the positive (or negative) side of the cut, respectively. 
Make it negative for same same direction moveouts (positive for 
positive hardstop).

tht1 implemented below is an opposite-sense hard stop.  when that's 
working, theta direction will be specified on a per cobra basis, not 
on a full instrument basis, as is implemented here."""

BIN_WIDTH = 2 * np.pi / 100
"""Angular bin width."""

BENCH_CALIBRATION_FILE = "updatedMotorMapsFromThisRun2.xml"
"""The path to the bench calibration XML file."""


def getBenchCalibrationData(fileName):
    """Returns some bench calibration data stored in a XML file.

    Parameters
    ----------
    fileName: str
        The path to the XML calibration file.
    
    Returns
    -------
    object
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
        - mapRangeTht: The theta angle map range in radians (nx2).
        - mapRangePhi: The phi angle map range in radians (nx2).
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
    mids = np.empty(nDataContainers, dtype="int")
    pids = np.empty(nDataContainers, dtype="int")
    centers = np.empty(nDataContainers, dtype="complex")
    tht0 = np.empty(nDataContainers)
    tht1 = np.empty(nDataContainers)
    phiIn = np.empty(nDataContainers)
    phiOut = np.empty(nDataContainers)
    pixelScale = np.empty(nDataContainers)
    L1 = np.empty(nDataContainers)
    L2 = np.empty(nDataContainers)
    angularStep = np.zeros(nDataContainers)
    mapRangeTht = np.zeros((nDataContainers, 2))
    mapRangePhi = np.zeros((nDataContainers, 2))
    S1Pm = np.empty((nDataContainers, nSteps))
    S1Nm = np.empty((nDataContainers, nSteps))
    S2Pm = np.empty((nDataContainers, nSteps))
    S2Nm = np.empty((nDataContainers, nSteps))

    for i in range(nDataContainers):
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
        pixelScale[i] = float(kinematics.find("Pixel_scale").text) / 1000.0

        if kinematics.find("Link1_Link_Length") is not None:
            L1[i] = float(kinematics.find("Link1_Link_Length").text)
        else:
            L1[i] = 25.0

        if kinematics.find("Link2_Link_Length") is not None:
            L2[i] = float(kinematics.find("Link2_Link_Length").text)
        else:
            L2[i] = 25.0

        # Save the motor information
        if nSteps > 0:
            # Get the angular step used in the measurements
            slowCalTable = dataContainers[i].find("SLOW_CALIBRATION_TABLE")
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
    centers: object
        A numpy array with the cobras central positions. If none, the cobras
        positions from a configuration file will be used.
    useRealMaps: bool
        If true, the cobra motor maps will be loaded from a calibration file.
        No motor maps will be used otherwise.
    useRealLinks: bool
        If true, the cobra link properties will be loaded from a calibration 
        file. Some approximations will be used otherwise.

    Returns
    -------
    object
        The bench geometry object, containing the following elements:
        - center: The cobra central positions in mm or pixel units (n).
        - tht0: The hard stop angle for same sense move-out in radians (n).
        - tht1: The hard stop angle for opposite sense move-out in radians (n).
        - phiIn: ... in radians (n).
        - phiOut: ... in radians (n).
        - L1: The link1 lengths in mm or pixel units (n).
        - L2: The link2 lengths in mm or pixel units (n).
        - S1Pm: The joint1 forward motor steps required to move an angular step (nxm).
        - S1Nm: The joint1 reverse motor steps required to move an angular step (nxm).
        - S2Pm: The joint2 forward motor steps required to move an angular step (nxm).
        - S2Nm: The joint2 reverse motor steps required to move an angular step (nxm).        
        - rf: The length transformation factor from mm to the L1 and L2 length units (n).
        - distCobras: The separation between two consecutive cobras in length units (n).
        - minDist: 2 times rf (n).
        - rMin: The minimum radius that the cobra can reach in length units (n).
        - rMax: The maximum radius that the cobra can reach in length units (n).
        - home0: The home position for the same sense move-out in length units (n).
        - home1: The home position for the opposite sense move-out in length units (n).
        - thtOverlap: The theta overlap ??? (n)
        - NN: A more useful way to the express the the nearest neighbor map.
        - dA: Parameter for populating the annular patrol areas uniformly (n).
        - rRange: Parameter for populating the annular patrol areas uniformly (n).
        - alpha: The motor noise normalization factor.
        - beta: The motor noise exponent.
        - thteps: Theta epsilon value.
        - binWidth: Angular bind width.
        - mids: The module ids from the calibration file (c).
        - pids: The positioner ids from the calibration file (c).
        - mapRangeTht: The theta angle map range in radians (cx2).
        - mapRangePhi: The phi angle map range in radians (cx2).
        
    """
    # Check if we should read the bench calibration data from an XML file
    if useRealMaps or useRealLinks or centers is None:
        calibrationData = getBenchCalibrationData(BENCH_CALIBRATION_FILE)
    else:
        calibrationData = None        

    # Check if we should use the provided cobra center locations
    if centers is not None:
        # Get the total number of cobras 
        nCobras = len(centers)
        
        # The centers should be in millimeters, so the length transformation 
        # factor is 1
        lengthTransformationFactor = np.ones(nCobras)
        
        # Simulate some theta ranges for each cobra
        thtRange = np.deg2rad(385.0)
        tht0 = 2 * np.pi * np.random.random(nCobras)
        tht1 = (tht0 + thtRange) % (2 * np.pi)
        
        # Check if we should use realistic link properties
        if useRealLinks:
            # Randomize the container indices
            nDataContainers = len(calibrationData["L1"])
            indices = np.random.randint(nDataContainers, size=(4, nCobras))

            # Assign random link properties to each cobra
            phiIn = calibrationData["phiIn"][indices[0]] 
            phiOut = calibrationData["phiOut"][indices[1]] 
            L1 = calibrationData["L1"][indices[2]] * calibrationData["pixelScale"][indices[2]]
            L2 = calibrationData["L2"][indices[3]] * calibrationData["pixelScale"][indices[3]]
        else: 
            # Create some random values based on some realistic values
            deltaPhi = np.pi - 2 * KEEP_OUT_ANGLE
            phiOut = -KEEP_OUT_ANGLE * (1.0 + 0.2 * np.random.random(nCobras))
            phiIn = phiOut - deltaPhi
            L1 = np.full(nCobras, LINK_LENGTH)
            L2 = np.full(nCobras, LINK_LENGTH)

        # Check if we should use realistic motor maps
        if useRealMaps:
            # Randomize the data container indices
            nDataContainers = calibrationData["S1Pm"].shape[0]
            indices = np.random.randint(nDataContainers, size=(4, nCobras))
            
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
        
        # The centers are in pixel units, so the length transformation 
        # factor is the inverse of the pixel scale
        lengthTransformationFactor = 1.0 / calibrationData["pixelScale"]
    
    # Calculate the maximum and minimum radius that the cobras can reach
    rMin = np.abs(L1 + L2 * np.exp(1j * phiIn))
    rMax = np.abs(L1 + L2 * np.exp(1j * phiOut))

    # By default phi home is phiIn
    phiHome = phiIn.copy()

    # For the realities of the test bench, having phi too far in is inconvenient.  
    # Use a larger value so that the theta arm does not go crazy.
    phiHome += 0.5

    # Read in phiHome from elsewhere
    # TBD (See MATLAB code)
        
    # Calculate the home positions (0 = same sense, 1 = opposite sense)
    home0 = centers + L1 * np.exp(1j * tht0) + L2 * np.exp(1j * (tht0 + phiHome))
    home1 = centers + L1 * np.exp(1j * tht1) + L2 * np.exp(1j * (tht1 + phiHome))

    # Calculate the theta overlap
    thtOverlap = ((tht1 - tht0 + np.pi) % (2 * np.pi)) - np.pi

    # Calculate the cobras distance matrix
    distanceMatrix = np.abs(centers[:, np.newaxis] - centers)
    
    # Obtain the nearest neighbors map 
    distCenters = COBRAS_SEPARATION * np.median(L1 + L2) / (2 * LINK_LENGTH)
    nnMap = np.logical_and(distanceMatrix > 1e-9, distanceMatrix < 1.5 * distCenters)

    # Save the nearest neighbors results on a more useful structure  
    (row, col) = np.where(nnMap)
    NN = {}
    NN["row"] = row
    NN["col"] = col
    NN["xy"] = (centers[row] + centers[col]) / 2
    
    # Calculate some parameters for populating the annular patrol areas uniformly
    dA = 1.0 / ((rMax / rMin) ** 2 - 1.0) 
    rRange = np.sqrt(rMax ** 2 - rMin ** 2)

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
    bench["rf"] = lengthTransformationFactor
    bench["distCobras"] = COBRAS_SEPARATION * lengthTransformationFactor
    bench["minDist"] = 2 * lengthTransformationFactor
    bench["rMin"] = rMin
    bench["rMax"] = rMax
    bench["home0"] = home0
    bench["home1"] = home1
    bench["thtOverlap"] = thtOverlap
    bench["NN"] = NN
    bench["dA"] = dA
    bench["rRange"] = rRange
    bench["alpha"] = ALPHA
    bench["beta"] = BETA
    bench["thteps"] = THT_EPS
    bench["binWidth"] = BIN_WIDTH

    if calibrationData is not None:
        bench["mids"] = calibrationData["mids"]
        bench["pids"] = calibrationData["pids"]
        bench["mapRangeTht"] = calibrationData["mapRangeTht"]
        bench["mapRangePhi"] = calibrationData["mapRangePhi"]
    else:
        bench["mids"] = None
        bench["pids"] = None
        bench["mapRangeTht"] = None
        bench["mapRangePhi"] = None
        
    return bench


def plotBenchGeometry(bench):
    """Plots the bench geometry.

    Parameters
    ----------
    bench: object
        The bench object. 
    
    """
    # Create the figure
    plotUtils.createNewFigure("Bench geometry", "x position", "y position")

    # Set the axes limits
    benchCenter = np.mean(bench["center"])
    benchRadius = np.max(np.abs(bench["center"] - benchCenter) + bench["rMax"])
    limRange = 1.05 * benchRadius * np.array([-1, 1]) 
    xLim = benchCenter.real + limRange
    yLim = benchCenter.imag + limRange
    plotUtils.setAxesLimits(xLim, yLim)
    
    # Plot the cobra patrol areas using ring shapes
    plotUtils.addRings(bench["center"], bench["rMin"], bench["rMax"], facecolors="blue", alpha=0.15)

    # add the stage 1 theta hard stops.
    plotUtils.addLines(bench['center'], bench['center'] + bench['rMax'] * np.exp(1j * bench['tht0']), linestyles='dashed' )    
    plotUtils.addLines(bench['center'], bench['center'] + bench['rMax'] * np.exp(1j * bench['tht1']), linestyles='dashdot')

if __name__ == "__main__":
    # Get the bench from the calibration file
    bench = defineBenchGeometry(None, 1, 1)
        
    # Print the data in the console
    for key in bench:
        print("Bench data: " + key)
        print(bench[key])
    
    # Plot the bench geomtry
    plotBenchGeometry(bench)
    plotUtils.pauseExecution()

