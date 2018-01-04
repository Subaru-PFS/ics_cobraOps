"""

Cobras calibration product class.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np
import xml.etree.ElementTree as ElementTree

from ics.cobraOps.AttributePrinter import AttributePrinter


class CobrasCalibrationProduct(AttributePrinter):
    """
    
    Class describing a cobras calibration product.
    
    """
    
    def __init__(self, fileName):
        """Constructs a new cobras calibration product using the information
        contained in an XML calibration file.
        
        Parameters
        ----------
        fileName: object
            The path to the XML calibration file.
        
        Returns
        -------
        object
            The cobras calibration product.
        
        """
        # Load the XML calibration file
        calibrationFileRootElement = ElementTree.parse(fileName).getroot()
        
        # Get all the data container elements
        dataContainers = calibrationFileRootElement.findall("ARM_DATA_CONTAINER")
        
        # The number of cobras is equal to the number of data containers
        self.nCobras = len(dataContainers)
        
        # Create some of the calibration data arrays
        self.moduleIds = np.empty(self.nCobras, dtype="int")
        self.positionerIds = np.empty(self.nCobras, dtype="int")
        self.centers = np.empty(self.nCobras, dtype="complex")
        self.tht0 = np.empty(self.nCobras)
        self.tht1 = np.empty(self.nCobras)
        self.phiIn = np.empty(self.nCobras)
        self.phiOut = np.empty(self.nCobras)
        self.L1 = np.empty(self.nCobras)
        self.L2 = np.empty(self.nCobras)
        
        # Check if the data containers have information about the motor maps
        slowCalTable = dataContainers[0].find("SLOW_CALIBRATION_TABLE")
        
        if slowCalTable is not None and slowCalTable.find("Joint1_fwd_stepsizes") is not None:
            # The number of motor map steps is saved in the first element of
            # the arrays
            self.motorMapSteps = int(slowCalTable.find("Joint1_fwd_stepsizes").text.split(",")[0])
            
            # Create the cobra motor map arrays
            self.angularSteps = np.empty(self.nCobras)
            self.S1Pm = np.empty((self.nCobras, self.motorMapSteps))
            self.S1Nm = np.empty((self.nCobras, self.motorMapSteps))
            self.S2Pm = np.empty((self.nCobras, self.motorMapSteps))
            self.S2Nm = np.empty((self.nCobras, self.motorMapSteps))
            self.F1Pm = np.empty((self.nCobras, self.motorMapSteps))
            self.F1Nm = np.empty((self.nCobras, self.motorMapSteps))
            self.F2Pm = np.empty((self.nCobras, self.motorMapSteps))
            self.F2Nm = np.empty((self.nCobras, self.motorMapSteps))
        
        # Fill the cobras calibration arrays
        for i in range(self.nCobras):
            # Save some of the data header information
            header = dataContainers[i].find("DATA_HEADER")
            self.moduleIds[i] = int(header.find("Module_Id").text)
            self.positionerIds[i] = int(header.find("Positioner_Id").text)
            
            # Save some of the kinematics information
            kinematics = dataContainers[i].find("KINEMATICS")
            self.centers[i] = float(kinematics.find("Global_base_pos_x").text) + float(kinematics.find("Global_base_pos_y").text) * 1j
            self.tht0[i] = np.deg2rad(float(kinematics.find("CCW_Global_base_ori_z").text))
            self.tht1[i] = np.deg2rad(float(kinematics.find("CW_Global_base_ori_z").text))
            self.phiIn[i] = np.deg2rad(float(kinematics.find("Joint2_CCW_limit_angle").text)) - np.pi
            self.phiOut[i] = np.deg2rad(float(kinematics.find("Joint2_CW_limit_angle").text)) - np.pi
            self.L1[i] = float(kinematics.find("Link1_Link_Length").text)
            self.L2[i] = float(kinematics.find("Link2_Link_Length").text)
            
            # Transform the length properties from pixel units to mm
            pixelScale = float(kinematics.find("Pixel_scale").text) / 1000.0
            self.centers[i] *= pixelScale
            self.L1[i] *= pixelScale
            self.L2[i] *= pixelScale
            
            # Save the motor calibration information
            if hasattr(self, "motorMapSteps"):
                # Get the angular step used in the measurements
                slowCalTable = dataContainers[i].find("SLOW_CALIBRATION_TABLE")
                fastCalTable = dataContainers[i].find("FAST_CALIBRATION_TABLE")
                angularPositions = slowCalTable.find("Joint1_fwd_regions").text.split(",")[2:-1]
                angularStep = float(angularPositions[1]) - float(angularPositions[0])
                
                # Get the cobra motors speeds in degrees per step
                slowJoint1Fwd = slowCalTable.find("Joint1_fwd_stepsizes").text.split(",")[2:-1]
                slowJoint1Rev = slowCalTable.find("Joint1_rev_stepsizes").text.split(",")[2:-1]
                slowJoint2Fwd = slowCalTable.find("Joint2_fwd_stepsizes").text.split(",")[2:-1]
                slowJoint2Rev = slowCalTable.find("Joint2_rev_stepsizes").text.split(",")[2:-1]
                fastJoint1Fwd = fastCalTable.find("Joint1_fwd_stepsizes").text.split(",")[2:-1]
                fastJoint1Rev = fastCalTable.find("Joint1_rev_stepsizes").text.split(",")[2:-1]
                fastJoint2Fwd = fastCalTable.find("Joint2_fwd_stepsizes").text.split(",")[2:-1]
                fastJoint2Rev = fastCalTable.find("Joint2_rev_stepsizes").text.split(",")[2:-1]
                
                # Calculate the motor steps required to move that angular step
                self.S1Pm[i] = angularStep / np.array(list(map(float, slowJoint1Fwd)))
                self.S1Nm[i] = angularStep / np.array(list(map(float, slowJoint1Rev)))
                self.S2Pm[i] = angularStep / np.array(list(map(float, slowJoint2Fwd)))
                self.S2Nm[i] = angularStep / np.array(list(map(float, slowJoint2Rev)))
                self.F1Pm[i] = angularStep / np.array(list(map(float, fastJoint1Fwd)))
                self.F1Nm[i] = angularStep / np.array(list(map(float, fastJoint1Rev)))
                self.F2Pm[i] = angularStep / np.array(list(map(float, fastJoint2Fwd)))
                self.F2Nm[i] = angularStep / np.array(list(map(float, fastJoint2Rev)))
                
                # Save the angular step in radians
                self.angularSteps[i] = np.deg2rad(angularStep)
        
        # HACK SOLUTION TO BAD PHI HOMES
        self.phiIn = np.maximum(-np.pi, self.phiIn)
        self.phiOut = np.minimum(self.phiOut, 0.0)
