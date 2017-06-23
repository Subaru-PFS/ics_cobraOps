"""

Some utility methods related with the cobra trajectories.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

import numpy as np
import time as time

import cobraUtils as cobraUtils
import benchUtils as benchUtils
import targetUtils as targetUtils
import plotUtils as plotUtils

           
THT_EPS = 1e-12
"""This is an epsilon to make sure that values near zero in theta are
interpreted on the positive (or negative) side of the cut, respectively. 
Make it negative for same same direction moveouts (positive for 
positive hardstop)."""


def generateTrajectories(finalPositions, bench):
    """Generates the cobra trajectories starting from their home positions.

    Parameters
    ----------
    finalPositions: object
        A complex numpy array with the cobras final positions.
    bench: object
        The bench geometry to use.

    Returns
    -------
    Object
        ...

    """
    # Extract some of the bench geometry information
    centers = bench["center"]
    L1 = bench["L1"]
    L2 = bench["L2"]
    home0 = bench["home0"]
    home1 = bench["home1"]
    tht0 = bench["tht0"]
    nCobras = len(centers)

    # Get the cobra rotation angles for the initial home positions and the final positions
    (startThtP, startPhiP) = targetUtils.getCobraRotationAngles(home0 - centers, L1, L2)
    (startThtN, startPhiN) = targetUtils.getCobraRotationAngles(home1 - centers, L1, L2)
    (finishTht, finishPhi) = targetUtils.getCobraRotationAngles(finalPositions - centers, L1, L2)

    # Check if the bench geometry contains the cobra motor maps information
    if bench["S1Pm"] is not None:
        # Extract the cobra motor maps, reversing the index order for the negative movement maps
        mapThtP = bench["S1Pm"]
        mapThtN = np.fliplr(bench["S1Nm"])
        mapPhiP = bench["S2Pm"]
        mapPhiN = np.fliplr(bench["S2Nm"])

        # Get the number of motor step bins in each cobra rotation angle
        nBinsTht = mapThtP.shape[1]
        nBinsPhi = mapPhiP.shape[1]

        # Get the bin step size in radians
        binWidth = bench["binWidth"]

        # Check if we are starting from the overlap region for negative movements
        STARTS_IN_OVERLAP = np.mod(startThtN - tht0, 2 * np.pi) < bench["thtOverlap"] + THT_EPS

        # Calculate the starting bins for the positive and negative movements
        startBinThtP = (np.mod(startThtP + THT_EPS - tht0, 2 * np.pi) - THT_EPS - bench["mapRangeTht"][0, 0]) / binWidth
        startBinThtN = (bench["mapRangeTht"][0, 1] - (np.mod(startThtN - tht0, 2 * np.pi) + 2 * np.pi * STARTS_IN_OVERLAP)) / binWidth
        startBinPhiP = (startPhiP - bench["mapRangePhi"][0, 0]) / binWidth
        startBinPhiN = (nBinsPhi - 1) - startBinPhiP
        
        # Calculate the starting indices
        startIndThtP = np.floor(startBinThtP).astype(int)
        startIndThtN = np.floor(startBinThtN).astype(int)
        startIndPhiP = np.floor(startBinPhiP).astype(int)
        startIndPhiN = np.floor(startBinPhiN).astype(int)
        startIndThtP[startIndThtP < 0] = 0
        startIndThtN[startIndThtN < 0] = 0
        startIndPhiP[startIndPhiP < 0] = 0
        startIndPhiN[startIndPhiN < 0] = 0

        #  Check if we are finishing at the overlap region for negative movements
        FINISHES_IN_OVERLAP = np.mod(finishTht + THT_EPS - tht0, 2 * np.pi) < bench["thtOverlap"]
            
        # Calculate the finishing bins for the positive and negative moves
        finishBinThtP = (np.mod(finishTht + THT_EPS - tht0, 2 * np.pi) - THT_EPS - bench["mapRangeTht"][0, 0]) / binWidth
        finishBinThtN = (nBinsTht - 1) - (finishBinThtP + 2 * np.pi * FINISHES_IN_OVERLAP / binWidth)
        finishBinPhiP = (finishPhi - bench["mapRangePhi"][0, 0]) / binWidth
        finishBinPhiN = (nBinsPhi - 1) - finishBinPhiP

        # Calculate the finishing indices
        finishIndThtP = np.floor(finishBinThtP).astype(int)
        finishIndThtN = np.floor(finishBinThtN).astype(int)
        finishIndPhiP = np.floor(finishBinPhiP).astype(int)
        finishIndPhiN = np.floor(finishBinPhiN).astype(int)
        finishIndThtP[finishIndThtP < 0] = 0 
        finishIndThtN[finishIndThtN < 0] = 0 
        finishIndPhiP[finishIndPhiP < 0] = 0 
        finishIndPhiN[finishIndPhiN < 0] = 0 

        # Calculate the residuals at the starting and finishing bin
        startOverCountThtP = startBinThtP - startIndThtP
        startOverCountThtN = startBinThtN - startIndThtN
        startOverCountPhiP = startBinPhiP - startIndPhiP
        startOverCountPhiN = startBinPhiN - startIndPhiN
        finishOverCountThtP = finishIndThtP + 1 - finishBinThtP
        finishOverCountThtN = finishIndThtN + 1 - finishBinThtN
        finishOverCountPhiP = finishIndPhiP + 1 - finishBinPhiP
        finishOverCountPhiN = finishIndPhiN + 1 - finishBinPhiN

        # Calculate the number of steps required to reach the finish position
        nStepsThtP = np.zeros(nCobras)
        nStepsThtN = np.zeros(nCobras)
        nStepsPhiP = np.zeros(nCobras)
        nStepsPhiN = np.zeros(nCobras)
        
        for i in range(nCobras):
            nStepsThtP[i] = np.sum(mapThtP[i, startIndThtP[i]:finishIndThtP[i] + 1]) \
                - mapThtP[i, startIndThtP[i]] * startOverCountThtP[i] \
                - mapThtP[i, finishIndThtP[i]] * finishOverCountThtP[i]
            nStepsThtN[i] = np.sum(mapThtN[i, startIndThtN[i]:finishIndThtN[i] + 1]) \
                - mapThtN[i, startIndThtN[i]] * startOverCountThtN[i] \
                - mapThtN[i, finishIndThtN[i]] * finishOverCountThtN[i]
            nStepsPhiP[i] = np.sum(mapPhiP[i, startIndPhiP[i]:finishIndPhiP[i] + 1]) \
                - mapPhiP[i, startIndPhiP[i]] * startOverCountPhiP[i] \
                - mapPhiP[i, finishIndPhiP[i]] * finishOverCountPhiP[i]
            nStepsPhiN[i] = np.sum(mapPhiN[i, startIndPhiN[i]:finishIndPhiN[i] + 1]) \
                - mapPhiN[i, startIndPhiN[i]] * startOverCountPhiN[i] \
                - mapPhiN[i, finishIndPhiN[i]] * finishOverCountPhiN[i]

        nStepsThtP[nStepsThtP < 0] = 0
        nStepsThtN[nStepsThtN < 0] = 0
        nStepsPhiP[nStepsPhiP < 0] = 0
        nStepsPhiN[nStepsPhiN < 0] = 0
        
        # Calculate the maximum number of steps for each cobra
        nStepsMax = np.max(np.vstack((nStepsThtP, nStepsThtN, nStepsPhiP, nStepsPhiN)), axis=0)
    
        # Add some error to the calibration motor maps
        fractionalBinError = bench["alpha"] * binWidth ** (bench["beta"] - 1)
        noisyMapThtP = mapThtP * calculateMapErrorFactor(fractionalBinError, mapThtP.shape)
        noisyMapThtN = mapThtN * calculateMapErrorFactor(fractionalBinError, mapThtN.shape)
        noisyMapPhiP = mapPhiP * calculateMapErrorFactor(fractionalBinError, mapPhiP.shape)
        noisyMapPhiN = mapPhiN * calculateMapErrorFactor(fractionalBinError, mapPhiN.shape)
    
        # Add an sticky bin at the maps end to prevent over-runs
        noisyMapThtP = np.hstack((noisyMapThtP, 1e9 * np.ones((nCobras, 1))))
        noisyMapThtN = np.hstack((noisyMapThtN, 1e9 * np.ones((nCobras, 1))))
        noisyMapPhiP = np.hstack((noisyMapPhiP, 1e9 * np.ones((nCobras, 1))))
        noisyMapPhiN = np.hstack((noisyMapPhiN, 1e9 * np.ones((nCobras, 1))))
        
        # Calculate the cobra trajectories
        stepsPerTime = 50
        maxStepsBins = np.ceil(np.max(nStepsMax) / stepsPerTime).astype(int) + 1
        thtP = np.empty((nCobras, maxStepsBins))
        thtN = np.empty((nCobras, maxStepsBins))
        phiP = np.empty((nCobras, maxStepsBins))
        phiN = np.empty((nCobras, maxStepsBins))
        trajBinsThtP = np.empty(nCobras, dtype="int")      
        trajBinsThtN = np.empty(nCobras, dtype="int")      
        trajBinsPhiP = np.empty(nCobras, dtype="int")      
        trajBinsPhiN = np.empty(nCobras, dtype="int")      
        
        for i in range(10):
            # THT P (tht moving out)
            stepOffset = noisyMapThtP[i, startIndThtP[i]] * startOverCountThtP[i]            
            finishInd = np.where(nStepsThtP[i] < np.cumsum(noisyMapThtP[i, startIndThtP[i]:]) - stepOffset)[0][0] + startIndThtP[i]
            stepCtr = np.append([0], np.cumsum(noisyMapThtP[i, startIndThtP[i]:finishInd + 1]) - stepOffset)
            binCtr = np.append([startBinThtP[i]], startIndThtP[i] + np.arange(1, len(stepCtr)))
            trajSteps = np.append(np.arange(0, nStepsThtP[i], stepsPerTime), [nStepsThtP[i]])
            thtP[i][0:len(trajSteps)] = tht0[i] + binWidth * np.interp(trajSteps, stepCtr, binCtr)
            trajBinsThtP[i] = len(trajSteps)

            # THT N (opposite sense)
            stepOffset = noisyMapThtN[i, startIndThtN[i]] * startOverCountThtN[i]            
            finishInd = np.where(nStepsThtN[i] < np.cumsum(noisyMapThtN[i, startIndThtN[i]:]) - stepOffset)[0][0] + startIndThtN[i]
            stepCtr = np.append([0], np.cumsum(noisyMapThtN[i, startIndThtN[i]:finishInd + 1]) - stepOffset)
            binCtr = np.append([startBinThtN[i]], startIndThtN[i] + np.arange(1, len(stepCtr)))
            trajSteps = np.append(np.arange(0, nStepsThtN[i], stepsPerTime), [nStepsThtN[i]])
            thtN[i][0:len(trajSteps)] = tht0[i] + bench["mapRangeTht"][0, 1] - binWidth * np.interp(trajSteps, stepCtr, binCtr)
            trajBinsThtN[i] = len(trajSteps)
            # print(nStepsThtP[i])
            # print(nStepsThtN[i])
            # print(str(tht0[i]) + " --> " + str(finishTht[i]) + " " + str(finishTht[i] % (2 * np.pi)))
            # print("a " + str(thtP[i, trajBinsThtP[i] - 1]) + " " + str(thtP[i, trajBinsThtP[i] - 1] % (2 * np.pi)))
            # print("a " + str(thtN[i, trajBinsThtN[i] - 1]) + " " + str(thtN[i, trajBinsThtN[i] - 1] % (2 * np.pi)))
            # print("b " + str(finishTht[i]) + " " + str(finishTht[i]% (2*np.pi)))
            
            # Phi P (tht moving out)
            """
            stepOffset = noisyMapThtP[i, startIndThtP[i]] * startOverCountThtP[i]            
            finishInd = np.where(nStepsThtP[i] < np.cumsum(noisyMapThtP[i, startIndThtP[i]:]) - stepOffset)[0][0] + startIndThtP[i]
            stepCtr = np.append([0], np.cumsum(noisyMapThtP[i, startIndThtP[i]:finishInd + 1]) - stepOffset)
            binCtr = np.append([startBinThtP[i]], startIndThtP[i] + np.arange(1, len(stepCtr)))
            trajSteps = np.append(np.arange(0, nStepsThtP[i], stepsPerTime), [nStepsThtP[i]])
            thtP[i][0:len(trajSteps)] = startThtP[i] + binWidth * np.interp(trajSteps, stepCtr, binCtr)
            trajBinsThtP[i] = len(trajSteps)
            """


def calculateMapErrorFactor(fractionalBinError, mapShape):
    """Calculates the map error factor.

    Parameters
    ----------
    fractionalBinError: float
        The fractional bin error.
    mapShape: tuple
        The map shape (nCobras x nBins).

    Returns
    -------
    Object
        Numpy array with the calculated map error factor.

    """
    # Initialize the position and the map factor arrays
    nBins = mapShape[0] * mapShape[1]
    x = np.zeros(nBins)
    factor = np.zeros(nBins)
    
    # Loop until all the bins have reached a position value equal or larger than one
    toMove = x < 1
     
    while np.any(toMove):
        # Calculate the movement for each bin that we still need to move
        if fractionalBinError == 0:
            dx = np.ones(np.sum(toMove))
        else:
            dx = np.random.normal(loc=1, scale=fractionalBinError, size=np.sum(toMove))
        
        # Update the factor array
        factor[toMove] += np.clip((1.0 - x[toMove]) / dx, 0, 1)

        # Update the position array
        x[toMove] += dx
        
        # Recalculate the bins that we need to move in the next loop step
        toMove = x < 1

    # Change the factor array shape to the input map shape
    factor.shape = mapShape
    
    return factor



if __name__ == "__main__":
    # Define the target density to use
    targetDensity = 2
    
    # Get the cobras central positions for the full PFI
    start = time.time()
    centers = cobraUtils.getPFICenters()
    print("Number of cobras: " + str(len(centers)))

    # Define the bench geometry
    bench = benchUtils.defineBenchGeometry(centers, True, True)
    
    # Create a random sample of targets
    targetPositions = targetUtils.generateTargets(targetDensity, bench)
    print("Number of simulated targets: " + str(len(targetPositions)))

    # Assign the target to the cobras and get the cobra positions
    (assignedTargets, cobraPositions) = targetUtils.assignTargets(targetPositions, bench)
    
    # Get the cobras for which the collision could not solved
    (problematicCobras, nearbyProblematicCobras) = targetUtils.getProblematicCobras(cobraPositions, bench)
    print("Number of unsolved collisions: " + str(len(problematicCobras) / 2))
    print("Total computation time (s): " + str(time.time() - start))
    
    # Generate the trajectories
    generateTrajectories(cobraPositions, bench) 

    # Plot the cobra-target associations
    start = time.time()
    targetUtils.plotCobraTargetAssociations(cobraPositions, problematicCobras, assignedTargets, targetPositions, bench)
    print("Plotting time (s): " + str(time.time() - start))
    plotUtils.pauseExecution()


