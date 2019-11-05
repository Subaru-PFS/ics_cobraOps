# Main modules

## cobraConstants.py

Defines some constants related with the cobras and the simulation code, like the default cobra link lengths, the number of cobras per module, or the default motor maps angular step size.

## CobrasCalibrationProduct.py

Defines the `CobrasCalibrationProduct` class. This class is used to read (and save for testing) the XML files with the cobras calibration information: cobra centers, links lengths, motor maps, etc.

## MotorMapGoup.py

Defines the `MotorMapGroup` class. This class is used to represent the motor map properties of a given group of cobras. These properties can be taken from a `CobrasCalibrationProduct` instance, or can be set to defaults: constant motor maps.

The class contains a method (`calculateSteps`) that calculates the number of motor steps necessary to move the cobras from their home positions to the target positions.

This class still needs some improvements, since it doesn't take yet into account the motor map errors and the movement from other positions that are not the home positions.

## CobraGroup.py

Defines the `CobraGroup` class. This class is used to represent the properties of a given group of cobras. These properties could be taken from a `CobrasCalibrationProduct` instance. Otherwise it will use default (perfect) properties: expected centers and link lengths, and constant motor maps.

The class contains a `MotorMapGroup` instance plus several methods to calculate the cobras patrol areas, the home positions, the fiber and elbow positions for a given configuration, and the theta and phi angles for each cobra.

## Bench.py

Defines the `Bench` class. This class is used to represent the PFI bench. It consists in a group of cobras plus some additional bench properties, like the bench center, the radius, and the cobra associations (array with the nearest cobras to a given cobra, which are those than can collide with each other).

This class has several methods that are useful to deal with cobra collisions. For example, `calculateCobraAssociationCollisions()` gives you the collisions between cobra association for a given fiber configuration. The fibers could be at the target positions (end positions), but they could also be positions along the cobra trajectories.

## TargetGroup.py

Defines the `TargetGroup` class. This class is used to represent the properties (xy coordinates and source id) of a group of PFS targets.

The class can deal with NULL targets, which are the targets of unassigned cobras.

Probably in the future we should add a target priority property...

## TargetSelector.py

Defines the `TargetSelector` class. This class is used to select (assign) targets to a given `Bench` instance. It takes as inputs a `TargetGroup` instance and a `Bench` instance and generates a new `TargetGroup` instance containing the targets assigned to each cobra in the the Bench. In some cases the these targets will be NULL targets.

This class is meant to be extended, by implementing the `run()` and `selectTargets()` methods. We have currently two example subclasses:
 * `DistanceTargetSelector.py`, which selects targets based on their distance to the cobra centers (similar to what the MATLAB code was doing)
 * `RandomTargetSelector.py`, which selects targets randomly (from the subset of targets that can be reach by each cobra).
 
In the future we can have a `NetflowSelector` subclass that selects targets based on the Netflow algorithm.

The `TargetSelector` class has a method to avoid end point collisions. This method can be run optionally, and if it's used it will reassign targets to cobras until the end point collisions are minimized. However, there could be cases when the collisions cannot be avoided because the two colliding cobras have only one possible target each. The MATLAB code was leaving one of the cobras unassigned, while in the python code we leave the two cobras assigned, which will generate an end collision at the end. We can change this when we decide what is the best thing to do (maybe look at a sky position?).  Also, the Netflow implementation might solve this problem automatically...

## TrajectoryGroup.py

Defines the `TrajectoryGroup` class. This class is used to represent the properties of a given group of cobra trajectories: final fiber positions, trajectory steps, theta and phi movement directions (positive or negative), theta and phi movement strategies (early or late), fiber and elbow positions at each step in the trajectory.

The class uses the motor maps inside the `Bench` instance to get the motor steps that are needed to reach the final target positions. It contains a method to calculate the cobra collisions along the trajectories: `calculateCobraAssociationCollisions()`.

## CollisionSimulator.py

Defines the `CollisionSimulator` class. This class is used to simulate a PFS observation for a given `Bench` instance and a `TargetGroup` instance (one target for each cobra, the result of running an specific `TargetSelector`).

This class contains all the logic designed to avoid trajectory collisions. When it's run, it calculates the cobra trajectories for the default theta and phi movement directions (positive or negative) and strategies (early or late) and detects trajectory collisions. If a collision is found, it changes the movement directions and strategies of the involved cobras until the collisions are minimized. Again, in some cases this cannot be avoided.

After running an instance of this class, one can access the finally adopted theta and phi movement directions and strategies, the cobras trajectories (an instance from the `TragetoryGroup` class) and all the unsolved end point and trajectory collisions.


# Utility modules

## AttributePrinter.py

Defines the `AttributePrinter` class. This class overrides the `__str__` method, printing all the instance attributes in their string representation. Any class that subclass it will inherit this property.

## plotUtils.py

This module contains several methods to create plots using matplotlib.

## targetUtils.py

This module contains methods to generate different target distributions.
