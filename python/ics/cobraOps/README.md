# Main modules

## cobraConstants.py

Defines some constants related with the cobras and the simulation code, like the default cobra link lengths, the number of cobras per module, or the default motor maps angular step size.

## CobraGroup.py

Defines the `CobraGroup` class. This class is used to represent the properties of a given group of cobras. These properties are taken from a Cobras calibration product.

The class contains several methods to calculate the cobras patrol areas, the home positions, the fiber and elbow positions for a given configuration, and the theta and phi angles for each cobra.

## Bench.py

Defines the `Bench` class. This class is used to represent the PFI bench. It consists in a group of cobras plus some additional bench properties, like the bench center, the radius, and the cobra associations (array with the nearest cobras to a given cobra, which are those than can collide with each other).

This class has several methods that are useful to deal with cobra collisions.

## TargetGroup.py

Defines the `TargetGroup` class. This class is used to represent the properties (xy coordinates, source id and priority) of a group of PFS targets.

The class can deal with NULL targets, which are the targets of unassigned cobras.

## TargetSelector.py

Defines the `TargetSelector` abstract class. This class is used to select (assign) targets to a given `Bench` instance. It takes as inputs a `TargetGroup` instance and a `Bench` instance and generates a new `TargetGroup` instance containing the targets assigned to each cobra in the the Bench. In some cases these targets will be NULL targets.

This class is meant to be extended, by implementing the `run()` and `selectTargets()` methods. We have currently three example subclasses:
 * `DistanceTargetSelector.py`, which selects targets based on their distance to the cobra centers (similar to what the MATLAB code was doing).
 * `PriorityTargetSelector.py`, which selects targets based on their scientific priority.
 * `RandomTargetSelector.py`, which selects targets randomly (from the subset of targets that can be reached by each cobra).

In the future we can have a `NetflowTargetSelector` subclass that selects targets based on the Netflow algorithm.

## CollisionSimulator2.py

Defines the `CollisionSimulator2` class. This class is used to simulate a PFS observation for a given `Bench` instance and a `TargetGroup` instance (one target for each cobra, the result of running an specific `TargetSelector`).


# Utility modules

## AttributePrinter.py

Defines the `AttributePrinter` class. This class overrides the `__str__` method to print all instance attributes followed by their respective values:

```
attribute_1 = value1
attribute_2 = value2
...
attribute_n = value_n
```

Any class that subclass it will inherit this method.

## plotUtils.py

This module contains several methods to create plots using matplotlib.

## targetUtils.py

This module contains methods to generate different target distributions.
