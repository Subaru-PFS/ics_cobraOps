"""

Defines some constants related with the cobras and the simulation code.

Consult the following papers for more detailed information:

  http://adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  http://adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  http://adsabs.harvard.edu/abs/2016arXiv160801075T

"""

MODULE_FIRST_LINE_LENGTH = 29
"""The number of cobras in the first line of a module."""

MODULE_SECOND_LINE_LENGTH = 28
"""The number of cobras in the second line of a module."""

MODULES_PER_SECTOR = 14
"""The number of modules in one PFI sector."""
 
COBRAS_SEPARATION = 8.0
"""The default separation between two consecutive cobras in mm."""

COBRA_LINK_LENGTH = 2.375
"""The default cobra link length in mm."""

COBRA_LINK_RADIUS = 1.0
"""The default cobra link radius in mm."""

PHI_SAFETY_ANGLE = 0.1
"""Phi safety angle value in radians (~ 5.7 deg)."""

HOMES_THETA_DISTANCE = 6.72
"""The default homes theta distance in radians (~385 deg)."""

MOTOR_MAP_ANGULAR_STEP = 0.06283
"""The default angular step in radians used to measure the cobra motor maps
(~3.6 deg)."""

MOTOR_STEPS_PER_ANGULAR_STEP = 50.0
"""The default number of motor steps needed to cover a motor angular step."""

MOTOR_NOISE_ALPHA = 0.07
"""The default cobras motor noise normalization factor."""

MOTOR_NOISE_BETA = 0.50
"""The default cobras motor noise exponent."""

NULL_TARGET_INDEX = -1
"""Integer value used to indicate the index of a null target."""

NULL_TARGET_POSITION = 0j
"""Complex value used to indicate the position of a null target."""

NULL_TARGET_ID = "NULL"
"""Unicode value used to indicate the id of a null target."""

TRAJECTORY_STEPS = 110
"""The total number of time steps in a cobra trajectory."""

TRAJECTORY_STEP_WIDTH = 50
"""The trajectory step width in units of motor steps."""
