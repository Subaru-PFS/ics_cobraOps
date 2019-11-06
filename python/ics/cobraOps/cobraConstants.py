"""

Defines some constants related with the cobras and the simulation code.

Consult the following papers for more detailed information:

  https://ui.adsabs.harvard.edu/abs/2012SPIE.8450E..17F
  https://ui.adsabs.harvard.edu/abs/2014SPIE.9151E..1YF
  https://ui.adsabs.harvard.edu/abs/2016arXiv160801075T
  https://ui.adsabs.harvard.edu/abs/2018SPIE10707E..28Y
  https://ui.adsabs.harvard.edu/abs/2018SPIE10702E..1CT

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

MOTOR1_STEP_SIZE = 0.06
"""The default stage 1 motor step size in degrees."""

MOTOR2_STEP_SIZE = 0.12
"""The default stage 2 motor step size in degrees."""

MOTOR_NOISE_ALPHA = 0.07
"""The default cobras motor noise normalization factor."""

MOTOR_NOISE_BETA = 0.50
"""The default cobras motor noise exponent."""

BLACK_DOT_RELATIVE_POSITION = 0 + 2.35j 
"""The default cobra back dot position relative to the cobra center in mm."""

BLACK_DOT_RADIUS = 1.375
"""The default cobra back dot radius in mm."""

NULL_TARGET_INDEX = -1
"""Integer value used to indicate the index of a null target."""

NULL_TARGET_POSITION = 0j
"""Complex value used to indicate the position of a null target."""

NULL_TARGET_ID = "NULL"
"""Unicode value used to indicate the id of a null target."""

NULL_TARGET_PRIORITY = -1
"""Float value used to indicate the priority of a null target."""
