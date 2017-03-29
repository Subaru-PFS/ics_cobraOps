# PFS survey target test data set 

## Introduction
More realistic target list is crucial to develop every software packages in the PFS system, especially Exposure Targeting Software (ETS) and cobra collision studies. This package includes preliminary target lists of each scientific component (cosmology, galaxy evolution, and galactic archaeology) in the coming PFS SSP survey. Most of targets in the catalog are taken from real astronomical data including HSC imaging data but will be updated after reflecting results of the target selection in each science working group (WG). **Note that DO NOT distribute the target lists to others outside PFS collaboration and DO NOT use them for other studies than PFS related things, because some unpublished data are included in the lists.**

## Description

### Target lists
* data/pfs_preliminary_target_cosmology.dat
* data/pfs_preliminary_target_galaxy.dat
* data/pfs_preliminary_target_archaeology.dat

These are a target list of each science WG (cosmology, galaxy evolution, and galactic archaeology). Each file contains the following information:

~~~~~
#  1  ID
#  2  R.A.          [deg.]
#  3  Dec.          [deg.]
#  4  Exposure Time [sec.]
#  5  Priority      [1(highest)-15(lowest)]
#  6  Magnitude     [AB mag]
#  7  Redshift
#  8  Object Type
~~~~~

This file can be a direct input to the ETS package.

### ETS and its outputs

* data/pfs_ets_output_cosmology.dat
* data/pfs_ets_output_galaxy.dat
* data/pfs_ets_output_archaeology.dat
* docs/coordinate_notes.pdf

These are outputs of the current version ETS. You can get the git repository on:
https://github.com/Subaru-PFS/ets_fiber_assigner

Please clone or download and see README.md for the package itself and its usage. **Please compile the source code first by typing "make" if you use the code for your own studies.** The output file can be taken by running the code as follows:

~~~~~
ets_demo assigner=naive input=data/pfs_preliminary_target_cosmology.dat n_exposures=10 output=data/pfs_ets_output_cosmology.dat time=2017-10-15T10:00:00Z
ets_demo assigner=naive input=data/pfs_preliminary_target_galaxy.dat n_exposures=10 output=data/pfs_ets_output_galaxy.dat time=2017-02-15T10:00:00Z
ets_demo assigner=naive input=data/pfs_preliminary_target_archaeology.dat n_exposures=10 output=data/pfs_ets_output_archaeology.dat time=2017-05-15T10:00:00Z
~~~~~

The output file contains the following information:

~~~~~
#  1  Target ID
#  2  Fiber  ID
#  3  X on PFI
#  4  Y on PFI
#  5  R.A.     [deg.]
#  6  Dec.     [deg.]
~~~~~

The results of each configuration is separated by two lines to the next set of configuration like this:

~~~~~
Exposure 0: duration 900s, AZ: 146.098, ALT 61.3219 PA: -3.2
  Target     Fiber          X          Y         RA        DEC
C017618         1   -8.32059    1.29269   33.97465   -4.50895
C017599         2  -14.20081   -4.12232   33.95721   -4.52689
...
C011682      2393  196.31990 -115.09705   34.62019   -4.82737
C011709      2394  191.18638 -106.76909   34.60433   -4.80356
Exposure 1: duration 900s, AZ: 146.108, ALT 61.3369 PA: -3.2
  Target     Fiber          X          Y         RA        DEC
C017565         1   -8.63457    1.06628   33.96151   -4.49960
C017613         2  -13.87411   -2.80201   33.94580   -4.51260
...
~~~~~

**If you are interested in only the first fiber configuration, just extract the first block from the output file.** Please refer to a note (docs/coordinates_notes.pdf) by Yuki Moritani (Kavli IPMU) for details on the coordinate system. The coordinate transformation is under discussion and will probably be updated in the future.


## Contact
Kiyoto Yabe (Kavli IPMU)
e-mail: kiyoto.yabe@ipmu.jp


