"""

Collection of unit tests for the TargetGroup class.

"""

import pytest
import numpy as np
import matplotlib.pyplot as plt

from ics.cobraOps.TargetGroup import TargetGroup


class TestTargetGroup():
    """A collection of tests for the TargetGroup class.

    """

    def test_constructor(self, finalPositions, targetIds, targetPriorities):
        # Create a new TargetGroup instance using only the target positions
        targets = TargetGroup(finalPositions)

        # Check that the number of targets and the positions are correct
        assert targets.nTargets == len(finalPositions)
        assert np.all(targets.positions == finalPositions)

        # Check that the default ids and priorities are correct
        assert np.all(targets.ids == ["0", "1", "2", "3"])
        assert np.all(targets.priorities == [1, 1, 1, 1])

        # Create a new TargetGroup instance with all the inputs
        targets = TargetGroup(finalPositions, targetIds, targetPriorities)

        # Check that the ids, positions and priorities are correct
        assert np.all(targets.ids == targetIds)
        assert np.all(targets.positions == finalPositions)
        assert np.all(targets.priorities == targetPriorities)

        # Create a new TargetGroup instance with one NULL target
        targetIds[1] = TargetGroup.NULL_TARGET_ID
        targets = TargetGroup(finalPositions, targetIds, targetPriorities)

        # Check that the ids, positions and priorities are correct
        notNull = targetIds != TargetGroup.NULL_TARGET_ID
        assert np.all(targets.ids[notNull] == targetIds[notNull])
        assert np.all(targets.positions[notNull] == finalPositions[notNull])
        assert np.all(targets.priorities[notNull] == targetPriorities[notNull])
        assert np.all(targets.ids[~notNull] == TargetGroup.NULL_TARGET_ID)
        assert np.all(
            targets.positions[~notNull] == TargetGroup.NULL_TARGET_POSITION)
        assert np.all(
            targets.priorities[~notNull] == TargetGroup.NULL_TARGET_PRIORITY)

    def test_constructor_exception(self, finalPositions, targetIds,
                                   targetPriorities):
        # Check that the constructor raises a exception when the input array
        # lengths do not coincide
        with pytest.raises(ValueError):
            TargetGroup(finalPositions, targetIds[1:], targetPriorities)

        with pytest.raises(ValueError):
            TargetGroup(finalPositions, targetIds, targetPriorities[2:])

    def test_fromFile_method(self, finalPositions, targetIds,
                            targetPriorities, tmpdir):
        # Write the targets information into a file
        fileName = str(tmpdir.join("myTestTargets.dat"))

        with open(fileName, "w+") as f:
            for p, i, pr in zip(finalPositions, targetIds, targetPriorities):
                f.write("%s, %s, %s, %s\n" % (np.real(p), np.imag(p), i, pr))

        # Create a new TargetGroup instance from the file
        targets = TargetGroup.fromFile(fileName)

        # Check that the ids, positions and priorities are correct
        assert np.all(targets.ids == targetIds)
        assert np.all(targets.positions == finalPositions)
        assert np.all(targets.priorities == targetPriorities)

    def test_saveToFile_method(self, finalPositions, targetIds,
                            targetPriorities, tmpdir):
        # Create a new TargetGroup instance with all the inputs
        targets = TargetGroup(finalPositions, targetIds, targetPriorities)

        # Write the targets information into a file
        fileName = str(tmpdir.join("myTestTargets.dat"))
        targets.saveToFile(fileName)

        # Read back the targets file
        targets = TargetGroup.fromFile(fileName)

        # Check that the ids, positions and priorities are correct
        assert np.all(targets.ids == targetIds)
        assert np.all(targets.positions == finalPositions)
        assert np.all(targets.priorities == targetPriorities)

    def test_select_method(self, finalPositions, targetIds, targetPriorities):
        # Create a new TargetGroup instance with all the inputs
        targets = TargetGroup(finalPositions, targetIds, targetPriorities)

        # Select some of the targets
        indices = [1, 2]
        selectedTargets = targets.select(indices)

        # Check that the result makes sense
        assert selectedTargets.nTargets == len(indices)
        assert np.all(selectedTargets.ids == targetIds[indices])
        assert np.all(selectedTargets.positions == finalPositions[indices])
        assert np.all(selectedTargets.priorities == targetPriorities[indices])

    def test_addToFigure_method(self, finalPositions, targetIds,
                                targetPriorities):
        # Create a new TargetGroup instance with all the inputs
        targets = TargetGroup(finalPositions, targetIds, targetPriorities)

        # Plot the targets
        plt.figure()
        targets.addToFigure()
        plt.show(block=False)
