"""

Collection of unit tests for the TargetSelector._outsideCircles method.

"""

import numpy as np

from ics.cobraOps.TargetSelector import TargetSelector


def bruteForce(points, centers, radii):
    return np.all(np.abs(points[:, np.newaxis] - centers) > radii, axis=1)


def test_matchesBruteForce():
    rng = np.random.default_rng(0)
    points = 50 * (rng.random(2000) + 1j * rng.random(2000))
    centers = 50 * (rng.random(300) + 1j * rng.random(300))
    radii = 1 + rng.random(300)
    assert np.array_equal(TargetSelector._outsideCircles(points, centers, radii),
                          bruteForce(points, centers, radii))


def test_pointsOnTheCircleEdge():
    centers = np.array([0j, 10 + 0j])
    radii = np.array([1.25, 1.2276])
    angles = np.exp(1j * np.linspace(0, 2 * np.pi, 7))
    points = (centers + radii * angles[:, np.newaxis]).ravel()
    assert np.array_equal(TargetSelector._outsideCircles(points, centers, radii),
                          bruteForce(points, centers, radii))


def test_scalarZeroAndNegativeRadius():
    points = np.array([0j, 1e-12 + 0j, 1 + 1j])
    centers = np.array([0j])
    for radius in (0, 0.5, -1):
        assert np.array_equal(
            TargetSelector._outsideCircles(points, centers, radius),
            bruteForce(points, centers, radius))


def test_emptyInputs():
    noPoints = TargetSelector._outsideCircles(
        np.zeros(0, dtype=complex), np.array([0j]), 1)
    assert noPoints.shape == (0,)
    noCircles = TargetSelector._outsideCircles(
        np.array([0j, 1j]), np.zeros(0, dtype=complex), 1)
    assert np.all(noCircles)
