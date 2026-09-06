import numpy as np
import pytest

from tudatpy.plotting import covariance_ellipsoid


def test_covariance_ellipsoid_axis_lengths_and_center():
    covariance = np.diag([4.0, 1.0, 0.25])
    center = np.array([10.0, -2.0, 3.0])

    x, y, z = covariance_ellipsoid(covariance, center, sigma=2.0, resolution=101)

    assert x.shape == (101, 101)
    assert np.max(x) == pytest.approx(center[0] + 4.0, abs=2.0e-3)
    assert np.min(x) == pytest.approx(center[0] - 4.0, abs=2.0e-3)
    assert np.max(y) == pytest.approx(center[1] + 2.0, abs=2.0e-3)
    assert np.min(y) == pytest.approx(center[1] - 2.0, abs=2.0e-3)
    assert np.max(z) == pytest.approx(center[2] + 1.0)
    assert np.min(z) == pytest.approx(center[2] - 1.0)


@pytest.mark.parametrize(
    "covariance",
    (
        np.eye(2),
        np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        np.diag([1.0, 1.0, -1.0]),
    ),
)
def test_covariance_ellipsoid_rejects_invalid_covariance(covariance):
    with pytest.raises(ValueError):
        covariance_ellipsoid(covariance)
