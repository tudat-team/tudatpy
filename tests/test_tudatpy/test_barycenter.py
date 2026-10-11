import numpy as np
import pytest

from tudatpy.dynamics import environment, environment_setup


def test_barycenter_registration():
    settings = environment_setup.BodyListSettings("SSB", "ECLIPJ2000")
    states = {
        "Jupiter": np.array([7.8e11, 0, 0, 0, 13000, 0]),
        "Io": np.array([7.8e11 + 4.2e8, 0, 0, 0, 30000, 0]),
    }
    for name, state in states.items():
        settings.add_empty_settings(name)
        settings.get(name).ephemeris_settings = environment_setup.ephemeris.constant(state)
    bodies = environment_setup.create_system_of_bodies(settings)
    center = bodies.add_barycenter("Jupiter_barycenter", ["Jupiter", "Io"], [1000, 1])

    assert isinstance(center, environment.Barycenter)
    assert not isinstance(center, environment.Body)
    assert bodies.does_barycenter_exist(center.name)
    assert not bodies.does_body_exist(center.name)
    assert bodies.get_barycenter(center.name) is center
    assert bodies.list_of_barycenters == [center.name]
    np.testing.assert_allclose(center.weights, np.array([1000, 1]) / 1001)
    expected = (1000 * states["Jupiter"] + states["Io"]) / 1001
    np.testing.assert_allclose(center.ephemeris.cartesian_state(0), expected)
    for name in center.members:
        np.testing.assert_allclose(
            center.ephemeris.cartesian_state(0) + center.get_member_offset(name, 0),
            states[name],
        )
    with pytest.raises(RuntimeError):
        bodies.add_barycenter("overlap", ["Io"], [1])
    with pytest.raises(RuntimeError):
        bodies.get(center.name)
