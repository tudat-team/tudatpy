import numpy as np
import pytest
from astropy import units as u
from astropy.time import Time

from tudatpy.data.sbdb import SBDBquery


class FakeSBDBValue:
    # Small fake object that behaves like an SBDB value with a .value attribute
    def __init__(self, value):
        self.value = value


@pytest.fixture
def encke():
    # Create a fake SBDB query for comet 2P/Encke
    result = SBDBquery.__new__(SBDBquery)

    result.MPCcode = "2P"
    result.query = {
        "object": {
            "fullname": "2P/Encke",
        },
        "orbit": {
            "model_pars": {
                "A1": FakeSBDBValue(1.0e-8),
                "A2": FakeSBDBValue(-2.0e-9),
                "A3": FakeSBDBValue(3.0e-10),
                "DT": FakeSBDBValue(1.25),
            },
            "first_obs": "1986-01-17",
            "last_obs": "2024-01-01",
            "elements": {
                "q": FakeSBDBValue(0.339),
                "tp": FakeSBDBValue(2460000.5),
            },
        },
    }

    return result


def test_nongrav_params(encke):
    expected = np.array(
        [
            (1.0e-8 * u.au / u.day**2).to(u.m / u.s**2).value,
            (-2.0e-9 * u.au / u.day**2).to(u.m / u.s**2).value,
            (3.0e-10 * u.au / u.day**2).to(u.m / u.s**2).value,
        ]
    )
    assert encke.nongrav_params == pytest.approx(expected)


def test_nongrav_params_missing(encke):
    del encke.query["orbit"]["model_pars"]["A2"]
    del encke.query["orbit"]["model_pars"]["A3"]
    result = encke.nongrav_params
    assert result[1] == 0
    assert result[2] == 0


def test_Dt(encke):
    assert encke.Dt == pytest.approx((1.25 * u.day).to(u.s).value)


def test_Dt_missing(encke):
    del encke.query["orbit"]["model_pars"]["DT"]
    with pytest.raises(ValueError, match="Asymmetry parameter DT is not available for object"):
        encke.Dt


def test_first_obs(encke):
    expected = Time("1986-01-17").jd
    assert encke.first_obs == pytest.approx(expected)


def test_last_obs(encke):
    expected = Time("2024-01-01").jd
    assert encke.last_obs == pytest.approx(expected)


def test_first_obs_missing(encke):
    del encke.query["orbit"]["first_obs"]
    with pytest.raises(ValueError, match="first observation is not available"):
        encke.first_obs


def test_perihelion(encke):
    expected = (0.339 * u.au).to(u.m).value
    assert encke.perihelion == pytest.approx(expected)


def test_perihelion_missing(encke):
    del encke.query["orbit"]["elements"]["q"]
    with pytest.raises(ValueError, match="Perihelion distance is not available for object"):
        encke.perihelion


def test_time_perihelion(encke):
    assert encke.time_perihelion == 2460000.5


def test_time_perihelion_missing(encke):
    del encke.query["orbit"]["elements"]["tp"]
    with pytest.raises(ValueError, match="Perihelion time is not available"):
        encke.time_perihelion
