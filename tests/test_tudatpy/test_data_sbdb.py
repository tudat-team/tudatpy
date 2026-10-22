import numpy as np
import pytest
from astropy import units as u
from astropy.time import Time

from tudatpy.data.sbdb import SBDBquery
from tudatpy.astro.time_representation import DateTime


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
                "A1": 1.0e-8 * u.au / u.day**2,
                "A2": -2.0e-9 * u.au / u.day**2,
                "A3": 3.0e-10 * u.au / u.day**2,
                "DT": 1.25 * u.day,
            },
            "first_obs": "1986-01-17",
            "last_obs": "2024-01-01",
            "elements": {
                "q": 0.339 * u.au,
                "tp": 2460000.5 * u.day,
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
    encke.query["orbit"]["model_pars"].pop("A2")
    encke.query["orbit"]["model_pars"].pop("A3")
    result = encke.nongrav_params
    assert result[1] == 0
    assert result[2] == 0


def test_Dt(encke):
    assert encke.Dt == pytest.approx((1.25 * u.day).to(u.s).value)


def test_Dt_missing(encke):
    encke.query["orbit"]["model_pars"].pop("DT")
    with pytest.raises(ValueError, match="Asymmetry parameter DT is not available for object"):
        encke.Dt


def test_first_obs(encke):
    expected = DateTime.from_iso_string("1986-01-17T00:00:00.00").to_epoch()
    assert encke.first_obs == pytest.approx(expected)


def test_last_obs(encke):
    expected = DateTime.from_iso_string("2024-01-01T00:00:00.00").to_epoch()
    assert encke.last_obs == pytest.approx(expected)


def test_first_obs_missing(encke):
    encke.query["orbit"].pop("first_obs")
    with pytest.raises(ValueError, match="first observation is not available"):
        encke.first_obs


def test_perihelion(encke):
    expected = (0.339 * u.au).to(u.m).value
    assert encke.perihelion == pytest.approx(expected)


def test_perihelion_missing(encke):
    encke.query["orbit"]["elements"].pop("q")
    with pytest.raises(ValueError, match="Perihelion distance is not available for object"):
        encke.perihelion


def test_time_perihelion(encke):
    assert encke.time_perihelion == DateTime.from_julian_day(2460000.5).to_epoch()


def test_time_perihelion_missing(encke):
    encke.query["orbit"]["elements"].pop("tp")
    with pytest.raises(ValueError, match="Perihelion time is not available"):
        encke.time_perihelion
