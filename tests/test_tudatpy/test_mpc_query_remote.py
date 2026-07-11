import pytest
from tudatpy.data_access.tracking.mpc import BatchMPC
import requests
from astroquery.mpc import MPC


def _mpc_server_reachable(timeout: float = 10.0) -> tuple[bool, str]:
    """Lightweight reachability probe for the MPC web service.

    Only needed to check that the internet works.
    For this reaosn, it hits the small, static observatory-codes page rather than the
    observations DB endpoint, so the check stays fast and doesn't
    exercise the query logic that the real test below is exercising.
    """
    try:
        response = requests.head(MPC.OBSERVATORY_CODES_URL, timeout=timeout, allow_redirects=True)
        response.raise_for_status()
    except requests.exceptions.RequestException as exc:
        return False, str(exc)
    return True, ""


@pytest.fixture(scope="session")
def require_mpc_server():
    """Skip (not fail) remote MPC tests if the service or network is down."""
    ok, reason = _mpc_server_reachable()
    if not ok:
        pytest.skip(f"MPC server unreachable, skipping remote test: {reason}")


@pytest.mark.remote_data
class TestBatchMPCGetObservationsRemote:
    def test_get_observations_returns_data_for_known_object(self, require_mpc_server):
        batch = BatchMPC()
        batch.get_observations([1])  # 1 Ceres has a long, stable observation history

        assert batch.size > 0
        assert {"number", "epoch", "RA", "DEC", "observatory"}.issubset(batch.table.columns)
        assert (
            batch.epoch_start < 0
        )  # first Ceres observation was recorded before Jan, 01, 01, 2000
        assert batch.epoch_end >= batch.epoch_start

    def test_get_observations_invalid_code_raises_value_error(self, require_mpc_server):
        batch = BatchMPC()
        with pytest.raises(ValueError):
            batch.get_observations(["not_a_real_object_designation_zzz"])
