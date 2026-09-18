import pandas as pd
import pytest

from tudatpy.data_input.environment_data.sbdb import SBDBbatch


@pytest.mark.remote_data(service="JPL SBDB")
def test_sbdb_batch_downloads_pages_and_filters():
    """Check that a live SBDB catalogue query can select Vesta."""
    batch = SBDBbatch(fields=["pdes", "name"])

    assert batch.get(4).iloc[0]["name"] == "Vesta"


def test_sbdb_batch_csv_round_trip(tmp_path):
    """Check that CSV loading preserves leading zeros in designations."""
    path = tmp_path / "sbdb.csv"
    pd.DataFrame({"pdes": ["00001"], "name": ["Ceres"]}).to_csv(path, index=False)

    batch = SBDBbatch(path)

    assert batch.get("00001").iloc[0]["name"] == "Ceres"
