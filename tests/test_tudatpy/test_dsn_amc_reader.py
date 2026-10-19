import pytest

from tudatpy.data_input.environment_data import get_valid_dsn_amc_zenith_delays, read_dsn_amc_file


def test_dsn_amc_reader(tmp_path):
    # Juno PJ6 records (DSS-25): a nominal record, and a record with non-zero surface flag
    path = tmp_path / "amc.txt"
    path.write_text(
        "# Path delay data file\n"
        "P 1495173540 0024.985 0132.403 0004.261 0206.209 -00002 010 0 0 0 0 0 0 009.6 023\n"
        "P 1495173570 0024.985 0132.403 0004.263 0206.209 -00002 010 0 0 0 3 0 0 009.6 023\n"
    )
    records = read_dsn_amc_file(str(path))
    assert len(records) == 2
    assert records[0].utc_time == 548445540.0
    assert records[0].zenith_wet_delay == pytest.approx(0.04261)
    assert records[0].zenith_dry_delay == pytest.approx(2.06209)
    assert records[0].algorithm_id == "010"
    assert records[1].surface_flag == 3

    zenith_delays = get_valid_dsn_amc_zenith_delays(records)
    assert list(zenith_delays) == [548445540.0]
    assert zenith_delays[548445540.0] == pytest.approx([2.06209, 0.04261])
