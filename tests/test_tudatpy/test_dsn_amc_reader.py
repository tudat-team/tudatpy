import pytest

from tudatpy.data_input.environment_data import read_dsn_amc_file


def test_native_amc_reader_exposure(tmp_path):
    path = tmp_path / "amc.txt"
    path.write_text(
        "# Path delay data file\n"
        "P 1495173540 0024.985 0132.403 0004.261 0206.209 -00002 "
        "010 0 0 0 3 0 0 009.6 023\n"
    )
    data = read_dsn_amc_file(str(path))
    assert data.file_name == str(path)
    assert len(data.records) == 1
    record = data.records[0]
    assert record.utc_time == 548445540.0
    assert record.zenith_wet_delay == pytest.approx(0.04261)
    assert record.zenith_dry_delay == pytest.approx(2.06209)
    assert record.algorithm_id == "010"
    assert record.surface_flag == 3
    assert record.cloud_liquid == -2.0
    with pytest.raises(AttributeError):
        record.surface_flag = 0


def test_native_amc_reader_reports_bad_line(tmp_path):
    path = tmp_path / "malformed.txt"
    path.write_text("# Path delay data file\nP incomplete\n")
    with pytest.raises(RuntimeError, match="line 2"):
        read_dsn_amc_file(str(path))
