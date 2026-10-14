"""Binding smoke only; EOP parsing and scientific checks live in C++ tests."""

from pathlib import Path

from tudatpy.data_input.environment_data import earth_orientation

FIXTURES = Path(__file__).parents[1] / "test_tudat/src/io/earth_orientation_fixtures"


def test_raw_earth_orientation_namespace_and_array_binding():
    data = earth_orientation.read_earth_orientation_data(str(FIXTURES / "modern_c04.txt"))
    record = next(iter(data.records.values()))
    assert isinstance(data, earth_orientation.EarthOrientationData)
    assert isinstance(record, earth_orientation.EarthOrientationRecord)
    assert isinstance(record.polar_motion, list)
    assert len(record.polar_motion) == 2
    assert isinstance(record.polar_motion_rates, list)


def test_optional_earth_orientation_values_bind_to_none():
    data = earth_orientation.read_earth_orientation_data(
        str(FIXTURES / "finals2000a.txt"), earth_orientation.EarthOrientationFormat.finals2000a
    )
    predicted = next(reversed(data.records.values()))
    assert predicted.length_of_day is None
    assert predicted.celestial_pole_offsets is None
    assert predicted.celestial_pole_offset_uncertainty is None
