"""Load both EOP formats through Python with their numerical and optional fields intact."""

from pathlib import Path

import numpy as np

from tudatpy.data_input import resource_paths
from tudatpy.data_input.environment_data import earth_orientation

FIXTURES = (
    Path(resource_paths.get_test_data_path()) / "pr1092" / "src/io/earth_orientation_fixtures"
)


def test_earth_orientation_formats_preserve_units_and_missing_values():
    """Read C04 values in SI units and preserve missing prediction fields from finals2000A."""
    data = earth_orientation.read_earth_orientation_data(str(FIXTURES / "modern_c04.txt"))
    record = next(iter(data.records.values()))

    # The published arcsecond polar motion is exposed as radians alongside the UTC offset in seconds.
    np.testing.assert_allclose(
        record.polar_motion[0], -0.134064 * np.pi / 648000, rtol=1e-15, atol=0.0
    )
    assert record.ut1_minus_utc == 0.3968898

    finals = earth_orientation.read_earth_orientation_data(
        str(FIXTURES / "finals2000a.txt"), earth_orientation.EarthOrientationFormat.finals2000a
    )
    predicted = next(reversed(finals.records.values()))

    # Absent predicted LOD and pole offsets remain None rather than supplied zeros.
    assert predicted.length_of_day is None
    assert predicted.celestial_pole_offsets is None
