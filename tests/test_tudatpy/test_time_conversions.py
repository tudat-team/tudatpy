from datetime import datetime
import pytest
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime

# ---------------------------------------------------------------------------
# Test cases for Python <-> DateTime conversion
# ---------------------------------------------------------------------------
# Covers: J2000, Unix epoch, pre-epoch, leap days, DST spring-forward gap,
# DST fall-back ambiguity, microsecond precision, year/month boundaries.
PYTHON_DATETIME_CASES = [
    datetime(2000, 1, 1, 12, 0, 0),  # J2000 epoch
    datetime(2000, 1, 1, 0, 0, 0),  # Y2K midnight
    datetime(1970, 1, 1, 0, 0, 0),  # Unix epoch
    datetime(1969, 12, 31, 23, 59, 59),  # one second before Unix epoch
    datetime(1960, 1, 1, 13, 7, 23, 891234),  # pre-1970 with microseconds
    datetime(1901, 5, 17, 8, 34, 30),  # early 20th century
    datetime(2000, 2, 29, 8, 34, 30, 234567),  # leap day Y2K
    datetime(2004, 2, 29, 18, 34, 30, 234567),  # leap day regular
    datetime(1900, 1, 1, 0, 0, 0),  # year 1900 (not a leap year)
    datetime(2023, 6, 20, 0, 5, 23, 281765),  # regular date with microseconds
    datetime(2023, 12, 31, 23, 59, 59, 999999),  # last microsecond of year
    datetime(2000, 1, 1, 0, 0, 0, 1),  # single microsecond
    datetime(2000, 12, 17, 12, 0, 0),  # near solstice, noon
    # DST spring-forward (Amsterdam CET→CEST): 2:00 → 3:00, so 2:xx doesn't exist locally
    datetime(2012, 3, 25, 1, 59, 22),  # last moment in CET
    datetime(2012, 3, 25, 2, 59, 22),  # non-existent in local time (DST gap)
    datetime(2012, 3, 25, 3, 59, 22),  # first moment in CEST
    # DST fall-back (Amsterdam CEST→CET): 3:00 → 2:00, so 2:xx occurs twice locally
    datetime(2012, 10, 28, 2, 30, 0),  # ambiguous local time
]

TUDAT_DATETIME_MICROSECONDS_OVERFLOW = [
    (
        DateTime(
            2024, 12, 31, 23, 59, 59.9999996
        ),  # probe (with microsecond overflow, last day of year)
        datetime(2025, 1, 1, 0, 0, 0, 0),  # expected result
    ),
    (
        DateTime(
            2024, 12, 30, 22, 51, 59.9999996
        ),  # probe (with microsecond overflow, random day of year)
        datetime(2024, 12, 30, 22, 52, 0, 0),
    ),
]


@pytest.mark.parametrize("dt", PYTHON_DATETIME_CASES)
def test_from_python_datetime_preserves_fields(dt):
    """from_python_datetime must read calendar fields literally, independent of local timezone."""
    tudat_dt = DateTime.from_python_datetime(dt)

    assert tudat_dt.year == dt.year
    assert tudat_dt.month == dt.month
    assert tudat_dt.day == dt.day
    assert tudat_dt.hour == dt.hour
    assert tudat_dt.minute == dt.minute
    expected_seconds = dt.second + dt.microsecond / 1e6
    assert float(tudat_dt.seconds) == pytest.approx(expected_seconds, abs=1e-9)


@pytest.mark.parametrize("dt", PYTHON_DATETIME_CASES)
def test_to_python_datetime_preserves_fields(dt):
    """to_python_datetime must return a naive datetime with the same calendar fields."""
    tudat_dt = DateTime(
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second + dt.microsecond / 1e6
    )
    result = tudat_dt.to_python_datetime()

    assert result.year == dt.year
    assert result.month == dt.month
    assert result.day == dt.day
    assert result.hour == dt.hour
    assert result.minute == dt.minute
    assert result.second == dt.second
    assert result.microsecond == pytest.approx(dt.microsecond, abs=1)


@pytest.mark.parametrize("dt", PYTHON_DATETIME_CASES)
def test_python_datetime_round_trip(dt):
    """Python datetime -> DateTime -> Python datetime must recover the original fields."""
    result = DateTime.from_python_datetime(dt).to_python_datetime()

    assert result.year == dt.year
    assert result.month == dt.month
    assert result.day == dt.day
    assert result.hour == dt.hour
    assert result.minute == dt.minute
    assert result.second == dt.second
    assert result.microsecond == pytest.approx(dt.microsecond, abs=1)


def test_datetime_conversions():

    python_datetime = datetime.fromisoformat("2023-06-20T00:05:23.281765")
    tudat_datetime = time_representation.DateTime.from_python_datetime(python_datetime)
    tudat_datetime_string = tudat_datetime.iso_string(False, 12)
    python_datetime_reconstructed = tudat_datetime.to_python_datetime()

    while tudat_datetime_string[-1] == "0":
        tudat_datetime_string = tudat_datetime_string[:-1]

    assert tudat_datetime_string == str(python_datetime)

    assert str(python_datetime_reconstructed) == str(python_datetime)

    assert "2023-06-20 00:05:23.281765" == str(python_datetime)

    julian_day = 2443494.5
    tudat_datettime = time_representation.DateTime.from_python_datetime(
        time_representation.DateTime.from_julian_day(julian_day).to_python_datetime()
    )
    assert julian_day == pytest.approx(tudat_datettime.to_julian_day(), abs=1e-9)


@pytest.mark.parametrize(
    ("dt", "expected"),
    TUDAT_DATETIME_MICROSECONDS_OVERFLOW,
)
def test_to_python_datetime_microseconds_overflow(dt, expected):
    """Verify that to_python_datetime correctly handles microsecond overflow."""
    result = dt.to_python_datetime()
    assert result == expected
