import math
from datetime import datetime
from tudatpy.astro import time_conversion


def test_datetime_conversions():

    python_datetime = datetime.fromisoformat("2023-06-20T00:05:23.281765")
    tudat_datetime = time_conversion.datetime_to_tudat(python_datetime)
    tudat_datetime_string = tudat_datetime.iso_string(False, 12)
    python_datetime_reconstructed = time_conversion.datetime_to_python(tudat_datetime)

    while tudat_datetime_string[-1] == "0":
        tudat_datetime_string = tudat_datetime_string[:-1]

    assert tudat_datetime_string == str(python_datetime)

    assert str(python_datetime_reconstructed) == str(python_datetime)

    assert "2023-06-20 00:05:23.281765" == str(python_datetime)

    julian_day = 2443494.5
    tudat_datettime = time_conversion.datetime_to_tudat(
        time_conversion.julian_day_to_calendar_date(julian_day)
    )
    assert julian_day == tudat_datettime.julian_day()

    # pre_1970_date = datetime(1960,1,1,13,7,23,891234)
    # julian_day_to_check = time_conversion.calendar_date_to_julian_day(pre_1970_date)
    #
    # assert math.floor( julian_day_to_check ) == 2436935
    # seconds_in_day_to_check = ( julian_day_to_check - math.floor( julian_day_to_check ) ) * 86400.0
    # manual_seconds_in_day = 1 * 3600 + 7 * 60 + 23.891234
    # assert math.fabs( manual_seconds_in_day - seconds_in_day_to_check ) < 1.0E-4
    # print('DATE TIME TEST SUCCESFUL')
