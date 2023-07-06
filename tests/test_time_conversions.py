# import sys
# sys.path.insert(0, "/home/dominic/Tudat/tudat-bundle/tudat-bundle/cmake-build-default/tudatpy")

from datetime import  datetime
from tudatpy.kernel.astro import time_conversion


def test_datetime_conversions():

    python_datetime = datetime.fromisoformat('2023-06-20T00:05:23.281765')
    tudat_datetime = time_conversion.datetime_to_tudat( python_datetime )
    tudat_datetime_string = tudat_datetime.iso_string( False )
    python_datetime_reconstructed = time_conversion.datetime_to_python( tudat_datetime )

    while (tudat_datetime_string[-1] == '0'):
        tudat_datetime_string = tudat_datetime_string[:-1]

    assert tudat_datetime_string == str( python_datetime )

    assert str( python_datetime_reconstructed ) == str( python_datetime )

    assert '2023-06-20 00:05:23.281765' == str(python_datetime)

    julian_day = 2443494.5
    tudat_datettime = time_conversion.datetime_to_tudat(time_conversion.julian_day_to_calendar_date(julian_day))
    assert julian_day == tudat_datettime.julian_day( )


