import datetime
import numpy
import typing
__all__ = ['DateTime', 'TAI_to_TT', 'TCB_to_TDB', 'TCG_to_TT', 'TDB_to_TCB', 'TT_to_TAI', 'TT_to_TCG', 'TT_to_TDB_approximate', 'TimeScaleConverter', 'TimeScales', 'add_days_to_datetime', 'add_seconds_to_datetime', 'calculate_seconds_in_current_julian_day', 'calendar_date_to_days_since_epoch', 'calendar_date_to_julian_day', 'calendar_date_to_julian_day_since_epoch', 'date_time_from_epoch', 'date_time_from_iso_string', 'datetime_to_python', 'datetime_to_tudat', 'default_time_scale_converter', 'epoch_from_date_time_components', 'epoch_from_date_time_iso_string', 'get_days_in_month', 'is_leap_year', 'julian_day_to_calendar_date', 'julian_day_to_modified_julian_day', 'julian_day_to_seconds_since_epoch', 'modified_julian_day_to_julian_day', 'seconds_since_epoch_to_julian_centuries_since_epoch', 'seconds_since_epoch_to_julian_day', 'seconds_since_epoch_to_julian_years_since_epoch', 'tai_scale', 'tdb_scale', 'tt_scale', 'ut1_scale', 'utc_scale']

class DateTime:
    """
		"""
    day: int
    hour: int
    minute: int
    month: int
    seconds: float
    year: int

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, year: int, month: int, day: int, hour: int=12, minute: int=0, seconds: float=0.0) -> None:
        ...

    def day_of_year(self) -> int:
        ...

    def epoch(self) -> float:
        ...

    def iso_string(self, add_T: bool=False, number_of_digits_seconds: int=15) -> str:
        ...

    def julian_day(self) -> float:
        ...

    def modified_julian_day(self) -> float:
        ...

class TimeScaleConverter:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def convert_time(self, input_scale: TimeScales, output_scale: TimeScales, input_value: float, earth_fixed_position: numpy.ndarray=...) -> float:
        ...

    def get_time_difference(self, input_scale: TimeScales, output_scale: TimeScales, input_value: float, earth_fixed_position: numpy.ndarray=...) -> float:
        ...

class TimeScales:
    """Members:
	
	tai_scale
	
	tt_scale
	
	tdb_scale
	
	utc_scale
	
	ut1_scale
	"""
    __members__: typing.ClassVar[dict[str, TimeScales]]
    tai_scale: typing.ClassVar[TimeScales]
    tdb_scale: typing.ClassVar[TimeScales]
    tt_scale: typing.ClassVar[TimeScales]
    ut1_scale: typing.ClassVar[TimeScales]
    utc_scale: typing.ClassVar[TimeScales]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

def TAI_to_TT(arg0: float) -> float:
    """Convert time from the TAI scale to the TT scale.
	
	The TAI scale is the International Atomic Time, and the TT scale is the Terrestrial Time.
	
	:param TAI_time:
			Time in seconds since J2000, in the TAI time scale.
	:return:
			Time in seconds since J2000, in the TT time scale.
	"""

def TCB_to_TDB(arg0: float) -> float:
    """Convert time from the TCB scale to the TDB scale.
	
	The TCB scale is the Barycentric Coordinate Time, and the TDB scale is the Barycentric Dynamical Time.
	
	:param TCB_time:
			Time in seconds since J2000, in the TCB time scale.
	:return:
			Time in seconds since J2000, in the TDB time scale.
	"""

def TCG_to_TT(arg0: float) -> float:
    """Convert time from the TCG scale to the TT scale.
	
	The TCG scale is the Geocentric Coordinate Time, and the TT scale is the Terrestrial Time.
	
	:param TCG_time:
			Time in seconds since J2000, in the TCG time scale.
	:return:
			Time in seconds since J2000, in the TT time scale.
	"""

def TDB_to_TCB(arg0: float) -> float:
    """Convert time from the TBD scale to the TCB scale.
	
	The TDB scale is the Barycentric Dynamical Time, and the TCB scale is the Barycentric Coordinate Time.
	
	Inverse function of :func:`TCB_to_TDB`.
	
	
	:param TDB_time:
			Time in seconds since J2000, in the TDB time scale.
	:return:
			Time in seconds since J2000, in the TCB time scale.
	"""

def TT_to_TAI(arg0: float) -> float:
    """Convert time from the TT scale to the TAI scale.
	
	The TT scale is the Terrestrial Time, and the TAI scale is the International Atomic Time.
	
	Inverse function of :func:`TAI_to_TT`.
	
	
	:param TT_time:
			Time in seconds since J2000, in the TT time scale.
	:return:
			Time in seconds since J2000, in the TAI time scale.
	"""

def TT_to_TCG(arg0: float) -> float:
    """Convert time from the TT scale to the TCG scale.
	
	The TT scale is the Terrestrial Time, and the TCG scale is the Geocentric Coordinate Time.
	
	Inverse function of :func:`TCG_to_TT`.
	
	
	:param TT_time:
			Time in seconds since J2000, in the TT time scale.
	:return:
			Time in seconds since J2000, in the TCG time scale.
	"""

def TT_to_TDB_approximate(TT_time: float) -> float:
    """Approximately convert time from the TT scale to the TDB scale.
	
	The TT scale is the Terrestrial Time, and the TDB scale is the Barycentric Dynamical Time.
	
	:param TT_time:
			Time in seconds since J2000, in the TT time scale.
	:return:
			Time in seconds since J2000, in the TDB time scale.
	"""

def add_days_to_datetime(datetime: ..., days_to_add: float) -> ...:
    ...

def add_seconds_to_datetime(datetime: ..., seconds_to_add: float) -> ...:
    ...

def calculate_seconds_in_current_julian_day(julian_day: float) -> float:
    """Determine the number of seconds that have elapsed in the given Julian day.
	
	
	
	:param julian_day:
			Date in Julian days (number of days since January 1st 4713 BC).
	:return:
			Number of seconds that have passed in the given Julian day.
	"""

def calendar_date_to_days_since_epoch(calendar_date: datetime.datetime, days_since_julian_day_zero: float=2451545.0) -> float:
    ...

def calendar_date_to_julian_day(calendar_date: datetime.datetime) -> float:
    """Convert a calendar date to Julian days.
	
	
	
	:param calendar_date:
			Datetime, using the default Python library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored.
	:return:
			Date in Julian days since January 1st 4713 BC.
	"""

def calendar_date_to_julian_day_since_epoch(calendar_date: datetime.datetime, days_since_julian_day_zero: float=2451545.0) -> float:
    ...

def date_time_from_epoch(epoch: float) -> DateTime:
    ...

def date_time_from_iso_string(iso_datetime: float) -> DateTime:
    ...

def datetime_to_python(datetime: ...) -> datetime.datetime:
    ...

def datetime_to_tudat(datetime: datetime.datetime) -> ...:
    ...

def default_time_scale_converter() -> ...:
    ...

def epoch_from_date_time_components(year: int, month: int, day: int, hour: int, minute: int, seconds: float) -> float:
    ...

def epoch_from_date_time_iso_string(iso_datetime: str) -> float:
    ...

def get_days_in_month(month: int, year: int) -> int:
    """Get the number of days in the month of a given year.
	
	
	
	:param month:
			Calendar month.
	:param year:
			Calendar year.
	:return:
			Number of days in the month of the given year.
	"""

def is_leap_year(year: int) -> bool:
    """Assess wether a year is a leap year or not.
	
	
	
	:param year:
			Calendar year.
	:return:
			A value of True means that the year is a leap year.
	"""

def julian_day_to_calendar_date(julian_day: float) -> datetime.datetime:
    """Convert Julian days to a calendar date.
	
	Inverse function of :func:`calendar_date_to_julian_day`.
	
	:param julian_day:
			Date in Julian days since January 1st 4713 BC.
	:return:
			Calendar date using the Datetime Python library, containing the date and time corresponding to the Julian date input.
	"""

def julian_day_to_modified_julian_day(julian_day: float) -> float:
    """Convert a Julian day to a Modified Julian day.
	
	
	
	:param julian_day:
			Date in Julian days (number of days since January 1st 4713 BC).
	:return:
			Date in modified Julian days (number of days since November 17th 1858).
	"""

def julian_day_to_seconds_since_epoch(julian_day: float, days_since_julian_day_zero: float=2451545.0) -> float:
    """Convert Julian days to seconds since a given epoch.
	
	
	
	:param julian_day:
			Date in Julian days since January 1st 4713 BC.
	:param epoch_since_julian_day_zero:
			Epoch since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
	:return:
			Seconds since the Julian date and the given epoch.
	"""

def modified_julian_day_to_julian_day(modified_julian_day: float) -> float:
    """Convert a Modified Julian day to a Julian day.
	
	Inverse function of :func:`julian_day_to_modified_julian_day`.
	
	:param modified_julian_day:
			Date in modified Julian days (number of days since November 17th 1858).
	:return:
			Date in Julian days (number of days since January 1st 4713 BC).
	"""

def seconds_since_epoch_to_julian_centuries_since_epoch(seconds_since_epoch: float) -> float:
    """Convert the number of seconds since a given (unspecified) epoch to Julian centuries since the same epoch.
	
	
	
	:param seconds_since_epoch:
			Seconds elapsed since a given (unspecified) epoch.
	:return:
			Julian centuries since the specified epoch.
	
			Since this is a float, not a integer, meaning that the fraction of the century is also included.
	"""

def seconds_since_epoch_to_julian_day(seconds_since_epoch: float, days_since_julian_day_zero: float=2451545.0) -> float:
    ...

def seconds_since_epoch_to_julian_years_since_epoch(seconds_since_epoch: float) -> float:
    """Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch.
	
	
	
	:param seconds_since_epoch:
			Seconds elapsed since a given (unspecified) epoch.
	:return:
			Julian years since the specified epoch.
	
			Since this is a float, not a integer, meaning that the fraction of the year is also included.
	"""
tai_scale: TimeScales
tdb_scale: TimeScales
tt_scale: TimeScales
ut1_scale: TimeScales
utc_scale: TimeScales