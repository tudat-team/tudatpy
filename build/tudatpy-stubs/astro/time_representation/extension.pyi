import datetime
import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['DateTime', 'TAI_to_TT', 'TCB_to_TDB', 'TCG_to_TT', 'TDB_to_TCB', 'TDB_to_TT', 'TT_to_TAI', 'TT_to_TCG', 'TT_to_TDB', 'TT_to_TDB_approximate', 'Time', 'TimeScaleConverter', 'TimeScales', 'add_days_to_datetime', 'add_seconds_to_datetime', 'calculate_seconds_in_current_julian_day', 'calendar_date_to_days_since_epoch', 'calendar_date_to_julian_day', 'calendar_date_to_julian_day_since_epoch', 'date_time_components_to_epoch', 'date_time_components_to_epoch_time_object', 'date_time_from_epoch', 'date_time_from_iso_string', 'datetime_to_python', 'datetime_to_tudat', 'default_time_scale_converter', 'epoch_from_date_time_components', 'epoch_from_date_time_iso_string', 'get_days_in_month', 'is_leap_year', 'iso_string_to_epoch', 'iso_string_to_epoch_time_object', 'julian_day_to_calendar_date', 'julian_day_to_modified_julian_day', 'julian_day_to_python_datetime', 'julian_day_to_seconds_since_epoch', 'modified_julian_day_to_julian_day', 'python_datetime_to_days_since_epoch', 'python_datetime_to_julian_day', 'seconds_since_epoch_to_julian_centuries_since_epoch', 'seconds_since_epoch_to_julian_day', 'seconds_since_epoch_to_julian_years_since_epoch', 'tai_scale', 'tdb_scale', 'tt_scale', 'ut1_scale', 'utc_scale', 'year_and_days_in_year_to_calendar_date']

class DateTime:
    """Class to store a calendar date and time of day, with high resolution.
    
    Class to store a calendar date and time of day, with high resolution compared to Python datetime.datetime. This class
    stores the seconds as a ``long double`` variable in the C++ implementation, corresponding to about
    16 or 19 digits of precision (depending on the compiler used). In either case, this will be sufficient for sub-femtosecond
    resolution. In addition, this class allows easy conversion to typical time representations in astrodynamics (seconds since J2000,
    Julian day, and modified Julian day)."""

    @staticmethod
    def from_epoch(epoch: float) -> DateTime:
        """
         Creates a Tudat-native :class:`DateTime` object from the seconds since J2000.
        
         Parameters
         ----------
         epoch : float
             Seconds since J2000
        
         Returns
         -------
         DateTime
             Tudat ``DateTime`` object.
        
         Examples
         --------
         In this example, the datetime is constructed from an epoch in seconds since J2000.
         
         .. code-block:: python
         
             from tudatpy.astro.time_representation import DateTime
        
             epoch_et = 788961600.0
        
             dt = DateTime.from_epoch(epoch_et)
             print(dt) # prints 2025-01-01 00:00:00.000000000000000
        """

    @staticmethod
    def from_epoch_time_object(epoch: Time) -> DateTime:
        ...

    @staticmethod
    def from_iso_string(iso_time: str) -> DateTime:
        """
         Creates a Tudat-native :class:`DateTime` object from an ISO datetime string.
        
         Parameters
         ----------
         iso_datetime : str
             Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)
        
         Returns
         -------
         DateTime
             Tudat ``DateTime`` object.
        
         Examples
         --------
         In this example, the datetime is constructed from the iso string.
         
         .. code-block:: python
         
             from tudatpy.astro.time_representation import DateTime
        
             dt = DateTime.from_iso_string("2025-01-01T00:00:00.000")
             print(dt) # prints 2025-01-01 00:00:00.000000000000000
        """

    @staticmethod
    def from_julian_day(julian_day: float) -> DateTime:
        """
         Creates a Tudat-native :class:`DateTime` object from a Julian day.
        
         Parameters
         ----------
         julian_day : float
             Date in Julian days since January 1st 4713 BC.
        
         Returns
         -------
         DateTime
             Tudat ``DateTime`` object.
        
         Examples
         --------
         In this example, the DateTime is constructed from a Julian day since January 1st 4713 BC.
         
         .. code-block:: python
         
             from tudatpy.astro.time_representation import DateTime
        
             julian_day = 2451545.0 
        
             dt = DateTime.from_julian_day(julian_day)
             print(dt) # prints 2000-01-01 12:00:00.000000000000000
        """

    @staticmethod
    def from_modified_julian_day(modified_julian_day: float) -> DateTime:
        """
         Creates a Tudat-native :class:`DateTime` object from a modified Julian day.
        
         Parameters
         ----------
         modified_julian_day : float
             Date in modified Julian days (number of days since November 17th 1858).
        
         Returns
         -------
         DateTime
             Tudat ``DateTime`` object.
        
         Examples
         --------
         In this example, the DateTime is constructed from a modified Julian day since November 17th 1858.
         
         .. code-block:: python
         
             from tudatpy.astro.time_representation import DateTime
        
             modified_julian_day = 51544.5
        
             dt = DateTime.from_modified_julian_day(modified_julian_day)
             print(dt) # prints 2000-01-01 12:00:00.000000000000000
        """

    @staticmethod
    def from_python_datetime(datetime: datetime.datetime) -> DateTime:
        """
        Function to convert a Python `datetime.datetime` object to a Tudat :class:`DateTime` object. The Tudat-native alternative has the advantage of providing sub-femtosecond resolution, as opposed to the microsecond resolution of the Python version.
        
        .. warning::
        
            This function uses the C++ `std::chrono` library, which is limited in the time range it can represent. If the range is exceeded, the conversion will overflow and **NOT** throw an exception.
        
            The exact range is platform-dependent. On Windows, dates between 1970-01-01 and 3000-12-31 are allowed, on MacOS dates after 1900-01-01 are allowed and on Linux dates between 1678-01-01 and 2261-12-31 are allowed.
        
        Parameters
        ----------
        datetime : datetime.datetime
            Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
        Returns
        -------
        DateTime
            DateTime object defined in Tudat
        
        Examples
        --------
        In this example, the Tudat DateTime object is constructed from python native datetime object.
        
        .. code-block:: python
        
            from datetime import datetime
            from tudatpy.astro.time_representation import DateTime
        
            python_datetime = datetime(2025, 1, 1, 0, 0, 0)
        
            dt = DateTime.from_python_datetime(python_datetime)
            print(dt) # prints 2025-01-01 00:00:00.000000000000000
        """

    @staticmethod
    def from_year_and_day_of_year(year: int, day_of_year: int) -> DateTime:
        """
        Create the Tudat :class:`DateTime` from the year and the number of days in the year.
        
        Parameters
        ----------
        year : int
            Calendar year.
        day_of_year : int
            Number of days that have passed in the year.
        Returns
        -------
        DateTime
            Corresponding calendar date as a :class:`DateTime` object. Note: the hours, minutes and seconds in the object are set to 0 when calling this function.
        
        Examples
        --------
        In this example, the calendar date corresponding to when 122 days have passed in 2020 is computed.
        
        .. code-block:: python
        
            # Compute the calendar date when 122 days have passed in 2020
            currentDate = time_representation.DateTime.from_year_and_day_of_year(2020, 122)
            # Print the converted output
            print(currentDate)  # prints (2020, 5, 2, 0, 0)
        """

    def __init__(self, year: int, month: int, day: int, hour: int=12, minute: int=0, seconds: float=0.0) -> None:
        ...

    def __repr__(self) -> str:
        ...

    def __str__(self) -> str:
        ...

    def add_days(self, days_to_add: Time) -> DateTime:
        """
         Function to create a new Tudat :class:`DateTime` object by adding a number of days (86400 seconds) to an existing Tudat :class:`DateTime` object
        
         .. note::
            
            This method does not modify the original :class:`DateTime` object, but returns a new one with the added days.
        
         Parameters
         ----------
         days_to_add : float
             Number of days to add
         Returns
         -------
         DateTime
             Tudat-native Datetime object created by adding the given number of days to the original DateTime
         Examples
         --------
         In this example, 1 day is added to a DateTime object to construct a new DateTime.
         
         .. code-block:: python
         
             from tudatpy.astro.time_representation import DateTime
        
             dt = DateTime(2025, 1, 1, 0, 0, 0.0)
             dt_days_added = dt.add_days(1.0)
             print(f"Original dt: {dt}")
             print(f"dt with days added: {dt_days_added}")
             # prints:
             # Original dt: 2025-01-01 00:00:00.000000000000000
             # dt with days added: 2025-01-02 00:00:00.000000000000000   
        """

    def add_seconds(self, seconds_to_add: Time) -> DateTime:
        """
         Function to create a new Tudat :class:`DateTime` object by adding a number of seconds to an existing Tudat :class:`DateTime` object.
        
         .. note::
            
            This method does not modify the original :class:`DateTime` object, but returns a new one with the added seconds.
        
         Parameters
         ----------
         seconds_to_add : float
             Number of seconds to add
         Returns
         -------
         DateTime
             Tudat-native Datetime object created by adding the given number of seconds to the original DateTime
         
         Examples
         --------
         In this example, 86400 seconds are added to a DateTime object to construct a new DateTime.
         
         .. code-block:: python
         
             from tudatpy.astro.time_representation import DateTime
        
             dt = DateTime(2025, 1, 1, 0, 0, 0.0)
             dt_seconds_added = dt.add_seconds(86400.0)
             print(f"Original dt: {dt}")
             print(f"dt with seconds added: {dt_seconds_added}")
             # prints:
             # Original dt: 2025-01-01 00:00:00.000000000000000
             # dt with seconds added: 2025-01-02 00:00:00.000000000000000   
        """

    def day_of_year(self) -> int:
        """
         .. warning::
        
            This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_day_of_year` instead.
        
         Function to get the day number in the current year
        
        
         Returns
         -------
         int
             Day number in the current year
        """

    def epoch(self) -> float:
        """
         .. warning::
        
            This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_epoch` instead.
        
         Function to get the epoch in seconds since J2000 for the current date and time
        
        
         Returns
         -------
         float
             Current epoch in seconds since J2000
        """

    def epoch_time_object(self) -> Time:
        ...

    def iso_string(self, add_T: bool=False, number_of_digits_seconds: int=15) -> str:
        """
         .. warning::
        
            This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_iso_string` instead.
        
         Function to get the ISO-compatible string.
        
        
         Function to get the current date and time as an ISO-compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..") where the seconds may be provided with any number of digits. The 'T' entry separating the date from the time may be omitted by setting the ``add_T`` parameter to false
        
        
         Parameters
         ----------
         add_T : bool
            Boolean denoting whether to use a 'T' or a blank space to separate the date from the time
         number_of_digits_seconds : int, default = 15
            Number of digits to use after the decimal separator (trailing zeros will be truncated)
        
         Returns
         -------
         str
             ISO-compatible string representing the date and time
        """

    def julian_day(self) -> float:
        """
         .. warning::
        
            This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_julian_day` instead.
        
         Function to get the epoch as Julian day for the current date and time
        
        
         Returns
         -------
         float
             Current Julian day
        """

    def modified_julian_day(self) -> float:
        """
         .. warning::
        
            This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_modified_julian_day` instead.
        
         Function to get the epoch as modified Julian day for the current date and time
        
        
         Returns
         -------
         float
             Current modified Julian day
        """

    def to_day_of_year(self) -> int:
        """
         Function to get the day number in the current year
        
        
         Returns
         -------
         int
             Day number in the current year
        """

    def to_days_since_reference_julian_day(self, reference_julian_day: float=2451545.0) -> float:
        """
         Convert a DateTime to Julian days since a given reference Julian date.
        
        
         Parameters
         ----------
         reference_julian_day : float, default = constants.JULIAN_DAY_ON_J2000
             Reference epoch (in days) since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0) corresponding to the 1st of January 2000.
         Returns
         -------
         float
             Date in Julian days since the given reference epoch.
        
         Examples
         --------
         In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days since J2000 (the 1st of January 2000).
        
         .. code-block:: python
        
           from tudatpy.astro.time_representation import DateTime
        
           # Define the calendar date using datetime
           dt = DateTime(2022, 5, 21, 13, 52, 41)
           # Convert the calendar date to Julian days since J2000
           julian_date = dt.to_days_since_reference_julian_day()
           # Print the converted output
           print(julian_date)  # prints 8176.07825231459
        """

    def to_epoch(self) -> float:
        """
         Function to get the epoch in seconds since J2000 for the current date and time
        
        
         Returns
         -------
         float
             Current epoch in seconds since J2000
        """

    def to_epoch_time_object(self) -> Time:
        ...

    def to_iso_string(self, add_T: bool=False, number_of_digits_seconds: int=15) -> str:
        """
         Function to get the ISO-compatible string.
        
        
         Function to get the current date and time as an ISO-compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..") where the seconds may be provided with any number of digits. The 'T' entry separating the date from the time may be omitted by setting the ``add_T`` parameter to false
        
        
         Parameters
         ----------
         add_T : bool
         Boolean denoting whether to use a 'T' or a blank space to separate the date from the time
        
         number_of_digits_seconds : int, default = 15
         Number of digits to use after the decimal separator (trailing zeros will be truncated)
        
         Returns
         -------
         str
             ISO-compatible string representing the date and time
        """

    def to_julian_day(self) -> float:
        """
         Function to get the epoch as Julian day for the current date and time
        
        
         Returns
         -------
         float
             Current Julian day
        """

    def to_modified_julian_day(self) -> float:
        """
         Function to get the epoch as modified Julian day for the current date and time
        
        
         Returns
         -------
         float
             Current modified Julian day
        """

    def to_python_datetime(self) -> datetime.datetime:
        """
        Method to convert retrieve a Python datetime.datetime object from the Tudat :class:`DateTime` object. This is the inverse of the :meth:`~tudatpy.astro.time_representation.DateTime.from_python_datetime` method.
        
        .. note::
        
            The conversion uses the C++ `std::chrono` library, which is limited the time range it can represent. If the range is exceeded, the conversion will fail and throw an exception.
            
            The exact range is platform-dependent. On Windows, dates between 1970-01-01 and 3000-12-31 are allowed, on MacOS dates after 1900-01-01 are allowed and on Linux dates between 1678-01-01 and 2261-12-31 are allowed.
        
        Returns
        -------
        datetime.datetime
            Datetime object, using the Python datetime library
        """

    @property
    def day(self) -> int:
        """
         Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month
        
        
         :type: int
        """

    @day.setter
    def day(self, arg1: int) -> None:
        ...

    @property
    def hour(self) -> int:
        """
         Full hours into the current day (value must be 0-23)
        
        
         :type: int
        """

    @hour.setter
    def hour(self, arg1: int) -> None:
        ...

    @property
    def minute(self) -> int:
        """
         Full minutes into the current hour (value must be 0-59)
        
        
         :type: int
        """

    @minute.setter
    def minute(self, arg1: int) -> None:
        ...

    @property
    def month(self) -> int:
        """
         Calendar month (value must be 1-12)
        
        
         :type: int
        """

    @month.setter
    def month(self, arg1: int) -> None:
        ...

    @property
    def seconds(self) -> float:
        """
         Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.
        
        
         :type: float
        """

    @seconds.setter
    def seconds(self, arg1: float) -> None:
        ...

    @property
    def year(self) -> int:
        """
         Calendar year
        
        
         :type: int
        """

    @year.setter
    def year(self, arg1: int) -> None:
        ...

class Time:
    """
    Class for defining time with a resolution that is sub-femtosecond for very long periods of time.
    
    Using double or long double precision as a representation of time, the issue of reduced precision will 
    occur over long time periods. For instance, over a period of 10^8 seconds (about 3 years), double and 
    long double representations have resolution of about 10^-8 and 10^-11 s respectively, which is 
    insufficient for various applications. 
    
    This class uses an integer to represent the number of hours since an epoch, and a long double to 
    represent the number of seconds into the present hour. This provides a resolution of < 1 femtosecond, 
    over a range of 2147483647 hours (about 300,000 years), which is more than sufficient for practical 
    applications.
    
    The Time class supports standard arithmetic operations (+, -, *, /) with Time objects and floats, comparison operations, and 
    automatic conversion to/from floating-point types.
        """

    @typing.overload
    def __add__(self, arg0: Time) -> Time:
        ...

    @typing.overload
    def __add__(self, arg0: float) -> Time:
        ...

    @typing.overload
    def __eq__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __eq__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __eq__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __eq__(self, arg0: float) -> bool:
        ...

    def __float__(self) -> float:
        ...

    @typing.overload
    def __ge__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __ge__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __ge__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __gt__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __gt__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __gt__(self, arg0: float) -> bool:
        ...

    def __hash__(self) -> int:
        ...

    @typing.overload
    def __iadd__(self, arg0: Time) -> None:
        ...

    @typing.overload
    def __iadd__(self, arg0: float) -> None:
        ...

    def __imul__(self, arg0: float) -> None:
        ...

    @typing.overload
    def __init__(self, seconds_since_j2000: float) -> None:
        """
             Create a Time object from seconds since J2000.
             
             Parameters
             ----------
             seconds_since_j2000 : float
                 Number of seconds since J2000 epoch
             
             Returns
             -------
             Time
                 Time object initialized to specified seconds since J2000
             
             Examples
             --------
             >>> from tudatpy.kernel import Time
             >>> t = Time(3600.0)  # 1 hour after J2000
        """

    @typing.overload
    def __init__(self, full_periods: int, seconds_into_full_period: float) -> None:
        """
             Create a Time object from full periods (hours) and seconds into the current period.
             
             Parameters
             ----------
             full_periods : int
                 Number of full hours since epoch
             seconds_into_full_period : float
                 Number of seconds into current hour. Need not be in range [0, 3600];
                 the time representation is normalized automatically.
             
             Returns
             -------
             Time
                 Time object initialized to specified time
             
             Examples
             --------
             >>> from tudatpy.kernel import Time
             >>> t = Time(2, 1800.0)  # 2.5 hours after epoch
        """

    @typing.overload
    def __isub__(self, arg0: Time) -> None:
        ...

    @typing.overload
    def __isub__(self, arg0: float) -> None:
        ...

    def __itruediv__(self, arg0: float) -> None:
        ...

    @typing.overload
    def __le__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __le__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __le__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __lt__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __lt__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __lt__(self, arg0: float) -> bool:
        ...

    def __mul__(self, arg0: float) -> Time:
        ...

    @typing.overload
    def __ne__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __ne__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __ne__(self, arg0: float) -> bool:
        ...

    def __radd__(self, arg0: float) -> Time:
        ...

    def __rmul__(self, arg0: float) -> Time:
        ...

    def __rsub__(self, arg0: float) -> Time:
        ...

    @typing.overload
    def __sub__(self, arg0: Time) -> Time:
        ...

    @typing.overload
    def __sub__(self, arg0: float) -> Time:
        ...

    def __truediv__(self, arg0: float) -> Time:
        ...

    def to_float(self) -> float:
        """
            Converts the time to a float (double) representing seconds since J2000.
                
            Returns
            -------
            float
                Number of seconds since J2000
            
            Examples
            --------
            In this example, a Time object is converted back to seconds since J2000.
            
            .. code-block:: python
            
                from tudatpy.kernel import Time
                
                # Create a Time object from seconds since J2000
                t = Time(3600.0)  # 1 hour after J2000
                
                # Convert back to seconds
                seconds = t.to_float()
                print(seconds)  # prints 3600.0
        """

class TimeScaleConverter:
    """Class to convert between different time scales (TAI, TT, TDB, UTC, UT1)
    
    Class to convert between different time scales (TAI, TT, TDB, UTC, UT1), as per algorithms described in (for instance) IERS 2010 Conventions and USNO circular no. 179.
    The algorithms used for the conversion are (where there is any choice in models):
    
    * Conversion between TDB and TT uses the SOFA function ``iauDtdb`` function (equivalently to :func:`TT_to_TDB` and :func:`TDB_to_TT`)
    * Conversion between UTC and UT1 applies (semi-)diurnal and daily measured variations depending on settings during object creation (typically as per :func:`~default_time_scale_converter`)
    * Leap seconds as per the latest SOFA package"""

    def convert_time(self, input_scale: TimeScales, output_scale: TimeScales, input_value: float, earth_fixed_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]=...) -> float:
        ...

    def convert_time_object(self, input_scale: TimeScales, output_scale: TimeScales, input_value: Time, earth_fixed_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]=...) -> Time:
        ...

    def get_time_difference(self, input_scale: TimeScales, output_scale: TimeScales, input_value: float, earth_fixed_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]=...) -> float:
        ...

    def get_time_object_difference(self, input_scale: TimeScales, output_scale: TimeScales, input_value: Time, earth_fixed_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]=...) -> Time:
        ...

class TimeScales:
    """Enumeration of available time scales between which the :class:`~TimeScaleConverter` can automaticaly convert.
    
     
    
    Members:
    
      tai_scale : 
     
    
      tt_scale : 
     
    
      tdb_scale : 
     
    
      utc_scale : 
     
    
      ut1_scale : 
     """
    __members__: typing.ClassVar[dict[str, TimeScales]]
    tai_scale: typing.ClassVar[TimeScales]
    tdb_scale: typing.ClassVar[TimeScales]
    tt_scale: typing.ClassVar[TimeScales]
    ut1_scale: typing.ClassVar[TimeScales]
    utc_scale: typing.ClassVar[TimeScales]

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

def TAI_to_TT(TAI_time: float) -> float:
    """Convert time from the TAI scale to the TT scale.
    
    The TAI scale is the International Atomic Time, and the TT scale is the Terrestrial Time.
    
    Parameters
    ----------
    TAI_time : float
        Time in seconds since J2000, in the TAI time scale.
    Returns
    -------
    float
        Time in seconds since J2000, in the TT time scale."""

def TCB_to_TDB(TCB_time: float) -> float:
    """Convert time from the TCB scale to the TDB scale.
    
    The TCB scale is the Barycentric Coordinate Time, and the TDB scale is the Barycentric Dynamical Time.
    
    Parameters
    ----------
    TCB_time : float
        Time in seconds since J2000, in the TCB time scale.
    Returns
    -------
    float
        Time in seconds since J2000, in the TDB time scale.
    
    
    
    
    
    Examples
    --------
    In this example, the calendar date of the 17th of February 2022, at 15:41 and 2 seconds is first converted to Julian seconds since J2000.
    Then, this date and time is converted from the TCB scale to the TDB scale.
    
    .. code-block:: python
    
      # Define the date and time
      date = datetime.datetime(2022, 2, 17, 15, 41, 2)
      # Convert it in Julian days since J2000
      date_J2000 = time_representation.python_datetime_to_julian_day(date)
      # Convert it in Julian seconds since J2000
      date_J2000_sec = time_representation.julian_day_to_seconds_since_epoch(date_J2000)
      # Check the date from the TCB scale to the TDB scale
      date_TDB_scale = time_representation.TCB_to_TDB(date_J2000_sec)
      # Print the converted output
      print(date_TDB_scale)  # prints 698384439.9176273"""

def TCG_to_TT(TCG_time: float) -> float:
    """Convert time from the TCG scale to the TT scale.
    
    The TCG scale is the Geocentric Coordinate Time, and the TT scale is the Terrestrial Time.
    
    Parameters
    ----------
    TCG_time : float
        Time in seconds since J2000, in the TCG time scale.
    Returns
    -------
    float
        Time in seconds since J2000, in the TT time scale."""

def TDB_to_TCB(TDB_time: float) -> float:
    """Convert time from the TBD scale to the TCB scale.
    
    The TDB scale is the Barycentric Dynamical Time, and the TCB scale is the Barycentric Coordinate Time.
    
    Inverse function of :func:`TCB_to_TDB`.
    
    
    Parameters
    ----------
    TDB_time : float
        Time in seconds since J2000, in the TDB time scale.
    Returns
    -------
    float
        Time in seconds since J2000, in the TCB time scale."""

def TDB_to_TT(TDB_time: Time, earth_fixed_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]) -> Time:
    """Convert time from the TDB scale to the TT scale.
    
    Convert time from the TT scale to the TDB scale, using the `iauDtdb <https://www2.mpia-hd.mpg.de/~mathar/progs/sofa_api/group__SR.html#gaeab39417bb16e66c570232102a055f2d>`_ function in Sofa, which has an accuracy of 3 nanoseconds or better in the time span 1950-2050
    To call the Sofa function, we assume UTC=UT1 and TT=TDB to compute UT1 from the input TDB. For this function, UT1 is only used in this sofa function to compute the local solar time, which in turn is used to compute the time-modulation of the topocentric term
    (dependent on the ground station position). Since these terms have an amplitude at the microsecond level, the error induced by this assumption is negligible given the inherent quality of the Sofa model.
    
    Parameters
    ----------
    TDB_time : float
        Time in seconds since J2000, in the TDB time scale.
    earth_fixed_position : numpy.ndarray, default=numpy.array([0, 0, 0])
        Earth-fixed position (e.g. in ITRF) that is used for detailed conversion between TDB and TT (induces a signature at the microsecond level)
    Returns
    -------
    float
        Time in seconds since J2000, in the TDB time scale."""

def TT_to_TAI(TT_time: float) -> float:
    """Convert time from the TT scale to the TAI scale.
    
    The TT scale is the Terrestrial Time, and the TAI scale is the International Atomic Time.
    
    Inverse function of :func:`TAI_to_TT`.
    
    
    Parameters
    ----------
    TT_time : float
        Time in seconds since J2000, in the TT time scale.
    Returns
    -------
    float
        Time in seconds since J2000, in the TAI time scale."""

def TT_to_TCG(TT_time: float) -> float:
    """Convert time from the TT scale to the TCG scale.
    
    The TT scale is the Terrestrial Time, and the TCG scale is the Geocentric Coordinate Time.
    
    Inverse function of :func:`TCG_to_TT`.
    
    
    Parameters
    ----------
    TT_time : float
        Time in seconds since J2000, in the TT time scale.
    Returns
    -------
    float
        Time in seconds since J2000, in the TCG time scale."""

def TT_to_TDB(TT_time: Time, earth_fixed_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]) -> Time:
    """Convert time from the TT scale to the TDB scale.
    
    Convert time from the TT scale to the TDB scale, using the `iauDtdb <https://www2.mpia-hd.mpg.de/~mathar/progs/sofa_api/group__SR.html#gaeab39417bb16e66c570232102a055f2d>`_ function in Sofa, which has an accuracy of 3 nanoseconds or better in the time span 1950-2050
    To call the Sofa function, we assume UTC=UT1 for the UT1 input. For this function, UT1 is only used in this sofa function to compute the local solar time, which in turn is used to compute the time-modulation of the topocentric term
    (dependent on the ground station position). Since these terms have an amplitude at the microsecond level, the error induced by this assumption is negligible given the inherent quality of the Sofa model.
    
    Parameters
    ----------
    TT_time : float
        Time in seconds since J2000, in the TT time scale.
    earth_fixed_position : numpy.ndarray, default=numpy.array([0, 0, 0])
        Earth-fixed position (e.g. in ITRF) that is used for detailed conversion between TDB and TT (induces a signature at the microsecond level)
    Returns
    -------
    float
        Time in seconds since J2000, in the TDB time scale."""

def TT_to_TDB_approximate(TT_time: float) -> float:
    """Approximately convert time from the TT scale to the TDB scale, using the following first-order approximation:
    
    .. math::
        t_{\\text{TDB}} = t_{\\text{TT}} + 0.001657  * \\sin( 628.3076 T_{\\text{TT}} + 6.2401 );
    
    with :math:`t` the epoch in seconds since J2000, :math:`T` the epoch in centuries since J2000.
    
    The TT scale is the Terrestrial Time, and the TDB scale is the Barycentric Dynamical Time.
    
    Parameters
    ----------
    TT_time : float
        Time in seconds since J2000, in the TT time scale.
    Returns
    -------
    float
        Time in seconds since J2000, in the TDB time scale."""

def add_days_to_datetime(datetime: DateTime, days_to_add: Time) -> DateTime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.add_days()` instead.
    
    Function to create a new Tudat :class:`DateTime` object by adding a number of days (86400 seconds) to an existing Tudat :class:`DateTime` object
    
    
    Parameters
    ----------
    datetime : DateTime
        Tudat-native Datetime object to which a number of days are to be added
    days_to_add : float
        Number of days to add
    Returns
    -------
    DateTime
        Tudat-native Datetime object created by adding the given number of days to the original DateTime"""

def add_seconds_to_datetime(datetime: DateTime, seconds_to_add: Time) -> DateTime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.add_seconds()` instead.
    
    Function to create a new Tudat :class:`DateTime` object by adding a number of seconds to an existing Tudat :class:`DateTime` object
    
    
    Parameters
    ----------
    datetime : DateTime
        Tudat-native Datetime object to which a number of seconds are to be added
    seconds_to_add : float
        Number of seconds to add
    Returns
    -------
    DateTime
        Tudat-native Datetime object created by adding the given number of seconds to the original DateTime"""

def calculate_seconds_in_current_julian_day(julian_day: float) -> float:
    """Determine the number of seconds that have elapsed in the given Julian day.
    
    
    Parameters
    ----------
    julian_day : float
        Date in Julian days (number of days since January 1st 4713 BC).
    Returns
    -------
    float
        Number of seconds that have passed in the given Julian day.
    
    
    
    
    
    Examples
    --------
    In this example, the number of seconds that have elapsed at the Julian day `2451545.2` is computed.
    
    .. code-block:: python
    
      # Compute the number of seconds that have passed in the given Julian day
      seconds_passed = time_representation.calculate_seconds_in_current_julian_day(constants.JULIAN_DAY_ON_J2000)
      # Print the converted output
      print(seconds_passed)  # prints 43200.0"""

def calendar_date_to_days_since_epoch(calendar_date: datetime.datetime, days_since_julian_day_zero: float=2451545.0) -> float:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_days_since_reference_julian_day()` instead.
    
    Convert a calendar date to Julian days since a given epoch.
    
    
    Parameters
    ----------
    calendar_date : datetime.datetime
        Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored.
    days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
        Reference epoch (in days) since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0) corresponding to the 1st of January 2000.
    Returns
    -------
    float
        Date in Julian days since the given epoch.
    
    
    
    
    
    Examples
    --------
    In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days since J2000 (the 1st of January 2000).
    
    .. code-block:: python
    
      # Define the calendar date using datetime
      calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
      # Convert the calendar date to Julian days since J2000
      julian_date = time_representation.calendar_date_to_days_since_epoch(calendar_date)
      # Print the converted output
      print(julian_date)  # prints 8176.07825231459"""

def calendar_date_to_julian_day(calendar_date: datetime.datetime) -> float:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_julian_day()` instead.
    
    Convert a calendar date to Julian days.
    
    
    Parameters
    ----------
    calendar_date : datetime.datetime
        Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
    Returns
    -------
    float
        Julian day number (days since noon January 1st 4713 BC.)
    
    
    
    
    
    Examples
    --------
    In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days.
    
    .. code-block:: python
    
      # Define the calendar date using datetime
      calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
      # Convert the calendar date to Julian days since January 1st 4713 BC
      julian_date = time_representation.calendar_date_to_julian_day(calendar_date)
      # Print the converted output
      print(julian_date)  # prints 2459721.0782523146"""

def calendar_date_to_julian_day_since_epoch(calendar_date: datetime.datetime, days_since_julian_day_zero: float=2451545.0) -> float:
    ...

def date_time_components_to_epoch(year: int, month: int, day: int, hour: int, minute: int, seconds: float) -> float:
    """Computes the epoch as seconds since J2000 from the entries of the current date and time.
    
    Computes the epoch as seconds since J2000. This function is added for convenience, and is equivalent to ``DateTime(...).to_epoch()``.
    
    Parameters
    ----------
    year : int
        Calendar year
    
    month : int
        Calendar month (value must be 1-12)
    
    day : int
        Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month
    
    hour : int
        Full hours into the current day (value must be 0-23)
    
    minute : int
        Full minutes into the current hour (value must be 0-59)
    
    seconds : float
        Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.
    
    Returns
    -------
    float
        Time in seconds since J2000."""

def date_time_components_to_epoch_time_object(year: int, month: int, day: int, hour: int, minute: int, seconds: float) -> Time:
    ...

def date_time_from_epoch(epoch: Time) -> DateTime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_epoch` instead.
    
    
    Creates a Tudat-native :class:`DateTime` object from the seconds since J2000.
    
    
    Parameters
    ----------
    epoch : float
        Seconds since J2000
    
    Returns
    -------
    DateTime
        Tudat ``DateTime`` object."""

def date_time_from_iso_string(iso_string: str) -> DateTime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_iso_string` instead.
    
    Creates a Tudat-native :class:`DateTime` object from an ISO datetime string.
    
    
    Parameters
    ----------
    iso_datetime : str
        Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)
    
    Returns
    -------
    DateTime
        Tudat ``DateTime`` object."""

def datetime_to_python(datetime: DateTime) -> datetime.datetime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.to_python_datetime` instead.
    
    
    Function to convert a Tudat :class:`DateTime` object to a Python datetime.datetime object. This is the inverse of the :func:`datetime_to_tudat` function
    
    Parameters
    ----------
    datetime : DateTime
        Tudat-native Datetime object. Both the date and the time (hour, minutes, and seconds), can be specified, up to sub-femtosecond resolution.
    Returns
    -------
    datetime.datetime
        Datetime object, using the Python datetime library"""

def datetime_to_tudat(datetime: datetime.datetime) -> DateTime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_python_datetime` instead.
    
    Function to convert a Python datetime.datetime object to a Tudat :class:`DateTime` object. The Tudat-native alternative has the advantage of providing sub-femtosecond resolution, as opposed to the microsecond resolution of the Python version
    
    
    Parameters
    ----------
    datetime : datetime.datetime
        Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
    Returns
    -------
    DateTime
        DateTime object defined in Tudat"""

def default_time_scale_converter() -> TimeScaleConverter:
    """Function to create a time-scale converter object with default settings. In particular, it uses default settings for conversion between UT1 and UTC:
    
    * Corrections for semi-diurnal variations due to libration for a non-rigid Earth as per Table 5.1b of IERS Conventions 2010
    * Corrections diurnal and semidiurnal variations due to ocean tides as per Tables 8.2a and 8.2b of the IERS Conventions 2010
    * For epoch 01-01-1962 and later: linear interpolation (correcting for discontinuities during days with leap seconds) of daily corrections for UTC-UT1 from the eopc04_14_IAU2000.62-now.txt file in the tudat-resources directory
    * For epochs before 01-01-1962, where UTC-UT1 is not available from these files, we use the values if for :math:\\Delta T = UT1-TT` from `here <https://webspace.science.uu.nl/~gent0113/deltat/deltat.htm>`_ (back to year 1620). In this period, we set the approxomation UTC=UT1
    
    See :class:`~TimeScaleConverter` for specific functionality and options for time-scale conversions.
    
    Returns
    -------
    TimeScaleConverter
        Object to convert between different terrestrial time scales."""

def epoch_from_date_time_components(year: int, month: int, day: int, hour: int, minute: int, seconds: float) -> Time:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`date_time_components_to_epoch` instead.
    
    Computes the epoch as seconds since J2000 from the entries of the current date and time.
    
    Computes the epoch as seconds since J2000. This function is added for convenience, and creates a :class:`DateTime` object, and subsequently calls its ``epoch`` function
    
    Parameters
    ----------
    year : int
        Calendar year
    
    month : int
        Calendar month (value must be 1-12)
    
    day : int
        Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month
    
    hour : int
        Full hours into the current day (value must be 0-23)
    
    minute : int
        Full minutes into the current hour (value must be 0-59)
    
    seconds : float
        Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.
    
    Returns
    -------
    float
        Time in seconds since J2000."""

def epoch_from_date_time_iso_string(iso_datetime: str) -> Time:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use :func:`iso_string_to_epoch` instead.
    
    Computes the epoch as seconds since J2000 from an ISO datetime string.
    
    Computes the epoch as seconds since J2000. This function is added for convenience, and creates a :class:`DateTime` object, and subsequently calls its ``epoch`` function
    
    Parameters
    ----------
    iso_datetime : str
        Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)
    
    Returns
    -------
    float
        Time in seconds since J2000."""

def get_days_in_month(month: int, year: int) -> int:
    """Get the number of days in the month of a given year.
    
    
    Parameters
    ----------
    month : int
        Calendar month.
    year : int
        Calendar year.
    Returns
    -------
    int
        Number of days in the month of the given year.
    
    
    
    
    
    Examples
    --------
    In this example, the number of days in February for both 2021 and 2020 are computed.
    
    .. code-block:: python
    
      # Check the number of days in February 2021
      days_feb_2021 = time_representation.get_days_in_month(2, 2021)
      # Print the converted output
      print(days_feb_2021)  # prints 28
      # Check the number of days in February 2022
      days_feb_2020 = time_representation.get_days_in_month(2, 2020)
      # Print the converted output
      print(days_feb_2020)  # prints 29"""

def is_leap_year(year: int) -> bool:
    """Assess wether a year is a leap year or not.
    
    
    Parameters
    ----------
    year : int
        Calendar year.
    Returns
    -------
    bool
        A value of True means that the year is a leap year.
    
    
    
    
    
    Examples
    --------
    In this example, the first list should contains only `True`, and the second `False`, since the first list uses leap years and the second does not.
    
    .. code-block:: python
    
      # Check known leap years
      leap_years = [time_representation.is_leap_year(year) for year in [2020, 2016, 2000, 2400]]
      # Print the converted output
      print(leap_years)  # prints [True, True, True, True]
      # Check known non-leap years
      non_leap_years = [time_representation.is_leap_year(year) for year in [2021, 2022, 2100, 2001]]
      # Print the converted output
      print(non_leap_years)  # prints [False, False, False, False]"""

def iso_string_to_epoch(iso_datetime: str) -> float:
    """Computes the epoch as seconds since J2000 from an ISO datetime string.
    
    Computes the epoch as seconds since J2000. This function is added for convenience and is equivalent to ``DateTime.from_iso_string(...).to_epoch()``.
    
    Parameters
    ----------
    iso_datetime : str
        Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)
    
    Returns
    -------
    float
        Time in seconds since J2000."""

def iso_string_to_epoch_time_object(iso_datetime: str) -> Time:
    ...

def julian_day_to_calendar_date(julian_day: float) -> datetime.datetime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_julian_day(...).to_python_datetime()` instead.
    
    Convert Julian days to a calendar date.
    
    Inverse function of :func:`calendar_date_to_julian_day`.
    
    Parameters
    ----------
    julian_day : float
        Date in Julian days since January 1st 4713 BC.
    Returns
    -------
    datetime.datetime
        Datetime object, using the Python datetime library, containing the date and time corresponding to the Julian date input.
    
    
    
    
    
    Examples
    --------
    In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to a calendar date.
    
    .. code-block:: python
    
      # Define the Julian date in days since January 1st 4713 BC
      julian_date = 2459721.0783
      # Convert the Julian date to a calendar date
      calendar_date = time_representation.julian_day_to_calendar_date(julian_date)
      # Print the converted output
      print(calendar_date)  # prints datetime.datetime(2022, 5, 21, 13, 52, 45)"""

def julian_day_to_modified_julian_day(julian_day: float) -> float:
    """Convert a Julian day to a Modified Julian day.
    
    
    Parameters
    ----------
    julian_day : float
        Date in Julian days (number of days since January 1st 4713 BC).
    Returns
    -------
    float
        Date in modified Julian days (number of days since November 17th 1858).
    
    
    
    
    
    Examples
    --------
    In this example, the Julian date `2451545.0` (J2000) is converted to a modified Julian date.
    
    .. code-block:: python
    
      # Convert from Julian Days to Modified Julian Days
      MJD = time_representation.julian_day_to_modified_julian_day(constants.JULIAN_DAY_ON_J2000)
      # Print the converted output
      print(MJD)  # prints 51544.5"""

def julian_day_to_python_datetime(julian_day: float) -> datetime.datetime:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_julian_day(...).to_python_datetime()` instead.
    
    
    Convert Julian days to a calendar date.
    
    Inverse function of :func:`python_datetime_to_julian_day`.
    
    Parameters
    ----------
    julian_day : float
        Date in Julian days since January 1st 4713 BC.
    Returns
    -------
    datetime.datetime
        Datetime object, using the Python datetime library, containing the date and time corresponding to the Julian date input.
    
    
    
    
    
    Examples
    --------
    In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to a calendar date.
    
    .. code-block:: python
    
      # Define the Julian date in days since January 1st 4713 BC
      julian_date = 2459721.0783
      # Convert the Julian date to a calendar date
      calendar_date = time_representation.julian_day_to_python_datetime(julian_date)
      # Print the converted output
      print(calendar_date)  # prints datetime.datetime(2022, 5, 21, 13, 52, 45)"""

def julian_day_to_seconds_since_epoch(julian_day: float, days_since_julian_day_zero: float=2451545.0) -> float:
    """Convert Julian days to seconds since a given epoch.
    
    
    Parameters
    ----------
    julian_day : float
        Date in Julian days since January 1st 4713 BC.
    days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
        Reference epoch (in days since January 1st 4713 BC) since when the number of seconds have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
    Returns
    -------
    float
        Seconds since the Julian date and the given epoch.
    
    
    
    
    
    Examples
    --------
    In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to seconds since J2000 (January 1st 2000).
    
    .. code-block:: python
    
      # Define the Julian date in days since January 1st 4713 BC
      julian_date = 2459721.0783
      # Convert the Julian date to the number of seconds since J2000
      seconds_since_J2000 = time_representation.julian_day_to_seconds_since_epoch(julian_date)
      # Print the converted output
      print(seconds_since_J2000)  # prints 706413165.1200145"""

def modified_julian_day_to_julian_day(modified_julian_day: float) -> float:
    """Convert a Modified Julian day to a Julian day.
    
    Inverse function of :func:`julian_day_to_modified_julian_day`.
    
    Parameters
    ----------
    modified_julian_day : float
        Date in modified Julian days (number of days since November 17th 1858).
    Returns
    -------
    float
        Date in Julian days (number of days since January 1st 4713 BC).
    
    
    
    
    
    Examples
    --------
    In this example, the Modified Julian date `51544.5` ( corresponding to J2000) is converted to a modified Julian date.
    
    .. code-block:: python
    
      # Define J2000 in Modified Julian Days
      J2000_MJD = 51544.5
      # Convert from Modified Julian Days to Julian Days
      J2000 = time_representation.modified_julian_day_to_julian_day(J2000_MJD)
      # Print the converted output
      print(J2000)  # prints 2451545.0"""

def python_datetime_to_days_since_epoch(datetime: datetime.datetime, days_since_julian_day_zero: float=2451545.0) -> float:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_days_since_reference_julian_day()` instead.
    
    
    Convert a calendar date to Julian days since a given epoch.
    
    
    Parameters
    ----------
    datetime : datetime.datetime
        Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored.
    days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
        Reference epoch (in days) since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0) corresponding to the 1st of January 2000.
    Returns
    -------
    float
        Date in Julian days since the given epoch.
    
    
    
    
    
    Examples
    --------
    In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days since J2000 (the 1st of January 2000).
    
    .. code-block:: python
    
      # Define the calendar date using datetime
      calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
      # Convert the calendar date to Julian days since J2000
      julian_date = time_representation.python_datetime_to_days_since_epoch(calendar_date)
      # Print the converted output
      print(julian_date)  # prints 8176.07825231459"""

def python_datetime_to_julian_day(datetime: datetime.datetime) -> float:
    """.. warning::
    
       This function is deprecated and will be removed in a future version of Tudat. Use `DateTime.from_python_datetime(...).to_julian_day()` instead.
    
    Convert a calendar date to Julian days.
    
    
    Parameters
    ----------
    datetime : datetime.datetime
        Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
    Returns
    -------
    float
        Julian day number (days since noon January 1st 4713 BC.)
    
    
    
    
    
    Examples
    --------
    In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days.
    
    .. code-block:: python
    
      # Define the calendar date using datetime
      calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
      # Convert the calendar date to Julian days since January 1st 4713 BC
      julian_date = time_representation.python_datetime_to_julian_day(calendar_date)
      # Print the converted output
      print(julian_date)  # prints 2459721.0782523146"""

def seconds_since_epoch_to_julian_centuries_since_epoch(seconds_since_epoch: float) -> float:
    """Convert the number of seconds since a given (unspecified) epoch to Julian centuries since the same epoch.
    
    Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch. This is equivalent to converting a time interval in seconds to Julian centuries
    
    Parameters
    ----------
    seconds_since_epoch : float
        Seconds elapsed since a given (unspecified) epoch.
    Returns
    -------
    float
        Julian centuries since the specified epoch.
    
        Since this is a float, not a integer, meaning that the fraction of the century is also included.
    
    
    
    
    
    
    Examples
    --------
    In this example, `706413165.12` seconds since a given epoch are converted to Julian centuries since the same epoch.
    
    .. code-block:: python
    
      # Define the number of seconds elapsed
      seconds_since_epoch = 706413165.12
      # Convert the number of seconds to Julian centuries
      julian_centuries = time_representation.seconds_since_epoch_to_julian_centuries_since_epoch(seconds_since_epoch)
      # Print the converted output
      print(julian_centuries)  # prints 0.2238488240930869"""

def seconds_since_epoch_to_julian_day(seconds_since_epoch: float, days_since_julian_day_zero: float=2451545.0) -> float:
    """Convert seconds since a given reference epoch to a Julian day.
    
    Inverse function of :func:`julian_day_to_seconds_since_epoch`.
    
    Parameters
    ----------
    seconds_since_epoch : float
        Seconds since ``days_since_julian_day_zero`` which are to be converted to date in Julian days.
    days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
        Reference epoch (in days since January 1st 4713 BC) since when the number of seconds have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
    Returns
    -------
    float
        Date in Julian days since January 1st 4713 BC, as computed from the input parameters
    
    
    
    
    
    Examples
    --------
    In this example, an amount of seconds since J2000 (January 1st 2000) is converted to the Julian date (in days since January 1st 4713 BC).
    
    .. code-block:: python
    
      # Define the amount of seconds since January 1st 2000
      seconds_since_J2000 = 706413165.1200145
      # Convert the amount of seconds since J2000 to the Julian date
      julian_date = time_representation.seconds_since_epoch_to_julian_day(seconds_since_J2000)
      # Print the converted output
      print(julian_date)  # prints 2459721.0783"""

def seconds_since_epoch_to_julian_years_since_epoch(seconds_since_epoch: float) -> float:
    """Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch.
    
    Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch. This is equivalent to converting a time interval in seconds to Julian years
    
    Parameters
    ----------
    seconds_since_epoch : float
        Seconds elapsed since a given (unspecified) epoch.
    Returns
    -------
    float
        Julian years since the specified epoch.
    
        Since this is a float, not a integer, meaning that the fraction of the year is also included.
    
    
    
    
    
    
    Examples
    --------
    In this example, `706413165.12` seconds since a given epoch are converted to Julian years since the same epoch.
    
    .. code-block:: python
    
      # Define the number of seconds elapsed
      seconds_since_epoch = 706413165.12
      # Convert the number of seconds to Julian years
      julian_years = time_representation.seconds_since_epoch_to_julian_years_since_epoch(seconds_since_epoch)
      # Print the converted output
      print(julian_years)  # prints 22.38488240930869"""

def year_and_days_in_year_to_calendar_date(year: int, days_in_year: int) -> DateTime:
    """.. warning::
    
        This function is deprecated and will be removed in a future version of Tudat. Use :func:`DateTime.from_year_and_day_of_year` instead.
    
    Create the Tudat :class:`DateTime` from the year and the number of days in the year.
    
    Parameters
    ----------
    year : int
        Calendar year.
    days_in_year : int
        Number of days that have passed in the year.
    Returns
    -------
    DateTime
        Corresponding calendar date as a :class:`DateTime` object. Note: the hours, minutes and seconds in the object are set to 0 when calling this function.
    
    Examples
    --------
    In this example, the calendar date corresponding to when 122 days have passed in 2020 is computed.
    
    .. code-block:: python
    
        # Compute the calendar date when 122 days have passed in 2020
        currentDate = time_representation.year_and_days_in_year_to_calendar_date(2020, 122)
        # Print the converted output
        print(currentDate)  # prints (2020, 5, 2, 0, 0)"""
tai_scale: TimeScales
tdb_scale: TimeScales
tt_scale: TimeScales
ut1_scale: TimeScales
utc_scale: TimeScales