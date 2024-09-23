from astropy.table.operations import vstack
import astropy.table.table
from astropy.table.table import Table
import astropy.time.core
from astropy.time.core import Time
from astropy.time.core import TimeDelta
from astropy import units as u
import astroquery.jplhorizons.core
import datetime as datetime
import math as math
import numpy
import numpy as np
import pandas as pd
import re as re
from tudatpy import constants
import tudatpy.numerical_simulation.environment_setup.ephemeris.expose_ephemeris
from tudatpy.numerical_simulation.environment_setup.ephemeris.expose_ephemeris import EphemerisSettings
from tudatpy.numerical_simulation.environment_setup.ephemeris.expose_ephemeris import tabulated
import tudatpy.numerical_simulation.environment_setup.expose_environment_setup
from tudatpy.numerical_simulation.environment_setup.expose_environment_setup import BodyListSettings
import typing
__all__ = ['BodyListSettings', 'EphemerisSettings', 'Horizons', 'HorizonsBatch', 'HorizonsQuery', 'Table', 'Time', 'TimeDelta', 'constants', 'datetime', 'jpl_horizons', 'math', 'np', 'pd', 're', 'tabulated', 'u', 'vstack']

class HorizonsBatch:

    def __init__(self, query_id_list: list[str], location: str, epoch_list: typing.Optional[list]=None, epoch_start: typing.Union[datetime.datetime, float, NoneType]=None, epoch_end: typing.Union[datetime.datetime, float, NoneType]=None, epoch_step: typing.Optional[str]=None, extended_query: bool=False) -> None:
        """
        Query object to retrieve Horizons data for multiple bodies.
                See documentation for `HorizonsQuery` for more extensive documentation.
                Note that the epochs requested must have data for all bodies queried.
                This class is useful for quickly creating ephemerides for many objects at once.
                For general purposes, use `HorizonsQuery instead`.
        
                Parameters
                ----------
                query_id : str
                    List of query terms to retrieve data for,
                    all queries behave as type default and the `query_type` behaviour
                    can not be set for the `HorizonsBatch` class.
                location : str
                    Coordinate centre for the data with syntax `site@body`.
                epoch_start : Union[datetime.datetime, float, None], optional
                    Starting date to retrieve data for.
                    Must be either a python datetime object or a float seconds since J2000 TDB.
                    Combined with `epoch_end` and
                    `epoch_step` can be used to retrieve a range of data.
                    If `epoch_list` is used, value must be None, by default None
                epoch_end : Union[datetime.datetime, float, None], optional
                    Final date to retrieve data for.
                    Must be either a python datetime object or a float seconds since J2000 TDB.
                    If the start stop, and step parameters dont result in an
                    integer multiple of values,the final date is not used.
                    Combined with `epoch_end` and
                    `epoch_step` can be used to retrieve a range of data.
                    If `epoch_list` is used, value must be None, by default None
                epoch_step : Union[str, None], optional
                    Step for the range of epochs to retrieve data for.
                    Can be be either a specific timestep such as `15m`,
                    or a number of partitions `10` for 10 equidistant steps within the range.
                    For the timestep, quantifiers `m` for minute, `h` for hour and `d` for day
                    can be used. Seconds are not available, consider using the `epoch_list`
                    parameter or a number of partitions instead.
                    Month and Year are normally available in JPL Horizons but are restricted
                    here because they may produce ambiguous results (leap years etc.).
                    If `epoch_list` is used, value must be None, by default None
                epoch_list : Union[list, None], optional
                    List of times in seconds since J2000 TDB. Can be used
                    to retrieve specific times instead of a range.
                    Must be None if start, end and step are set, by default None
                extended_query : bool, optional
                    Enables the retrieval of larger collections of data, by default False.
                
        """

    def add_batch_ephemerides(self, body_settings: tudatpy.numerical_simulation.environment_setup.expose_environment_setup.BodyListSettings, frame_origin: str, frame_orientation: str='ECLIPJ2000', aberations: str='geometric') -> None:
        """
        Uses the data queried to add ephemerides of the bodies querried
                to the body_settings. The names of the bodies added can be retrieved
                using the names property.
        
                Parameters
                ----------
                body_settings : environment_setup.BodyListSettings
                    Tudat body settings object.
                frame_origin : str
                    Global frame origin, should match the queries' location parameter.
                frame_orientation : str, optional
                    Reference Frame Orientation, equivalent to the astroquery refplane
                    parameter. Options are 'J2000' and 'ECLIPJ2000', by default "ECLIPJ2000".
                aberations : str, optional
                    Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
                    'apparent', by default "geometric".
        
                See the Horizons System Manual for more info:
        
                https://ssd.jpl.nasa.gov/horizons/manual.html#output
        
                
        """

    @property
    def names(self) -> typing.Optional[list[str]]:
        """
        Retrieves a list of names of the query objects. Returns `None`
                if `add_batch_ephemerides` has not been run yet.
        """

class HorizonsQuery:
    """This class provides an interface to JPL's Horizon System.
	JPL Horizons provides access to highly accurate ephemerides for many solar system objects,
	including asteriuds, comets, planets, moons and select spacecraft.
	The class extends astroquery's to cater to the needs of Tudat users,
	while maintaining compatibility with all of astroquery's features.
	
	There are some notable differences:
	
	Time input has been simplified to reduce ambiguities:
	- List of times are given in seconds since J2000 TDB.
	- Start can be given in datetime format or seconds since J2000 TDB.
	- Timesteps like months and years are not permitted.
	
	And some additional features:
	- Extended query allows data retrieval limits
	to be broken by automatically splitting up a query into multiple subqueries
	and combining the data.
	- Ephemeris settings can automatically be generated using Vectors API.
	
	
	
	"""
    query_limit: typing.ClassVar[int] = 90024
    query_limit_list: typing.ClassVar[int] = 50

    def __init__(self, query_id: str, location: str, query_type: str='default', epoch_start: typing.Union[datetime.datetime, float, NoneType]=None, epoch_end: typing.Union[datetime.datetime, float, NoneType]=None, epoch_step: typing.Optional[str]=None, epoch_list: typing.Optional[list]=None, extended_query: bool=False) -> None:
        """
        Query object to retrieve Horizons data for a singular body.
        
                Parameters
                ----------
                query_id : str
                    Query term for object to retrieve data for. The query may be ambiguous,
                    in which case suggestions for exact objects are presented.
                    The behaviour for this parameter changes based on the `query_type`
                    parameter.
                    The default behaviour is the same as the official Horizons API.
        
                    Here are some examples for the default behaviour:
        
                    - `Earth` - ambiguous, will
                    suggest 3 (Earth-Moon Barycentre) and 399 (Earth).
        
                    - `3` - will retrieve Earth-Moon Barycentre.
        
                    - `3;` - semi-colon searches for minor planets.            In this case it will search for the minor planet with MPC code 3: Juno.
        
                    - `-3` - A minus sign searches for spacecraft.
                    In this case the Mars Orbiter Mission.
        
                    See the Horizons System manual
                    for an extensive explanation on this parameter and the timespans
                    page for some examples:
                    https://ssd.jpl.nasa.gov/horizons/manual.html#select
        
                    https://ssd.jpl.nasa.gov/horizons/time_spans.html
        
                location : str
                    Coordinate centre for the data with syntax `site@body`.
        
                    In general, the syntax is: `site@body` but there are several shorthands.
                    The sight may be a specific location on the body, such as an observatory.
                    500 specifies the geocentre, and can be left out as a shorthand
                    (500@399 = @399).
                    Here are some examples:
        
                    - `@10` or `@Sun` - if no site is given, the geocenter will be taken.
        
                    - `500@399` or `@399` or `500` - Geocentric. The site defaults to `500`.
                    The body defaults to `@399`.
        
                    - `@0` or `@SSB` - The solar system barycentre.
        
                    - `0` - without the `@` symbol, this location is equivalent to `0@399`
                    which is the observatory with MPC code 0 on Earth and not the SSB.
                    In this case: Greenwich Observatory.
        
                    - `Greenwich` - Equivalent to `Greenwich@399`, Greenwich Observatory.
        
                    - `Earth` - Heaven on Earth Observatory, not equivalent to `@399`.
        
                    - `@Earth` - ambiguous, will suggest `@399` or `@3`.
        
                    - `Geocentric` - Equivalent to `500@399`.
        
                    As mistakes can easily made, it is highly recommended to consult
                    the manual and use specific codes for this parameter:
        
                    https://ssd.jpl.nasa.gov/horizons/manual.html#center
                query_type : str, optional
                    The query type constrains the search for the `query_id` parameter.
                    Can have values: 'default', 'None', 'smallbody', 'designation', 'name',
                    'asteroid_name', 'comet_name', by default "default"
        
                    While all objects can be found with the default parameter,
                    some queries have ambiguities. Specifying the query type eliminates the
                    specific syntax. For example searching Juno with `query_type` `spacecraft`
                    retrieves data for the spacecraft while setting `smallbody`
                    retrieves data for the asteroid Juno.
                    See the Horizons system Manual and astroquery documentation for more info:
        
                    https://ssd.jpl.nasa.gov/horizons/manual.html#select
        
                    https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html
        
                epoch_start : Union[datetime.datetime, float, None], optional
                    Starting date to retrieve data for.
                    Must be either a python datetime object or a float seconds since J2000 TDB.
                    Combined with `epoch_end` and
                    `epoch_step` can be used to retrieve a range of data.
                    If `epoch_list` is used, value must be None, by default None
                epoch_end : Union[datetime.datetime, float, None], optional
                    Final date to retrieve data for.
                    Must be either a python datetime object or a float seconds since J2000 TDB.
                    If the start stop, and step parameters dont result in an
                    integer multiple of values,the final date is not used.
                    Combined with `epoch_end` and
                    `epoch_step` can be used to retrieve a range of data.
                    If `epoch_list` is used, value must be None, by default None
                epoch_step : Union[str, None], optional
                    Step for the range of epochs to retrieve data for.
                    Can be be either a specific timestep such as `15m`,
                    or a number of partitions `10` for 10 equidistant steps within the range.
                    For the timestep, quantifiers `m` for minute, `h` for hour and `d` for day
                    can be used. Seconds are not available, consider using the `epoch_list`
                    parameter or a number of partitions instead.
                    Month and Year are normally available in JPL Horizons but are restricted
                    here because they may produce ambiguous results (leap years etc.).
                    If `epoch_list` is used, value must be None, by default None
                epoch_list : Union[list, None], optional
                    List of times in seconds since J2000 TDB. Can be used
                    to retrieve specific times instead of a range.
                    Must be None if start, end and step are set, by default None
                extended_query : bool, optional
                    Enables the retrieval of larger collections of data, by default False.
                    Horizons System has a limit on how much data can be output (90024) lines.
                    Additionally there is a limit on input length which may be exceeded if a
                    long list of times is given with the `epoch_list` option (50 epochs).
                    When enabled, the times are split in to multiple subqueries.
                    The data of these subqueries is then combined and presented as a single
                    query. Using a very large number of queries may make data retrieval slow,
                    especially when using the list option. Consider whether your use case
                    requires a large number of specific times,
                    or if a range can be used instead.
        
                    The epoch step partitions (`10` not `10h`) is currently unsupported with
                    the extended_query.
        
                
        """

    def _convert_time_to_astropy(self, time) -> astropy.time.core.Time:
        """
        Internal method to convert inputted times to astropy
        """

    def _format_time_list(self, times) -> list[astropy.time.core.Time]:
        """
        Internal method to format input time lists to a JPL accepted format
        """

    def _format_time_range(self, time: typing.Union[float, datetime.datetime]) -> str:
        """
        Internal method to format input time ranges to a JPL accepted format
        """

    def _infer_name(self) -> None:
        """
        Internal Method to extract the name, designation and MPC code for the
                object from the data.
        """

    def _interpret_timestep(self, timestep: str, start: typing.Union[datetime.datetime, float], end: typing.Union[datetime.datetime, float]):
        """
        Internal method to determine the length of the timestep in seconds
                This is to help partition the queries in the case of an extended query
        """

    def _validate_time_range(self, start, end) -> None:
        """
        Internal method to validate start/end user input
        """

    def cartesian(self, frame_orientation: str='ECLIPJ2000', aberations: str='geometric') -> numpy.ndarray:
        """
        Retrieve the cartesian state using the Horizons Vector API.
        
                Parameters
                ----------
                frame_orientation : str, optional
                    Reference Frame Orientation, equivalent to the astroquery refplane
                    parameter. Options are 'J2000'/'earth', 'ECLIPJ2000'/'ecliptic' and 'body'
                    , by default "ECLIPJ2000".
                    See astroquery documentation for information about the 'body' option:
        
                    https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html
                aberations : str, optional
                    Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
                    'apparent', by default "geometric".
        
                    See the Horizons System Manual for more info:
        
                    https://ssd.jpl.nasa.gov/horizons/manual.html#output
        
                Returns
                -------
                np.ndarray
                    returns an n by 7 array with the time in seconds since J2000 TDB,
                    and the cartesian position and velocities.
                
        """

    def create_ephemeris_tabulated(self, frame_origin: str, frame_orientation: str='ECLIPJ2000', aberations: str='geometric') -> tudatpy.numerical_simulation.environment_setup.ephemeris.expose_ephemeris.EphemerisSettings:
        """
        Create ephemeris settings for a body using Horizons Vector API.
        
                Parameters
                ----------
                frame_origin : str
                    Global frame origin, should match the queries' location parameter.
                frame_orientation : str, optional
                    Reference Frame Orientation, equivalent to the astroquery refplane
                    parameter. Options are 'J2000' and 'ECLIPJ2000', by default "ECLIPJ2000".
                aberations : str, optional
                    Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
                    'apparent', by default "geometric".
        
                    See the Horizons System Manual for more info:
        
                    https://ssd.jpl.nasa.gov/horizons/manual.html#output
        
                Returns
                -------
                environment_setup.ephemeris.EphemerisSettings
                    Ephemeris settings for the query's body.
        
                Examples
                ----------
                Add Ephemerides of JUICE to the body_settings
        
                >>> body_settings.add_empty_settings("JUICE")
                >>> body_settings.get("JUICE").ephemeris_settings = query.create_ephemeris_tabulated(
                        frame_origin=global_frame_origin,
                        frame_orientation=global_frame_orientation,
                    )
                
        """

    def ephemerides(self, reference_system: str='J2000', extra_precision: bool=False, *args, **kwargs) -> astropy.table.table.Table:
        """
        Implements the JPL Horizons ephemerides API and returns it in raw Astropy table format.
                Ephemerides API provides time-interpolated observer parameters such as right ascension and declination.
                Note that this means that values provided are not actual observations.
        
                A number of quantities are retrieved, their definitions can be found here:
                https://ssd.jpl.nasa.gov/horizons/manual.html#obsquan.
                By default all available quantities are retrieved.
        
                More parameters can be passed directly to the astroquery call. These can be passed as kwargs: kwargs=("refraction":True).
                Check the astroquery documentation for an overview:
                https://astroquery.readthedocs.io/en/latest/api/astroquery.jplhorizons.HorizonsClass.html#astroquery.jplhorizons.HorizonsClass.ephemerides
        
        
                Parameters
                ----------
                reference_system : str, optional
                    Coordinate reference system, value must be one of `ICRF`/`J2000` or B1950, by default "J2000"
                extra_precision : bool, optional
                    Enables extra precision in right ascension and declination values, by default False
        
                Returns
                -------
                astropy.table.Table
                    Unprocessed output in astropy table format.
        
                Raises
                ------
                ValueError
                    If time query has incorrect format or an incorrect reference system is chosen
                
        """

    def interpolated_observations(self, degrees: bool=False, reference_system: str='J2000', extra_precision: bool=True, *args, **kwargs) -> numpy.ndarray:
        """
        Retrieves interpolated Right Ascension and Declination from the Horizons ephemerides API.
                Note that these values are not real observations but instead interpolated
                values based on the Horizons ephemeris system.
        
                Parameters
                ----------
                degrees : bool, optional
                    return values in degrees if True, radians if False, by default false
                reference_system : str, optional
                    Coordinate reference system, value must be one of `ICRF`/`J2000` or B1950, by default "J2000"
                extra_precision : bool, optional
                    Enables extra precision in Right Ascension and Declination values, by default False
        
                Returns
                -------
                np.ndarray
                    Numpy array (N, 3) with time in seconds since J2000 TDB and the Right Ascension and Declination.
        
                Raises
                ------
                ValueError
                    If time query has incorrect format or an incorrect reference system is chosen
                
        """

    def vectors(self, frame_orientation: str='ECLIPJ2000', aberations: str='geometric') -> astropy.table.table.Table:
        """
        Retrieve Horizons Vectors api data in raw astropy format.
                For general purposes, use the `.cartesian()` method instead.
        
                Parameters
                ----------
                frame_orientation : str, optional
                    Reference Frame Orientation, equivalent to the astroquery refplane
                    parameter. Options are 'J2000'/'earth', 'ECLIPJ2000'/'ecliptic' and 'body'
                    , by default "ECLIPJ2000".
                    See astroquery documentation for information about the 'body' option:
        
                    https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html
                aberations : str, optional
                    Aberations to be accounted for. Options are: 'geometric', 'astrometric' and
                    'apparent', by default "geometric".
        
                    See the Horizons System Manual for more info:
        
                    https://ssd.jpl.nasa.gov/horizons/manual.html#output
        
                Returns
                -------
                astropy.table.Table
                    Unprocessed vectors API data in astropy Table format.
                
        """

    @property
    def MPC_number(self) -> typing.Optional[str]:
        """
        Retrieve the MPC (Minor Planet Centre) number of the object.
                The MPC number is infered from data retrieved and will return none
                if data has not been retrieved yet.
                The MPC number is only relevant to minor planets such as asteroids, TNOs and
                Near-Earth Asteroids.
        """

    @property
    def designation(self) -> typing.Optional[str]:
        """
        Retrieve the relevant designation of the query's object.
                The designation is infered from the data retrieved and will return none
                if data has not been retrieved yet.
                Minor planets and Comets will return their provisional designation
                (1898 DQ for Eros, 1982 HG1 for Halley).
                A comets' formal designation can often be retrieved using the `name` property
                Spacecraft and Major Planets/ Moons will return their JPL number
                (-28 for JUICE, 6 for Saturn Barycentre).
        """

    @property
    def name(self) -> typing.Optional[str]:
        """
        Retrieve the name of the query's object.
                The name is infered from the data retrieved and will return none
                if data has not been retrieved yet.
                Unnamed minor planets will use their designation instead.
                If a name can not be infered, the raw name from Horizons will be returned.
                Please consider raising an issue on the Tudat github in such cases.
        """

def jpl_horizons(horizons_query: str, horizons_location: str, frame_origin: str, frame_orientation: str='ECLIPJ2000', query_type: str='default', epoch_start: typing.Union[datetime.datetime, float, NoneType]=None, epoch_end: typing.Union[datetime.datetime, float, NoneType]=None, epoch_step: typing.Optional[str]=None, epoch_list: typing.Optional[list]=None, extended_query: bool=False, aberations: str='geometric'):
    """Factory function for creating ephemeris model settings from tabulated JPL Horizons vectors.
	
	JPL Horizons provides access to highly accurate ephemerides for many solar system objects,
	including asteriuds, comets, planets, moons and select spacecraft.
	
	This function is a wrapper for the tudatpy.data.horizons functionality.
	That api is not available on the api documentation yet.
	For now, visit the HorizonsQuery souce code for extensive documentation:
	https://github.com/tudat-team/tudatpy/blob/master/tudatpy/data/horizons.py
	
	For more information on the Horizons System, visit: https://ssd.jpl.nasa.gov/horizons/manual.html
	
	Examples
	----------
	Add Ephemerides of JUICE to the body_settings
	
		>>> juice_eph_settings = jpl_horizons(
				horizons_query="-121", #-121 is the query code for JUICE
				horizons_location="500@399", # Geocentre@Earth
				frame_origin="Earth", #tudat frame origin and orientation
				frame_orientation="ECLIPJ2000",
				epoch_start=datetime.datetime(2018, 10, 21),
				epoch_end=datetime.datetime(2023, 9, 1),
				epoch_step="1d",
				extended_query=True,
			)
	
		>>> body_settings.add_empty_settings("JUICE")
		>>> body_settings.get("JUICE").ephemeris_settings = juice_eph_settings
	
	"""
Horizons: astroquery.jplhorizons.core.HorizonsClass