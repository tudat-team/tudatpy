import numpy
import pybind11_stubgen.typing_ext
from ..math import interpolators
import typing
__all__ = ['OdfRawFileContents', 'SolarActivityContainer', 'SolarActivityData', 'TrackingDataType', 'TrackingTxtFileContents', 'day', 'doppler_bandwidth', 'doppler_base_frequency', 'doppler_measured_frequency', 'doppler_noise', 'downlink_frequency', 'dsn_receiving_station_nr', 'dsn_transmitting_station_nr', 'file_name', 'get_atmosphere_tables_path', 'get_earth_orientation_path', 'get_ephemeris_path', 'get_gravity_models_path', 'get_quadrature_path', 'get_resource_path', 'get_space_weather_path', 'get_spice_kernel_path', 'grail_antenna_file_reader', 'grail_mass_level_0_file_reader', 'grail_mass_level_1_file_reader', 'hour', 'light_time_measurement_accuracy', 'light_time_measurement_delay', 'minute', 'month', 'n_way_light_time', 'observation_body', 'observation_time_scale', 'observed_body', 'planet_nr', 'read_matrix_history_from_file', 'read_odf_file', 'read_solar_activity_data', 'read_tracking_txt_file', 'read_vector_history_from_file', 'receiving_station_name', 'residual_de405', 'second', 'set_dsn_weather_data_in_ground_stations', 'signal_to_noise', 'spacecraft_id', 'spacecraft_transponder_delay', 'spectral_max', 'tdb_spacecraft_j2000', 'time_tag_delay', 'transmitting_station_name', 'uplink_frequency', 'vx_planet_frame', 'vy_planet_frame', 'vz_planet_frame', 'x_planet_frame', 'y_planet_frame', 'year', 'z_planet_frame']

class OdfRawFileContents:
    """No documentation available."""

    def write_to_text_file(self, output_file: str) -> None:
        """
        No documentation available.
        """

class SolarActivityContainer:

    def __init__(self, solar_activity_data_map: dict[float, SolarActivityData]) -> None:
        ...

    def get_solar_activity_data(self, time: float) -> SolarActivityData:
        """
                Returns the nearest SolarActivityData (in UTC Julian days) for the given time in seconds since J2000.
        """

    def get_solar_activity_data_map(self) -> dict[float, SolarActivityData]:
        """
        Returns the full map of SolarActivityData.
        """

class SolarActivityData:
    """No documentation available."""

    @property
    def solar_radio_flux_107_observed(self) -> float:
        ...

class TrackingDataType:
    """No documentation available.
    
    Members:
    
      year : No documentation available.
    
      month : No documentation available.
    
      day : No documentation available.
    
      hour : No documentation available.
    
      minute : No documentation available.
    
      second : No documentation available.
    
      time_tag_delay : No documentation available.
    
      observation_time_scale : No documentation available.
    
      file_name : No documentation available.
    
      n_way_light_time : No documentation available.
    
      light_time_measurement_delay : No documentation available.
    
      light_time_measurement_accuracy : No documentation available.
    
      dsn_transmitting_station_nr : No documentation available.
    
      dsn_receiving_station_nr : No documentation available.
    
      observation_body : No documentation available.
    
      observed_body : No documentation available.
    
      spacecraft_id : No documentation available.
    
      planet_nr : No documentation available.
    
      tdb_spacecraft_j2000 : No documentation available.
    
      x_planet_frame : No documentation available.
    
      y_planet_frame : No documentation available.
    
      z_planet_frame : No documentation available.
    
      vx_planet_frame : No documentation available.
    
      vy_planet_frame : No documentation available.
    
      vz_planet_frame : No documentation available.
    
      residual_de405 : No documentation available.
    
      spacecraft_transponder_delay : No documentation available.
    
      uplink_frequency : No documentation available.
    
      downlink_frequency : No documentation available.
    
      signal_to_noise : No documentation available.
    
      spectral_max : No documentation available.
    
      doppler_measured_frequency : No documentation available.
    
      doppler_base_frequency : No documentation available.
    
      doppler_noise : No documentation available.
    
      doppler_bandwidth : No documentation available.
    
      receiving_station_name : No documentation available.
    
      transmitting_station_name : No documentation available."""
    __members__: typing.ClassVar[dict[str, TrackingDataType]]
    day: typing.ClassVar[TrackingDataType]
    doppler_bandwidth: typing.ClassVar[TrackingDataType]
    doppler_base_frequency: typing.ClassVar[TrackingDataType]
    doppler_measured_frequency: typing.ClassVar[TrackingDataType]
    doppler_noise: typing.ClassVar[TrackingDataType]
    downlink_frequency: typing.ClassVar[TrackingDataType]
    dsn_receiving_station_nr: typing.ClassVar[TrackingDataType]
    dsn_transmitting_station_nr: typing.ClassVar[TrackingDataType]
    file_name: typing.ClassVar[TrackingDataType]
    hour: typing.ClassVar[TrackingDataType]
    light_time_measurement_accuracy: typing.ClassVar[TrackingDataType]
    light_time_measurement_delay: typing.ClassVar[TrackingDataType]
    minute: typing.ClassVar[TrackingDataType]
    month: typing.ClassVar[TrackingDataType]
    n_way_light_time: typing.ClassVar[TrackingDataType]
    observation_body: typing.ClassVar[TrackingDataType]
    observation_time_scale: typing.ClassVar[TrackingDataType]
    observed_body: typing.ClassVar[TrackingDataType]
    planet_nr: typing.ClassVar[TrackingDataType]
    receiving_station_name: typing.ClassVar[TrackingDataType]
    residual_de405: typing.ClassVar[TrackingDataType]
    second: typing.ClassVar[TrackingDataType]
    signal_to_noise: typing.ClassVar[TrackingDataType]
    spacecraft_id: typing.ClassVar[TrackingDataType]
    spacecraft_transponder_delay: typing.ClassVar[TrackingDataType]
    spectral_max: typing.ClassVar[TrackingDataType]
    tdb_spacecraft_j2000: typing.ClassVar[TrackingDataType]
    time_tag_delay: typing.ClassVar[TrackingDataType]
    transmitting_station_name: typing.ClassVar[TrackingDataType]
    uplink_frequency: typing.ClassVar[TrackingDataType]
    vx_planet_frame: typing.ClassVar[TrackingDataType]
    vy_planet_frame: typing.ClassVar[TrackingDataType]
    vz_planet_frame: typing.ClassVar[TrackingDataType]
    x_planet_frame: typing.ClassVar[TrackingDataType]
    y_planet_frame: typing.ClassVar[TrackingDataType]
    year: typing.ClassVar[TrackingDataType]
    z_planet_frame: typing.ClassVar[TrackingDataType]

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

class TrackingTxtFileContents:
    """No documentation available."""

    def __init__(self, file_name: str, column_types: list[str], comment_symbol: str='#', value_separators: str=',:\t ') -> None:
        """
        No documentation available.
        """

    def add_metadata_str(self, tracking_data_type: TrackingDataType, str_value: str) -> None:
        """
        No documentation available.
        """

    def add_metadata_val(self, tracking_data_type: TrackingDataType, value: float) -> None:
        """
        No documentation available.
        """

    def get_available_datatypes(self) -> list[TrackingDataType]:
        """
        No documentation available.
        """

    @property
    def column_field_types(self) -> list[str]:
        """
        No documentation available.
        """

    @property
    def double_datamap(self) -> dict[TrackingDataType, list[float]]:
        """
        No documentation available.
        """

    @property
    def num_rows(self) -> int:
        """
        No documentation available.
        """

    @property
    def raw_datamap(self) -> dict[str, list[str]]:
        """
        No documentation available.
        """

def get_atmosphere_tables_path() -> str:
    """Get the path at which tudat atmosphere tables are located.
    
    Returns
    -------
    str
        Local path at which tudat atmosphere tables are located."""

def get_earth_orientation_path() -> str:
    """Get the path at which the Earth orientation resources used by tudat are located.
    
    Returns
    -------
    str
        Local path at which tudat Earth orientation resources are located."""

def get_ephemeris_path() -> str:
    """Get the path at which the ephemeris used by tudat are located.
    
    Returns
    -------
    str
        Local path at which the tudat ephemeris resources are located."""

def get_gravity_models_path() -> str:
    """Get the path at which tudat gravity models are located.
    
    Returns
    -------
    str
        Local path at which tudat gravity models are located."""

def get_quadrature_path() -> str:
    """Get the path at which the Gaussian quadrature resources are located.
    
    Returns
    -------
    str
        Local path at which tudat Gaussian quadrature resources are located."""

def get_resource_path() -> str:
    """Get the path at which tudat resources are located.
    
    Returns
    -------
    str
        Local path at which tudat resources are located."""

def get_space_weather_path() -> str:
    """Get the path at which tudat space weather is located.
    
    Returns
    -------
    str
        Local path at which tudat space weather is located."""

def get_spice_kernel_path() -> str:
    """Get the path at which the SPICE kernel used by tudat is located.
    
    Returns
    -------
    str
        Local path at which the SPICE kernel is located."""

def grail_antenna_file_reader(file_name: str) -> tuple[list[float], list[float]]:
    """No documentation available."""

def grail_mass_level_0_file_reader(file_name: str) -> dict[float, float]:
    """No documentation available."""

def grail_mass_level_1_file_reader(file_name: str, data_level: str='1b') -> dict[float, float]:
    """No documentation available."""

def read_matrix_history_from_file(matrix_rows: int, matrix_columns: int, file_name: str) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
    """Read a matrix history from a file.
    
    
    Parameters
    ----------
    matrix_rows : int
        Number of rows in the matrix at each epoch.
    matrix_columns : int
        Number of columns in the matrix at each epoch.
    file_name : str
        Name of the file containing the matrix history.
    Returns
    -------
    Dict[float, numpy.ndarray]
        Dictionary mapping epochs to the matrix at the given epoch."""

def read_odf_file(file_name: str) -> OdfRawFileContents:
    """No documentation available."""

def read_solar_activity_data(file_path: str) -> dict[float, SolarActivityData]:
    """Reads a space weather data file and produces a dictionary with solar activity data for a range of epochs. Data files can be obtained from http://celestrak.com/SpaceData and should follow the legacy format.
    
    :param file_path: Path to the space weather data file."""

def read_tracking_txt_file(file_name: str, column_types: list[str], comment_symbol: str='#', value_separators: str=',:\t ', ignore_omitted_columns: bool=False) -> TrackingTxtFileContents:
    ...

def read_vector_history_from_file(vector_size: int, file_name: str) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
    """Read a vector history from a file.
    
    
    Parameters
    ----------
    vector_size : int
        Size of the vector at each epoch.
    file_name : str
        Name of the file containing the vector history.
    Returns
    -------
    Dict[float, numpy.ndarray]
        Dictionary mapping epochs to the vector at the given epoch."""

def set_dsn_weather_data_in_ground_stations(bodies: ..., weather_file_names: list[str], interpolator_settings: interpolators.InterpolatorSettings=..., ground_stations_per_complex: dict[int, list[str]]={10: ['DSS-13', 'DSS-14', 'DSS-15', 'DSS-24', 'DSS-25', 'DSS-26', 'DSS-27'], 40: ['DSS-34', 'DSS-35', 'DSS-36', 'DSS-43', 'DSS-45'], 60: ['DSS-54', 'DSS-55', 'DSS-63', 'DSS-65']}, body_with_ground_stations_name: str='Earth') -> None:
    """No documentation available."""
day: TrackingDataType
doppler_bandwidth: TrackingDataType
doppler_base_frequency: TrackingDataType
doppler_measured_frequency: TrackingDataType
doppler_noise: TrackingDataType
downlink_frequency: TrackingDataType
dsn_receiving_station_nr: TrackingDataType
dsn_transmitting_station_nr: TrackingDataType
file_name: TrackingDataType
hour: TrackingDataType
light_time_measurement_accuracy: TrackingDataType
light_time_measurement_delay: TrackingDataType
minute: TrackingDataType
month: TrackingDataType
n_way_light_time: TrackingDataType
observation_body: TrackingDataType
observation_time_scale: TrackingDataType
observed_body: TrackingDataType
planet_nr: TrackingDataType
receiving_station_name: TrackingDataType
residual_de405: TrackingDataType
second: TrackingDataType
signal_to_noise: TrackingDataType
spacecraft_id: TrackingDataType
spacecraft_transponder_delay: TrackingDataType
spectral_max: TrackingDataType
tdb_spacecraft_j2000: TrackingDataType
time_tag_delay: TrackingDataType
transmitting_station_name: TrackingDataType
uplink_frequency: TrackingDataType
vx_planet_frame: TrackingDataType
vy_planet_frame: TrackingDataType
vz_planet_frame: TrackingDataType
x_planet_frame: TrackingDataType
y_planet_frame: TrackingDataType
year: TrackingDataType
z_planet_frame: TrackingDataType