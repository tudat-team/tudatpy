import tudatpy.math.interpolators.expose_interpolators
import tudatpy.numerical_simulation.environment.expose_environment
import typing
__all__ = ['DynamicCoefficientNames', 'OdfRawFileContents', 'SolarActivityData', 'StaticCoefficientNames', 'TrackingDataType', 'TrackingTxtFileContents', 'ca', 'cap', 'caq', 'car', 'clb', 'cll', 'cllp', 'cllq', 'cllr', 'cln', 'clnp', 'clnq', 'clnr', 'cm', 'cma', 'cmad', 'cmp', 'cmq', 'cmr', 'cn', 'cna', 'cnad', 'cnb', 'cnp', 'cnq', 'cnr', 'cy', 'cyb', 'cyp', 'cyq', 'cyr', 'day', 'doppler_bandwidth', 'doppler_base_frequency', 'doppler_measured_frequency', 'doppler_noise', 'downlink_frequency', 'dsn_receiving_station_nr', 'dsn_transmitting_station_nr', 'file_name', 'get_atmosphere_tables_path', 'get_earth_orientation_path', 'get_ephemeris_path', 'get_gravity_models_path', 'get_quadrature_path', 'get_resource_path', 'get_space_weather_path', 'get_spice_kernel_path', 'hour', 'light_time_measurement_accuracy', 'light_time_measurement_delay', 'minute', 'missile_DATCOM_data', 'month', 'n_way_light_time', 'observation_body', 'observation_time_scale', 'observed_body', 'planet_nr', 'read_matrix_history_from_file', 'read_odf_file', 'read_solar_activity_data', 'read_tracking_txt_file', 'read_vector_history_from_file', 'residual_de405', 'second', 'set_dsn_weather_data_in_ground_stations', 'signal_to_noise', 'spacecraft_id', 'spacecraft_transponder_delay', 'spectral_max', 'tdb_spacecraft_j2000', 'tdb_time_j2000', 'time_tag_delay', 'uplink_frequency', 'vlbi_station_name', 'vx_planet_frame', 'vy_planet_frame', 'vz_planet_frame', 'x_planet_frame', 'y_planet_frame', 'year', 'z_planet_frame']

class DynamicCoefficientNames:
    """Enumeration of Missile DATCOM dynamic aerodynamic coefficient types.
	
	
	:member cnq:
	:member cmq:
	:member caq:
	:member cyq:
	:member clnq:
	:member cllq:
	:member cnr:
	:member cmr:
	:member car:
	:member cyr:
	:member clnr:
	:member cllr:
	:member cnp:
	:member cmp:
	:member cap:
	:member cyp:
	:member clnp:
	:member cllp:
	:member cnad:
	:member cmad:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, DynamicCoefficientNames]]
    cap: typing.ClassVar[DynamicCoefficientNames]
    caq: typing.ClassVar[DynamicCoefficientNames]
    car: typing.ClassVar[DynamicCoefficientNames]
    cllp: typing.ClassVar[DynamicCoefficientNames]
    cllq: typing.ClassVar[DynamicCoefficientNames]
    cllr: typing.ClassVar[DynamicCoefficientNames]
    clnp: typing.ClassVar[DynamicCoefficientNames]
    clnq: typing.ClassVar[DynamicCoefficientNames]
    clnr: typing.ClassVar[DynamicCoefficientNames]
    cmad: typing.ClassVar[DynamicCoefficientNames]
    cmp: typing.ClassVar[DynamicCoefficientNames]
    cmq: typing.ClassVar[DynamicCoefficientNames]
    cmr: typing.ClassVar[DynamicCoefficientNames]
    cnad: typing.ClassVar[DynamicCoefficientNames]
    cnp: typing.ClassVar[DynamicCoefficientNames]
    cnq: typing.ClassVar[DynamicCoefficientNames]
    cnr: typing.ClassVar[DynamicCoefficientNames]
    cyp: typing.ClassVar[DynamicCoefficientNames]
    cyq: typing.ClassVar[DynamicCoefficientNames]
    cyr: typing.ClassVar[DynamicCoefficientNames]

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

class OdfRawFileContents:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def write_to_text_file(self, output_file: str) -> None:
        ...

class SolarActivityData:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class StaticCoefficientNames:
    """Enumeration of Missile DATCOM static aerodynamic coefficient types.
	
	
	:member cn:
	:member cm:
	:member ca:
	:member cy:
	:member cln:
	:member cll:
	:member cna:
	:member cma:
	:member cyb:
	:member cnb:
	:member clb:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, StaticCoefficientNames]]
    ca: typing.ClassVar[StaticCoefficientNames]
    clb: typing.ClassVar[StaticCoefficientNames]
    cll: typing.ClassVar[StaticCoefficientNames]
    cln: typing.ClassVar[StaticCoefficientNames]
    cm: typing.ClassVar[StaticCoefficientNames]
    cma: typing.ClassVar[StaticCoefficientNames]
    cn: typing.ClassVar[StaticCoefficientNames]
    cna: typing.ClassVar[StaticCoefficientNames]
    cnb: typing.ClassVar[StaticCoefficientNames]
    cy: typing.ClassVar[StaticCoefficientNames]
    cyb: typing.ClassVar[StaticCoefficientNames]

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

class TrackingDataType:
    """Members:
	
	year :
	
	month :
	
	day :
	
	hour :
	
	minute :
	
	second :
	
	time_tag_delay :
	
	observation_time_scale :
	
	file_name :
	
	n_way_light_time :
	
	light_time_measurement_delay :
	
	light_time_measurement_accuracy :
	
	dsn_transmitting_station_nr :
	
	dsn_receiving_station_nr :
	
	observation_body :
	
	observed_body :
	
	spacecraft_id :
	
	planet_nr :
	
	tdb_time_j2000 :
	
	tdb_spacecraft_j2000 :
	
	x_planet_frame :
	
	y_planet_frame :
	
	z_planet_frame :
	
	vx_planet_frame :
	
	vy_planet_frame :
	
	vz_planet_frame :
	
	residual_de405 :
	
	spacecraft_transponder_delay :
	
	uplink_frequency :
	
	downlink_frequency :
	
	signal_to_noise :
	
	spectral_max :
	
	doppler_measured_frequency :
	
	doppler_base_frequency :
	
	doppler_noise :
	
	doppler_bandwidth :
	
	vlbi_station_name :
	"""
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
    residual_de405: typing.ClassVar[TrackingDataType]
    second: typing.ClassVar[TrackingDataType]
    signal_to_noise: typing.ClassVar[TrackingDataType]
    spacecraft_id: typing.ClassVar[TrackingDataType]
    spacecraft_transponder_delay: typing.ClassVar[TrackingDataType]
    spectral_max: typing.ClassVar[TrackingDataType]
    tdb_spacecraft_j2000: typing.ClassVar[TrackingDataType]
    tdb_time_j2000: typing.ClassVar[TrackingDataType]
    time_tag_delay: typing.ClassVar[TrackingDataType]
    uplink_frequency: typing.ClassVar[TrackingDataType]
    vlbi_station_name: typing.ClassVar[TrackingDataType]
    vx_planet_frame: typing.ClassVar[TrackingDataType]
    vy_planet_frame: typing.ClassVar[TrackingDataType]
    vz_planet_frame: typing.ClassVar[TrackingDataType]
    x_planet_frame: typing.ClassVar[TrackingDataType]
    y_planet_frame: typing.ClassVar[TrackingDataType]
    year: typing.ClassVar[TrackingDataType]
    z_planet_frame: typing.ClassVar[TrackingDataType]

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

class TrackingTxtFileContents:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, file_name: str, column_types: list[str], comment_symbol: str='#', value_separators: str=',:\t ') -> None:
        ...

    @property
    def column_field_types(self) -> list[str]:
        ...

    @property
    def double_datamap(self) -> dict[..., list[float]]:
        ...

    @property
    def num_rows(self) -> int:
        ...

    @property
    def raw_datamap(self) -> dict[str, list[str]]:
        ...

class missile_DATCOM_data:
    """Class containing data and methods interfacing the Missile DATCOM software.
	
	This class is the main method that can be used to interface tudat with the Missile DATCOM software.
	It can be initialised with the output file from Missile DATCOM, and provides methods to convert these results
	into tudat-compatible data.
	
	.. note:: The Missile DATCOM software from which outputs can be interfaced to TUDAT is an entirely separate software from Tudat(Py).
			  Please refer to Missile DATCOM user manuals for information on how to use it. These can be accessed on the US Defence Technical
			  Information Center at accession numbers `ADA267447 <https://apps.dtic.mil/sti/citations/ADA267447>`_ and
			  `ADA503576 <https://apps.dtic.mil/sti/citations/ADA503576>`_.
	
	.. note:: The interfacing of Missile DATCOM to tudat assumes that aerodynamic coefficients are computed as a function of both
			  Mach number and angle of attack.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, file_name_and_path: str) -> None:
        """
        Class constructor.
        
        	Function used to construct and initialise the class. In essence, it can be used to read and extract the aerodynamic coefficients
        	computed by Missile DATCOM, and save them in different formats.
        
        
        	:param file_name_and_path:
        		Full path and file name of the `for004.dat` Missile DATCOM results output file.
        """

    def get_Reynolds_numbers(self) -> list[float]:
        """
        Get the list of Reynolds numbers at which Missile DATCOM has been run.
        	:return:
        		List of Reynolds numbers.
        """

    def get_angle_of_attacks(self) -> list[float]:
        """
        Get the list of angle of attacks at which Missile DATCOM has been run.
        	:return:
        		List of angle of attacks.
        """

    def get_dynamic_coefficient(self, mach_index: int, angle_of_attack_index: int, coefficient_index: ...) -> float:
        """
        Get a specific dynamic coefficient from the result database.
        
        	:param mach_index:
        		Index of the Mach number for which to get the static coefficient.
        	:param angle_of_attack_index:
        		Index of the angle of attack for which to get the static coefficient.
        	:param coefficient_index:
        		Type of the dynamic aerodynamic coefficient.
        	:return:
        		Dynamic aerodynamic coefficient.
        """

    def get_mach_numbers(self) -> list[float]:
        """
        Get the list of Mach numbers at which Missile DATCOM has been run.
        	:return:
        		List of Mach numbers.
        """

    def get_static_coefficient(self, mach_index: int, angle_of_attack_index: int, coefficient_index: ...) -> float:
        """
        Get a specific static coefficient from the result database.
        
        	:param mach_index:
        		Index of the Mach number for which to get the static coefficient.
        	:param angle_of_attack_index:
        		Index of the angle of attack for which to get the static coefficient.
        	:param coefficient_index:
        		Type of the static aerodynamic coefficient.
        	:return:
        		Static aerodynamic coefficient.
        """

    def write_all_coefficients_to_files(self, file_name_base: str, base_precision: int=15, exponent_width: int=2) -> None:
        """
        Write all the aerodynamic coefficients to CSV files.
        
        	:param file_name_base:
        		Full base path and name of the file that will be saved. The name of each aerodynamic coefficient will be included at the end of the file name.
        	:param base_precision:
        		Number of digits to represent the base of the floating-point number.
        	:param exponent_width:
        		Number of digits to represent the exponent of the floating-point number.
        """

    def write_force_and_moment_coefficients_to_files(self, file_name_base: str, base_precision: int=15, exponent_width: int=2) -> None:
        """
        Write the force and moment coefficients to a file in the format taken by the :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_from_files` function.
        
        	:param file_name_base:
        		Full base path and name of the file that will be saved. The name of each aerodynamic coefficient will be included at the end of the file name.
        	:param base_precision:
        		Number of digits to represent the base of the floating-point number.
        	:param exponent_width:
        		Number of digits to represent the exponent of the floating-point number.
        """

def get_atmosphere_tables_path() -> str:
    """Get the path at which tudat atmosphere tables are located.
	:return:
			Local path at which tudat atmosphere tables are located.
	"""

def get_earth_orientation_path() -> str:
    """Get the path at which the Earth orientation resources used by tudat are located.
	:return:
			Local path at which tudat Earth orientation resources are located.
	"""

def get_ephemeris_path() -> str:
    """Get the path at which the ephemeris used by tudat are located.
	:return:
			Local path at which the tudat ephemeris resources are located.
	"""

def get_gravity_models_path() -> str:
    """Get the path at which tudat gravity models are located.
	:return:
			Local path at which tudat gravity models are located.
	"""

def get_quadrature_path() -> str:
    """Get the path at which the Gaussian quadrature resources are located.
	:return:
			Local path at which tudat Gaussian quadrature resources are located.
	"""

def get_resource_path() -> str:
    """Get the path at which tudat resources are located.
	:return:
		Local path at which tudat resources are located.
	"""

def get_space_weather_path() -> str:
    """Get the path at which tudat space weather is located.
	:return:
			Local path at which tudat space weather is located.
	"""

def get_spice_kernel_path() -> str:
    """Get the path at which the SPICE kernel used by tudat is located.
	:return:
			Local path at which the SPICE kernel is located.
	"""

def read_matrix_history_from_file(matrix_rows: int, matrix_columns: int, file_name: str) -> dict[float, ..., -1, -1, 0, -1, ...]:
    """Read a matrix history from a file.
	
	:param matrix_rows:
			Number of rows in the matrix at each epoch.
	:param matrix_columns:
			Number of columns in the matrix at each epoch.
	:param file_name:
			Name of the file containing the matrix history.
	:return:
			Dictionary mapping epochs to the matrix at the given epoch.
	"""

def read_odf_file(file_name: str) -> OdfRawFileContents:
    ...

def read_solar_activity_data(file_path: str) -> dict[float, SolarActivityData]:
    """Reads a space weather data file and produces a dictionary with solar activity data for a range of epochs. Data files can be obtained from http://celestrak.com/SpaceData and should follow the legacy format.
	
	:param file_path: Path to the space weather data file.
	"""

def read_tracking_txt_file(file_name: str, column_types: list[str], comment_symbol: str='#', value_separators: str=',:\t ') -> ...:
    ...

def read_vector_history_from_file(vector_size: int, file_name: str) -> dict[float, ..., -1, 1, 0, -1, ...]:
    """Read a vector history from a file.
	
	:param vector_size:
			Size of the vector at each epoch.
	:param file_name:
			Name of the file containing the vector history.
	:return:
			Dictionary mapping epochs to the vector at the given epoch.
	"""

def set_dsn_weather_data_in_ground_stations(bodies: tudatpy.numerical_simulation.environment.expose_environment.SystemOfBodies, weather_file_names: list[str], interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=..., ground_stations_per_complex: dict[int, list[str]]={10: ['DSS-13', 'DSS-14', 'DSS-15', 'DSS-24', 'DSS-25', 'DSS-26', 'DSS-27'], 40: ['DSS-34', 'DSS-35', 'DSS-36', 'DSS-43', 'DSS-45'], 60: ['DSS-54', 'DSS-55', 'DSS-63', 'DSS-65']}, body_with_ground_stations_name: str='Earth') -> None:
    ...
ca: StaticCoefficientNames
cap: DynamicCoefficientNames
caq: DynamicCoefficientNames
car: DynamicCoefficientNames
clb: StaticCoefficientNames
cll: StaticCoefficientNames
cllp: DynamicCoefficientNames
cllq: DynamicCoefficientNames
cllr: DynamicCoefficientNames
cln: StaticCoefficientNames
clnp: DynamicCoefficientNames
clnq: DynamicCoefficientNames
clnr: DynamicCoefficientNames
cm: StaticCoefficientNames
cma: StaticCoefficientNames
cmad: DynamicCoefficientNames
cmp: DynamicCoefficientNames
cmq: DynamicCoefficientNames
cmr: DynamicCoefficientNames
cn: StaticCoefficientNames
cna: StaticCoefficientNames
cnad: DynamicCoefficientNames
cnb: StaticCoefficientNames
cnp: DynamicCoefficientNames
cnq: DynamicCoefficientNames
cnr: DynamicCoefficientNames
cy: StaticCoefficientNames
cyb: StaticCoefficientNames
cyp: DynamicCoefficientNames
cyq: DynamicCoefficientNames
cyr: DynamicCoefficientNames
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
residual_de405: TrackingDataType
second: TrackingDataType
signal_to_noise: TrackingDataType
spacecraft_id: TrackingDataType
spacecraft_transponder_delay: TrackingDataType
spectral_max: TrackingDataType
tdb_spacecraft_j2000: TrackingDataType
tdb_time_j2000: TrackingDataType
time_tag_delay: TrackingDataType
uplink_frequency: TrackingDataType
vlbi_station_name: TrackingDataType
vx_planet_frame: TrackingDataType
vy_planet_frame: TrackingDataType
vz_planet_frame: TrackingDataType
x_planet_frame: TrackingDataType
y_planet_frame: TrackingDataType
year: TrackingDataType
z_planet_frame: TrackingDataType