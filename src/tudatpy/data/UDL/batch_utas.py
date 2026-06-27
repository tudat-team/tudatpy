import json
import numpy as np
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics import environment
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import (
    model_settings,
    links,
)
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime


class StationPairObservations:
    """Time series data for a single station pair.

    Attributes
    ----------
    epochs : list[float]
        Observation epochs in TDB seconds since J2000.
    tdoa : list[float]
        Time Difference of Arrival observations in seconds.
    tdoa_uncertainties : list[float]
        TDOA measurement uncertainties in seconds.
    fdoa : list[float]
        Frequency Difference of Arrival observations in Hz.
    fdoa_uncertainties : list[float]
        FDOA measurement uncertainties in Hz.
    """

    def __init__(self):
        self.epochs = []
        self.tdoa = []
        self.tdoa_uncertainties = []
        self.fdoa = []
        self.fdoa_uncertainties = []

    def __len__(self):
        return len(self.epochs)


class UTASMetadata:
    """Strongly-typed metadata structure for UTAS observations.

    Attributes
    ----------
    target_id : str
        Identifier for the observed target (e.g., satellite catalog number).
    frequency : float
        Center frequency of the signal in Hz.
    bandwidth : float
        Signal bandwidth in Hz.
    sensor1_delay : float
        Signal arrival delay for sensor 1 in seconds.
    sensor2_delay : float
        Signal arrival delay for sensor 2 in seconds.
    data_mode : str
        Data classification: EXERCISE, REAL, SIMULATED, or TEST.
    origin : str
        Originating system identifier.
    source : str
        Data source name.
    ucts : int
        Uncorrelated track status flag.
    """

    def __init__(self):
        self.target_id = ""
        self.frequency = 0.0
        self.bandwidth = 0.0
        self.sensor1_delay = 0.0
        self.sensor2_delay = 0.0
        self.data_mode = ""
        self.origin = ""
        self.source = ""
        self.ucts = 0


class BatchUTAS:
    """Batch loader for UTAS format TDOA/FDOA observations.

    This is the main user-facing class for loading UTAS observations from JSON files
    and converting them to Tudat format for use in orbit determination.

    Only supports single-target data. If your input files contain observations of
    multiple targets, you must filter them beforehand and create separate BatchUTAS
    instances for each target. Multiple station pairs across files are supported
    (as long as all files observe the same target).

    Examples
    ----------
    Basic sequence of usage:

    Initialise and retrieve data:

    >>> batch = BatchUTAS(["observations_day1.json", "observations_day2.json"])
    >>> print(f"Target: {batch.target_id}")
    >>> print(f"Station pairs: {batch.station_pairs}")
    >>> print(f"Station names: {batch.station_names}")
    >>> print(f"Number of observations: {batch.num_observations}")

    Convert to Tudat format (creates ground stations automatically):

    >>> bodies = environment_setup.create_system_of_bodies(body_settings)
    >>> observation_collection = batch.to_tudat(bodies)

    Use custom target name instead of NORAD ID:

    >>> observation_collection = batch.to_tudat(bodies, target_name_override="MySatellite")
    """

    def __init__(self, file_paths: list[str]) -> None:
        """Construct from a list of JSON file paths.

        Parameters
        ----------
        file_paths : list[str]
            List of paths to UTAS JSON files. All files must contain observations
            of the same target. Different station pairs across files are supported.

        Raises
        ------
        RuntimeError
            If files contain multiple targets (lists all found targets in error message).
        """
        self._file_paths = file_paths
        self._metadata = UTASMetadata()
        self._metadata_initialized = False
        self._observations_by_station_pair: dict[tuple[str, str], StationPairObservations] = {}
        self._station_positions: dict[str, dict[str, float]] = {}
        self._found_targets: set[str] = set()

        self._parse_files()

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def target_id(self) -> str:
        """Target identifier (e.g., satellite catalog number)."""
        return self._metadata.target_id

    @property
    def num_observations(self) -> int:
        """Total number of observations across all station pairs."""
        return sum(len(obs) for obs in self._observations_by_station_pair.values())

    @property
    def station_pairs(self) -> list[tuple[str, str]]:
        """List of station pairs as (station1_id, station2_id) tuples."""
        return list(self._observations_by_station_pair.keys())

    @property
    def station_names(self) -> set[str]:
        """Set of unique station names across all station pairs."""
        return set(self._station_positions.keys())

    @property
    def num_station_pairs(self) -> int:
        """Number of unique station pairs."""
        return len(self._observations_by_station_pair)

    # =========================================================================
    # Public methods
    # =========================================================================

    def get_metadata(self) -> UTASMetadata:
        """Get full UTAS metadata.

        Returns
        -------
        UTASMetadata
            Metadata containing all UTAS-specific fields.
        """
        return self._metadata

    def get_observations_for_station_pair(
        self, station_pair: tuple[str, str]
    ) -> StationPairObservations:
        """Get observations for a specific station pair.

        Parameters
        ----------
        station_pair : tuple[str, str]
            The station pair as (station1_id, station2_id).

        Returns
        -------
        StationPairObservations
            Observations containing epochs, TDOA, FDOA and their uncertainties.

        Raises
        ------
        RuntimeError
            If the station pair is not found.

        Examples
        --------
        >>> obs = batch.get_observations_for_station_pair(("STATION_A", "STATION_B"))
        >>> print(f"Number of observations: {len(obs)}")
        >>> print(f"TDOA values: {obs.tdoa}")
        """
        if station_pair not in self._observations_by_station_pair:
            raise RuntimeError(
                f"Station pair {station_pair} not found. "
                f"Available pairs: {list(self._observations_by_station_pair.keys())}"
            )
        return self._observations_by_station_pair[station_pair]

    def get_all_observations_by_station_pair(self) -> dict[tuple[str, str], StationPairObservations]:
        """Get all observations organized by station pair.

        Returns
        -------
        dict[tuple[str, str], StationPairObservations]
            Dictionary mapping station pairs to their observation data.

        Examples
        --------
        >>> all_obs = batch.get_all_observations_by_station_pair()
        >>> for station_pair, obs in all_obs.items():
        ...     print(f"{station_pair}: {len(obs)} observations")
        """
        return self._observations_by_station_pair

    def to_tudat(
        self,
        bodies: environment.SystemOfBodies,
        station_body: str = "Earth",
        target_name_override: str = "",
    ) -> observations.ObservationCollection:
        """Convert to Tudat observation collection.

        This method performs all necessary setup:
        1. Ensures the station body has a compatible shape model
        2. Creates ground stations on the body
        3. Builds and returns the observation collection

        Parameters
        ----------
        bodies : SystemOfBodies
            System of bodies (will be modified to add ground stations).
        station_body : str, default="Earth"
            Name of the body on which to place ground stations.
        target_name_override : str, default=""
            Custom name for the target in link definitions. If empty, uses the
            target ID from the data (typically NORAD ID). Use this to match
            the body name in your simulation.

        Returns
        -------
        ObservationCollection
            Tudat observation collection containing TDOA and FDOA observation sets.

        Raises
        ------
        RuntimeError
            If station body has an incompatible shape model (must be OblateSpheroidBodyShapeModel).

        Examples
        --------
        Use custom target name instead of NORAD ID:

        >>> observation_collection = batch.to_tudat(bodies, target_name_override="MySatellite")
        """
        self._ensure_shape_model(bodies, station_body)
        self._create_ground_stations(bodies, station_body)
        return self._get_observation_collection(station_body, target_name_override)

    def ensure_shape_model(
        self,
        bodies: environment.SystemOfBodies,
        station_body: str = "Earth",
    ) -> None:
        """Ensure the station body has a compatible shape model.

        Creates an oblate spheroid shape model from SPICE if none exists.
        Called automatically by to_tudat().

        Parameters
        ----------
        bodies : SystemOfBodies
            System of bodies.
        station_body : str, default="Earth"
            Name of the body to check/modify.

        Raises
        ------
        RuntimeError
            If body has an incompatible (non-oblate-spheroid) shape model.
        """
        try:
            bodies.get_body(station_body)
        except Exception:
            bodies.add_body(environment.Body(), station_body)

        body = bodies.get_body(station_body)

        if body.get_shape_model() is None:
            shape_model = environment_setup.create_body_shape_model(
                environment_setup.from_spice_oblate_spherical_body_shape_settings(),
                station_body,
            )
            body.set_shape_model(shape_model)

    def create_ground_stations(
        self,
        bodies: environment.SystemOfBodies,
        station_body: str = "Earth",
    ) -> list[str]:
        """Create ground stations on the specified body.

        Called automatically by to_tudat().

        Parameters
        ----------
        bodies : SystemOfBodies
            System of bodies (modified in place).
        station_body : str, default="Earth"
            Body on which to create stations.

        Returns
        -------
        list[str]
            Names of the created stations.
        """
        station_names = []
        body = bodies.get_body(station_body)

        for station_name, position in self._station_positions.items():
            tudat_position = self._convert_to_tudat_geodetic(position)

            settings = environment_setup.ground_station.basic_station(
                station_name=station_name,
                station_nominal_position=tudat_position.tolist(),
            )

            if station_name not in body.ground_station_list:
                environment_setup.add_ground_station(body, settings)
                station_names.append(station_name)

        return station_names

    def get_link_definitions(
        self,
        station_body: str = "Earth",
        target_name_override: str = "",
    ) -> list:
        """Get the link definitions for all station pairs in this batch.

        Parameters
        ----------
        station_body : str, default="Earth"
            Name of body hosting ground stations.
        target_name_override : str, default=""
            Custom name for the target in link definitions. If empty, uses the
            target ID from the data (typically NORAD ID).

        Returns
        -------
        list[LinkDefinition]
            Link definitions with receiver, receiver2, and transmitter link ends,
            one per station pair.
        """
        target_name = target_name_override if target_name_override else self._metadata.target_id

        link_definitions = []
        for station_pair in self.station_pairs:
            link_ends = dict()
            link_ends[links.receiver] = links.body_reference_point_link_end_id(
                station_body, station_pair[0]
            )
            link_ends[links.receiver2] = links.body_reference_point_link_end_id(
                station_body, station_pair[1]
            )
            link_ends[links.transmitter] = links.body_origin_link_end_id(target_name)

            link_definitions.append(links.link_definition(link_ends))

        return link_definitions

    def get_observation_collection(
        self,
        station_body: str = "Earth",
        target_name_override: str = "",
    ) -> observations.ObservationCollection:
        """Get observation collection without modifying bodies.

        Use this if you've already created ground stations manually.

        Parameters
        ----------
        station_body : str, default="Earth"
            Name of body hosting ground stations.
        target_name_override : str, default=""
            Custom name for the target in link definitions. If empty, uses the
            target ID from the data (typically NORAD ID).

        Returns
        -------
        ObservationCollection
            Observation collection with TDOA and FDOA sets.
        """
        return self._get_observation_collection(station_body, target_name_override)

    # =========================================================================
    # Private methods
    # =========================================================================

    def _parse_files(self) -> None:
        """Parse all JSON files and populate internal data structures."""
        for path in self._file_paths:
            self._parse_file(path)

        if not self._observations_by_station_pair:
            raise RuntimeError("BatchUTAS: No observations found in provided files")

    def _parse_file(self, file_path: str) -> None:
        """Parse a single UTAS JSON file.

        Parameters
        ----------
        file_path : str
            Path to the JSON file.
        """
        with open(file_path, "r") as f:
            try:
                data = json.load(f)
            except json.JSONDecodeError as e:
                raise RuntimeError(f"BatchUTAS: JSON parse error in {file_path}: {e}")

        # Handle different JSON structures
        if isinstance(data, list):
            observations_list = data
        elif isinstance(data, dict) and "observations" in data and isinstance(data["observations"], list):
            observations_list = data["observations"]
        else:
            raise RuntimeError(
                f"BatchUTAS: Unexpected JSON structure in {file_path}. "
                "Expected a list of observations or a dict with 'observations' key."
            )

        if not observations_list:
            print(f"WARNING: BatchUTAS: File {file_path} contains no observations")
            return

        self._parse_observation_array(observations_list, file_path)

    def _parse_observation_array(self, observations_list: list, file_path: str) -> None:
        """Parse an array of observation records.

        Parameters
        ----------
        observations_list : list
            List of observation dictionaries.
        file_path : str
            Path to the file (for error messages).
        """
        first_obs = observations_list[0]

        # Extract target and validate single-target constraint
        target_id = self._get_string_or_number(first_obs, "satNo")
        self._found_targets.add(target_id)
        self._validate_single_target(target_id, file_path)

        # Extract station info for this file
        station1_id = self._get_required(first_obs, "origSensorId1", str)
        station2_id = self._get_required(first_obs, "origSensorId2", str)
        station_pair = (station1_id, station2_id)

        # Extract station positions
        station1_pos = {
            "altitude": self._get_required(first_obs, "senlat", float),
            "latitude": self._get_required(first_obs, "senlat", float),
            "longitude": self._get_required(first_obs, "senlon", float),
        }
        station2_pos = {
            "altitude": self._get_required(first_obs, "sen2lat", float),
            "latitude": self._get_required(first_obs, "sen2lat", float),
            "longitude": self._get_required(first_obs, "sen2lon", float),
        }

        # Store station positions
        self._station_positions[station1_id] = station1_pos
        self._station_positions[station2_id] = station2_pos

        # Initialize metadata from first file
        if not self._metadata_initialized:
            self._metadata.target_id = target_id
            self._metadata.frequency = self._get_required(first_obs, "frequency", float)
            self._metadata.bandwidth = self._get_optional(first_obs, "bandwidth", 0.0, float)
            self._metadata.sensor1_delay = self._get_optional(first_obs, "sensor1Delay", 0.0, float)
            self._metadata.sensor2_delay = self._get_optional(first_obs, "sensor2Delay", 0.0, float)
            self._metadata.data_mode = self._get_optional(first_obs, "dataMode", "", str)
            self._metadata.origin = self._get_optional(first_obs, "origin", "", str)
            self._metadata.source = self._get_optional(first_obs, "source", "", str)
            self._metadata.ucts = self._get_optional(first_obs, "ucts", 0, int)
            self._metadata_initialized = True

        # Get or create observation storage for this station pair
        if station_pair not in self._observations_by_station_pair:
            self._observations_by_station_pair[station_pair] = StationPairObservations()

        station_obs = self._observations_by_station_pair[station_pair]

        # Parse time series data
        for obs in observations_list:
            # Check target consistency
            obs_target = self._get_string_or_number(obs, "satNo")
            self._found_targets.add(obs_target)
            self._validate_single_target(obs_target, file_path)

            # Time
            ob_time = self._get_required(obs, "obTime", str)
            epoch = self._convert_iso_to_epoch(ob_time)
            station_obs.epochs.append(epoch)

            # TDOA
            station_obs.tdoa.append(self._get_required(obs, "tdoa", float))
            station_obs.tdoa_uncertainties.append(
                self._get_optional(obs, "tdoaUnc", 0.0, float)
            )

            # FDOA
            station_obs.fdoa.append(self._get_required(obs, "fdoa", float))
            station_obs.fdoa_uncertainties.append(
                self._get_optional(obs, "fdoaUnc", 0.0, float)
            )

    def _validate_single_target(self, new_target_id: str, file_path: str) -> None:
        """Validate that only a single target is present across all files.

        Parameters
        ----------
        new_target_id : str
            The target ID found in the current observation.
        file_path : str
            Path to the file (for error messages).

        Raises
        ------
        RuntimeError
            If multiple targets are detected.
        """
        if len(self._found_targets) > 1:
            targets_str = ", ".join(f"'{t}'" for t in sorted(self._found_targets))
            raise RuntimeError(
                f"BatchUTAS: Multiple targets detected. BatchUTAS only supports single-target data.\n"
                f"Found targets: {targets_str}\n"
                f"Please create separate BatchUTAS instances for each target by filtering input files.\n"
                f"Error occurred while parsing: {file_path}"
            )

    def _convert_iso_to_epoch(self, iso_time: str) -> float:
        """Convert an ISO time string to TDB seconds since J2000.

        Parameters
        ----------
        iso_time : str
            ISO format time string (e.g., "2023-01-01T00:00:00.000Z").

        Returns
        -------
        float
            Time in TDB seconds since J2000.
        """
        # Strip trailing 'Z' if present
        time_str = iso_time
        if time_str and time_str[-1] in ("Z", "z"):
            time_str = time_str[:-1]

        if len(time_str) < 19:
            raise RuntimeError(f"BatchUTAS: Invalid time format: {iso_time}")

        try:
            dt = DateTime.from_iso_string(time_str)
            time_utc = dt.epoch()

            time_scale_converter = time_representation.default_time_scale_converter()
            time_tdb = time_scale_converter.convert_time(
                input_scale=time_representation.utc_scale,
                output_scale=time_representation.tdb_scale,
                input_value=time_utc,
            )
            return time_tdb
        except Exception as e:
            raise RuntimeError(f"BatchUTAS: Failed to parse time '{time_str}': {e}")

    def _convert_to_tudat_geodetic(self, position: dict[str, float]) -> np.ndarray:
        """Convert a geodetic position to Tudat format (altitude[m], latitude[rad], longitude[rad]).

        Assumes input is in meters and degrees.

        Parameters
        ----------
        position : dict[str, float]
            Dictionary with 'altitude', 'latitude', 'longitude' keys.

        Returns
        -------
        np.ndarray
            3-element array [altitude[m], latitude[rad], longitude[rad]].
        """
        return np.array([
            position["altitude"],
            position["latitude"] * np.pi / 180.0,
            position["longitude"] * np.pi / 180.0,
        ])

    def _ensure_shape_model(self, bodies: environment.SystemOfBodies, station_body: str) -> None:
        """Ensure the station body has a compatible shape model.

        Parameters
        ----------
        bodies : SystemOfBodies
            System of bodies.
        station_body : str
            Name of body to check/modify.
        """
        try:
            bodies.get_body(station_body)
        except Exception:
            bodies.add_body(environment.Body(), station_body)

        body = bodies.get_body(station_body)

        if body.get_shape_model() is None:
            shape_model = environment_setup.create_body_shape_model(
                environment_setup.from_spice_oblate_spherical_body_shape_settings(),
                station_body,
            )
            body.set_shape_model(shape_model)

    def _create_ground_stations(self, bodies: environment.SystemOfBodies, station_body: str) -> None:
        """Create ground stations on the specified body.

        Parameters
        ----------
        bodies : SystemOfBodies
            System of bodies (modified in place).
        station_body : str
            Body on which to create stations.
        """
        body = bodies.get_body(station_body)

        for station_name, position in self._station_positions.items():
            tudat_pos = self._convert_to_tudat_geodetic(position)

            settings = environment_setup.ground_station.basic_station(
                station_name=station_name,
                station_nominal_position=tudat_pos.tolist(),
            )

            if station_name not in body.ground_station_list:
                environment_setup.add_ground_station(body, settings)

    def _get_observation_collection(
        self,
        station_body: str,
        target_name_override: str,
    ) -> observations.ObservationCollection:
        """Build the Tudat observation collection.

        Parameters
        ----------
        station_body : str
            Name of body hosting ground stations.
        target_name_override : str
            Custom name for the target in link definitions.

        Returns
        -------
        ObservationCollection
            Observation collection with TDOA and FDOA sets.
        """
        target_name = target_name_override if target_name_override else self._metadata.target_id

        observation_set_list = []

        for station_pair in self.station_pairs:
            station_obs = self._observations_by_station_pair[station_pair]

            # Build link ends
            link_ends = dict()
            link_ends[links.receiver] = links.body_reference_point_link_end_id(
                station_body, station_pair[0]
            )
            link_ends[links.receiver2] = links.body_reference_point_link_end_id(
                station_body, station_pair[1]
            )
            link_ends[links.transmitter] = links.body_origin_link_end_id(target_name)
            link_def = links.link_definition(link_ends)

            # Build observation data
            observation_times = []
            tdoa_observations = []
            fdoa_observations = []

            for i in range(len(station_obs)):
                observation_times.append(station_obs.epochs[i])

                tdoa_entry = np.array([[station_obs.tdoa[i]]])
                tdoa_observations.append(tdoa_entry)

                fdoa_entry = np.array([[station_obs.fdoa[i]]])
                fdoa_observations.append(fdoa_entry)

            # Create TDOA observation set
            tdoa_set = observations.create_single_observation_set(
                model_settings.differenced_time_of_arrival_type,
                link_def,
                tdoa_observations,
                observation_times,
                links.receiver,
            )
            observation_set_list.append(tdoa_set)

            # Create FDOA observation set
            fdoa_set = observations.create_single_observation_set(
                model_settings.differenced_frequency_of_arrival_type,
                link_def,
                fdoa_observations,
                observation_times,
                links.receiver,
            )
            observation_set_list.append(fdoa_set)

        return observations.ObservationCollection(observation_set_list)

    # =========================================================================
    # Static helper methods
    # =========================================================================

    @staticmethod
    def _get_required(obj: dict, key: str, expected_type: type):
        """Get a required field from a dictionary.

        Parameters
        ----------
        obj : dict
            The dictionary to extract from.
        key : str
            The key to look up.
        expected_type : type
            The expected type of the value.

        Returns
        -------
        The value from the dictionary.

        Raises
        ------
        RuntimeError
            If the key is not found or has the wrong type.
        """
        if key not in obj:
            raise RuntimeError(f"BatchUTAS: Required field '{key}' not found")

        val = obj[key]
        if not isinstance(val, expected_type):
            # Allow int where float is expected
            if expected_type is float and isinstance(val, (int, np.integer)):
                return float(val)
            raise RuntimeError(
                f"BatchUTAS: Field '{key}' has wrong type. "
                f"Expected {expected_type.__name__}, got {type(val).__name__}"
            )
        return val

    @staticmethod
    def _get_optional(obj: dict, key: str, default, expected_type: type):
        """Get an optional field from a dictionary.

        Parameters
        ----------
        obj : dict
            The dictionary to extract from.
        key : str
            The key to look up.
        default : any
            Default value if key is not found.
        expected_type : type
            The expected type of the value.

        Returns
        -------
        The value from the dictionary or the default.
        """
        if key not in obj:
            return default

        val = obj[key]
        if not isinstance(val, expected_type):
            if expected_type is float and isinstance(val, (int, np.integer)):
                return float(val)
            return default
        return val

    @staticmethod
    def _get_string_or_number(obj: dict, key: str) -> str:
        """Get a field that may be a string or number, returned as string.

        Parameters
        ----------
        obj : dict
            The dictionary to extract from.
        key : str
            The key to look up.

        Returns
        -------
        str
            The value as a string.

        Raises
        ------
        RuntimeError
            If the key is not found or the value is not a string or number.
        """
        if key not in obj:
            raise RuntimeError(f"BatchUTAS: Required field '{key}' not found")

        val = obj[key]
        if isinstance(val, str):
            return val
        elif isinstance(val, (int, np.integer)):
            return str(int(val))
        elif isinstance(val, (float, np.floating)):
            return str(float(val))
        else:
            raise RuntimeError(
                f"BatchUTAS: Field '{key}' must be string or number, got {type(val).__name__}"
            )
