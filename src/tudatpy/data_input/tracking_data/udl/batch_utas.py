"""Unified Data Library TDOA/FDOA tracking-data reader."""

import json

import numpy as np

from tudatpy.astro import element_conversion
from tudatpy.astro.time_representation import DateTime
from tudatpy.data_input.tracking_data import TrackingData, TrackingSupplementaryData
from tudatpy.dynamics import environment, environment_setup


class StationPairObservations:
    """Time series data for a single station pair.

    Attributes
    ----------
    epochs : list[float]
        Observation epochs in UTC seconds since J2000.
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
    """Metadata from the first UTAS observation in a batch.

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
    uct : bool
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
        self.uct = False


class BatchUTAS:
    """Batch loader for UTAS format TDOA/FDOA observations.

    This is the main user-facing class for loading UTAS observations from JSON files
    and converting them to Tudat format for use in orbit determination.

    Only supports single-target data. If your input files contain observations of
    multiple targets, you must filter them beforehand and create separate BatchUTAS
    instances for each target. Multiple station pairs across files are supported
    (as long as all files observe the same target).
    The reader expects correlated UTAS records containing ``satNo``, both TDOA and
    FDOA values, and explicit, static geodetic positions for both sensors.

    Examples
    ----------
    Basic sequence of usage:

    Initialise and retrieve data:

    >>> batch = BatchUTAS(["observations_day1.json", "observations_day2.json"])
    >>> print(f"Target: {batch.target_id}")
    >>> print(f"Station pairs: {batch.station_pairs}")
    >>> print(f"Station names: {batch.station_names}")
    >>> print(f"Number of observations: {batch.num_observations}")

    Convert to the common tracking-data format:

    >>> tracking_data, supplementary_data = batch.to_tracking_dataset()

    Use custom target name instead of NORAD ID:

    >>> tracking_data, _ = batch.to_tracking_dataset(spacecraft_name="MySatellite")
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
        """Total number of TDOA/FDOA record pairs across all station pairs."""
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
        """Get metadata from the first UTAS observation in the batch.

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

    def get_all_observations_by_station_pair(
        self,
    ) -> dict[tuple[str, str], StationPairObservations]:
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

    def to_tracking_dataset(
        self, station_body: str = "Earth", spacecraft_name: str | None = None
    ) -> tuple[list[TrackingData], list[TrackingSupplementaryData]]:
        """Convert to the common tracking-data containers.

        Parameters
        ----------
        station_body : str, default="Earth"
            Name of the body on which to place ground stations.
        spacecraft_name : str | None, default=None
            Custom name for the target in link definitions. If omitted, uses the
            target ID from the data (typically NORAD ID). Use this to match
            the body name in your simulation.

        Returns
        -------
        tuple[list[TrackingData], list[TrackingSupplementaryData]]
            TDOA and FDOA tracking data, and an empty supplementary-data list.

        Examples
        --------
        Use custom target name instead of NORAD ID:

        >>> tracking_data, _ = batch.to_tracking_dataset(spacecraft_name="MySatellite")
        """
        return self._get_tracking_data(station_body, spacecraft_name), []

    def create_ground_stations(
        self,
        bodies: environment.SystemOfBodies,
        station_body: str = "Earth",
    ) -> list[str]:
        """Create ground stations on the specified body.

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
        if body.shape_model is None:
            raise RuntimeError(
                f"BatchUTAS: Body '{station_body}' needs a shape model to create "
                "stations from geodetic coordinates"
            )

        for station_name, position in self._station_positions.items():
            tudat_position = self._convert_to_tudat_geodetic(position)

            settings = environment_setup.ground_station.basic_station(
                station_name=station_name,
                station_nominal_position=tudat_position.tolist(),
                station_position_element_type=element_conversion.geodetic_position_type,
            )

            if station_name not in body.ground_station_list:
                environment_setup.add_ground_station(body, settings)
                station_names.append(station_name)

        return station_names

    def get_link_definitions(
        self,
        station_body: str = "Earth",
        spacecraft_name: str | None = None,
    ) -> list:
        """Get the link definitions for all station pairs in this batch.

        Parameters
        ----------
        station_body : str, default="Earth"
            Name of body hosting ground stations.
        spacecraft_name : str | None, default=None
            Custom name for the target in link definitions. If omitted, uses the
            target ID from the data (typically NORAD ID).

        Returns
        -------
        list[list]
            Plain link definitions with receiver, receiver_2, and transmitter link
            ends, one per station pair.
        """
        target_name = spacecraft_name or self._metadata.target_id

        link_definitions = []
        for station_pair in self.station_pairs:
            link_definitions.append(
                [
                    ((target_name, ""), "transmitter"),
                    ((station_body, station_pair[0]), "receiver"),
                    ((station_body, station_pair[1]), "receiver_2"),
                ]
            )

        return link_definitions

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
        elif (
            isinstance(data, dict)
            and "observations" in data
            and isinstance(data["observations"], list)
        ):
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
        if station1_id == station2_id:
            raise RuntimeError("BatchUTAS: A station cannot be paired with itself")
        station_pair = (station1_id, station2_id)

        # Extract station positions
        station1_pos = {
            "altitude": self._get_required(first_obs, "senalt", float),
            "latitude": self._get_required(first_obs, "senlat", float),
            "longitude": self._get_required(first_obs, "senlon", float),
        }
        station2_pos = {
            "altitude": self._get_required(first_obs, "sen2alt", float),
            "latitude": self._get_required(first_obs, "sen2lat", float),
            "longitude": self._get_required(first_obs, "sen2lon", float),
        }

        # Store station positions, rejecting moving sensors because Tudat ground
        # stations are static unless an explicit motion model is provided.
        for station_id, position in (
            (station1_id, station1_pos),
            (station2_id, station2_pos),
        ):
            if station_id in self._station_positions:
                if self._station_positions[station_id] != position:
                    raise RuntimeError(
                        f"BatchUTAS: Position of station '{station_id}' differs across files"
                    )
            else:
                self._station_positions[station_id] = position

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
            if first_obs.get("uct") is not None:
                self._metadata.uct = self._get_required(first_obs, "uct", bool)
            else:
                self._metadata.uct = bool(self._get_optional(first_obs, "ucts", 0, int))
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

            obs_station_pair = (
                self._get_required(obs, "origSensorId1", str),
                self._get_required(obs, "origSensorId2", str),
            )
            if obs_station_pair != station_pair:
                raise RuntimeError(
                    "BatchUTAS: Multiple station pairs detected in one file. "
                    f"Expected {station_pair}, found {obs_station_pair} in {file_path}."
                )
            position_keys = ("senalt", "senlat", "senlon", "sen2alt", "sen2lat", "sen2lon")
            if any(
                self._get_required(obs, key, float) != self._get_required(first_obs, key, float)
                for key in position_keys
            ):
                raise RuntimeError(f"BatchUTAS: Station positions vary within {file_path}")

            # Time
            ob_time = self._get_required(obs, "obTime", str)
            epoch = self._convert_iso_to_epoch(ob_time)
            station_obs.epochs.append(epoch)

            # TDOA
            station_obs.tdoa.append(self._get_required(obs, "tdoa", float))
            station_obs.tdoa_uncertainties.append(self._get_optional(obs, "tdoaUnc", 0.0, float))

            # FDOA
            station_obs.fdoa.append(self._get_required(obs, "fdoa", float))
            station_obs.fdoa_uncertainties.append(self._get_optional(obs, "fdoaUnc", 0.0, float))

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
        """Convert an ISO time string to UTC seconds since J2000.

        Parameters
        ----------
        iso_time : str
            ISO format time string (e.g., "2023-01-01T00:00:00.000Z").

        Returns
        -------
        float
            Time in UTC seconds since J2000.
        """
        # Strip trailing 'Z' if present
        time_str = iso_time
        if time_str and time_str[-1] in ("Z", "z"):
            time_str = time_str[:-1]

        if len(time_str) < 19:
            raise RuntimeError(f"BatchUTAS: Invalid time format: {iso_time}")

        try:
            dt = DateTime.from_iso_string(time_str)
            return dt.epoch()
        except Exception as e:
            raise RuntimeError(f"BatchUTAS: Failed to parse time '{time_str}': {e}")

    def _convert_to_tudat_geodetic(self, position: dict[str, float]) -> np.ndarray:
        """Convert a geodetic position to Tudat format (altitude[m], latitude[rad], longitude[rad]).

        Assumes input altitude is in kilometres and angles are in degrees, as
        defined by the UDL schema.

        Parameters
        ----------
        position : dict[str, float]
            Dictionary with 'altitude', 'latitude', 'longitude' keys.

        Returns
        -------
        np.ndarray
            3-element array [altitude[m], latitude[rad], longitude[rad]].
        """
        if not -90.0 <= position["latitude"] <= 90.0:
            raise RuntimeError("BatchUTAS: Station latitude must be in [-90, 90] degrees")
        if not -180.0 <= position["longitude"] <= 180.0:
            raise RuntimeError("BatchUTAS: Station longitude must be in [-180, 180] degrees")

        return np.array(
            [
                1000.0 * position["altitude"],
                position["latitude"] * np.pi / 180.0,
                position["longitude"] * np.pi / 180.0,
            ]
        )

    def _get_tracking_data(
        self,
        station_body: str,
        spacecraft_name: str | None,
    ) -> list[TrackingData]:
        """Build the common Tudat tracking-data containers.

        Parameters
        ----------
        station_body : str
            Name of body hosting ground stations.
        spacecraft_name : str | None
            Custom name for the target in link definitions.

        Returns
        -------
        list[TrackingData]
            One TDOA and one FDOA container per station pair.
        """
        tracking_data_list = []

        for station_pair, link_ends in zip(
            self.station_pairs,
            self.get_link_definitions(station_body, spacecraft_name),
        ):
            station_obs = self._observations_by_station_pair[station_pair]

            # UDL uses sensor2 - sensor1, while Tudat's differenced models use
            # receiver - receiver_2. Negate the source values to match Tudat.
            for observable_type, values in (
                ("DifferencedTimeOfArrival", station_obs.tdoa),
                ("DifferencedFrequencyOfArrival", station_obs.fdoa),
            ):
                tracking_data_list.append(
                    TrackingData(
                        observable_type=observable_type,
                        link_ends=link_ends,
                        observations=[np.array([-value]) for value in values],
                        epochs=station_obs.epochs,
                        reference_link_end="receiver",
                        time_scale="UTC",
                    )
                )

        return tracking_data_list

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
            if (
                expected_type is float
                and isinstance(val, (int, np.integer))
                and not isinstance(val, (bool, np.bool_))
            ):
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
        if val is None:
            return default
        if not isinstance(val, expected_type):
            if (
                expected_type is float
                and isinstance(val, (int, np.integer))
                and not isinstance(val, (bool, np.bool_))
            ):
                return float(val)
            raise RuntimeError(
                f"BatchUTAS: Field '{key}' has wrong type. "
                f"Expected {expected_type.__name__}, got {type(val).__name__}"
            )
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
        if isinstance(val, (bool, np.bool_)):
            raise RuntimeError(f"BatchUTAS: Field '{key}' must be string or number")
        if isinstance(val, str) and val:
            return val
        elif isinstance(val, (int, np.integer)):
            return str(int(val))
        elif isinstance(val, (float, np.floating)) and np.isfinite(val) and val.is_integer():
            return str(int(val))
        else:
            raise RuntimeError(
                f"BatchUTAS: Field '{key}' must be string or number, got {type(val).__name__}"
            )


def read_utas_data(
    utas_file_names: list[str],
    spacecraft_name: str | None = None,
    station_body: str = "Earth",
) -> tuple[list[TrackingData], list[TrackingSupplementaryData]]:
    """Read single-target UDL TDOA/FDOA JSON files.

    The receiving stations must exist on ``station_body`` before the returned
    tracking data is converted to an observation collection. Use
    :class:`BatchUTAS` directly to create stations from the positions in the
    source records.

    Parameters
    ----------
    utas_file_names : list[str]
        Paths to UDL JSON files.
    spacecraft_name : str | None, default=None
        Name of the transmitting body. The UDL target ID is used if omitted.
    station_body : str, default="Earth"
        Body on which both receiving stations are located.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        TDOA and FDOA tracking data, and an empty supplementary-data list.
    """
    return BatchUTAS(utas_file_names).to_tracking_dataset(station_body, spacecraft_name)
