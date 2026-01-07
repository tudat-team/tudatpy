.. _unified_data_library:

``unified_data_library``
========================
Interface for loading observation data in the Unified Data Library (UDL) format.

This module provides classes for reading unified data library (UDL) observation data, currently only supporting JSON file formatting. 
Currently supported vendors include:

- UTAS (University of Tasmania)
   * Differenced Frequency of Arrival (FDOA)
   * Differenced Time of Arrival (TDOA)

Classes
-------

.. class:: UTASMetadata

   Strongly-typed metadata structure for UTAS observations.

   Contains UTAS-specific fields including target ID, signal parameters, and
   data provenance. Station information is stored separately as there may be
   multiple station pairs across files.

   .. attribute:: target_id
      :type: str

      Identifier for the observed target (satellite catalog number).

   .. attribute:: frequency
      :type: float

      Center frequency of the signal in Hz.

   .. attribute:: bandwidth
      :type: float

      Signal bandwidth in Hz.

   .. attribute:: sensor1_delay
      :type: float

      Signal arrival delay for sensor 1 in seconds.

   .. attribute:: sensor2_delay
      :type: float

      Signal arrival delay for sensor 2 in seconds.

   .. attribute:: data_mode
      :type: str

      Data classification: EXERCISE, REAL, SIMULATED, or TEST.

   .. attribute:: origin
      :type: str

      Originating system identifier.

   .. attribute:: source
      :type: str

      Data source name.

   .. attribute:: ucts
      :type: str

      Uncorrelated track status flag.


.. class:: BatchUTAS

   Batch loader for UTAS format TDOA/FDOA observations.

   This is the main user-facing class for loading UTAS observations from JSON files
   and converting them to Tudat format for use in orbit determination.

   **Important:** This class only supports single-target data. If your input files
   contain observations of multiple targets, you must filter them beforehand and
   create separate BatchUTAS instances for each target. You can also instantiate 
   multiple separate BatchUTAS classes, one per target. Multiple station pairs
   across files are supported (as long as all files observe the same target).

   **Example:**

   .. code-block:: python

      batch = BatchUTAS(["observations_day1.json", "observations_day2.json"])
      print(f"Target: {batch.target_id}")
      print(f"Station pairs: {batch.station_pairs}")
      print(f"Station names: {batch.station_names}")
      print(f"Number of observations: {batch.num_observations}")

      # Convert to Tudat format (creates ground stations automatically)
      observation_collection = batch.to_tudat(bodies)

   .. method:: __init__(file_paths)

      Construct from a list of JSON file paths.

      :param file_paths: List of paths to UTAS JSON files. All files must contain observations
                         of the same target. Different station pairs across files are supported.
      :type file_paths: list[str]
      :raises RuntimeError: If files contain multiple targets (lists all found targets in error message).

   .. method:: to_tudat(bodies, station_body="Earth", target_name_override="")

      Convert to Tudat observation collection.

      This method performs all necessary setup:

      1. Ensures the station body has a compatible shape model
      2. Creates ground stations on the body
      3. Builds and returns the observation collection

      :param bodies: System of bodies (will be modified to add ground stations).
      :type bodies: SystemOfBodies
      :param station_body: Name of the body on which to place ground stations.
      :type station_body: str
      :param target_name_override: Custom name for the target in link definitions. If empty, uses the
                                   target ID from the data (typically NORAD ID). Use this to match
                                   the body name in your simulation.
      :type target_name_override: str
      :returns: Tudat observation collection containing TDOA and FDOA observation sets.
      :rtype: ObservationCollection
      :raises RuntimeError: If station body has an incompatible shape model (must be OblateSpheroidBodyShapeModel).

   .. method:: ensure_shape_model(bodies, station_body="Earth")

      Ensure the station body has a compatible shape model.

      Creates an oblate spheroid shape model from SPICE if none exists.
      Called automatically by :meth:`to_tudat`.

      :param bodies: System of bodies.
      :type bodies: SystemOfBodies
      :param station_body: Name of the body to check/modify.
      :type station_body: str
      :raises RuntimeError: If body has an incompatible (non-oblate-spheroid) shape model.

   .. method:: create_ground_stations(bodies, station_body="Earth")

      Create ground stations on the specified body.

      Called automatically by :meth:`to_tudat`.

      :param bodies: System of bodies (modified in place).
      :type bodies: SystemOfBodies
      :param station_body: Body on which to create stations.
      :type station_body: str
      :returns: Names of the created stations.
      :rtype: list[str]

   .. method:: get_link_definitions(station_body="Earth", target_name_override="")

      Get the link definitions for all station pairs in this batch.

      :param station_body: Name of body hosting ground stations.
      :type station_body: str
      :param target_name_override: Custom name for the target in link definitions. If empty, uses the
                                   target ID from the data (typically NORAD ID).
      :type target_name_override: str
      :returns: Link definitions with receiver, receiver2, and transmitter link ends,
                one per station pair.
      :rtype: list[LinkDefinition]

   .. method:: get_observation_collection(station_body="Earth", target_name_override="")

      Get observation collection without modifying bodies.

      Use this if you have already created ground stations manually.

      :param station_body: Name of body hosting ground stations.
      :type station_body: str
      :param target_name_override: Custom name for the target in link definitions. If empty, uses the
                                   target ID from the data (typically NORAD ID).
      :type target_name_override: str
      :returns: Observation collection with TDOA and FDOA sets.
      :rtype: ObservationCollection

   .. method:: get_metadata()

      Get full UTAS metadata.

      :returns: Metadata containing all UTAS-specific fields.
      :rtype: UTASMetadata

   .. attribute:: target_id
      :type: str

      Target identifier (e.g., satellite catalog number). Read-only.

   .. attribute:: num_observations
      :type: int

      Total number of observations across all station pairs. Read-only.

   .. attribute:: station_pairs
      :type: list[tuple[str, str]]

      List of station pairs as (station1_id, station2_id) tuples. Read-only.

   .. attribute:: station_names
      :type: set[str]

      Set of unique station names across all station pairs. Read-only.

   .. attribute:: num_station_pairs
      :type: int

      Number of unique station pairs. Read-only.

   .. attribute:: epochs
      :type: numpy.ndarray

      Observation epochs in TDB seconds since J2000. Read-only.

   .. attribute:: tdoa_observations
      :type: numpy.ndarray

      TDOA observations in seconds. Read-only.

   .. attribute:: tdoa_uncertainties
      :type: numpy.ndarray

      TDOA uncertainties in seconds. Read-only.

   .. attribute:: fdoa_observations
      :type: numpy.ndarray

      FDOA observations in Hz. Read-only.

   .. attribute:: fdoa_uncertainties
      :type: numpy.ndarray

      FDOA uncertainties in Hz. Read-only.


.. class:: BatchUTAS_time_object

   Batch loader for UTAS format observations using Time object precision.

   This is an alternative version of :class:`BatchUTAS` that uses ``Time`` instead of ``double``
   for the time type.

   This class has the same interface as :class:`BatchUTAS`. See :class:`BatchUTAS` for
   detailed documentation of all methods and attributes.

