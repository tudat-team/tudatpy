/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "expose_station_history.h"
#include <pybind11/stl.h>
#include "tudat/io/readSinexStationHistory.h"
namespace py = pybind11;
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace station_history
{
namespace tio = tudat::input_output;
void expose_station_history( py::module& m )
{
    py::enum_< tio::SinexEpochMeaning >( m, "SinexEpochMeaning", R"doc(

Enumeration of meanings of a SINEX epoch. The zero epoch ``00:000:00000`` has a meaning that depends on the field in
which it is given; it is never converted to a date.

      )doc" )
            .value( "specified", tio::SinexEpochMeaning::specified, "A calendar epoch." )
            .value( "unspecified", tio::SinexEpochMeaning::unspecified, "Zero epoch, meaning that the epoch is not given." )
            .value( "header_start_default",
                    tio::SinexEpochMeaning::headerStartDefault,
                    "Zero epoch, meaning the start epoch of the file header." )
            .value( "header_end_default",
                    tio::SinexEpochMeaning::headerEndDefault,
                    "Zero epoch, meaning the end epoch of the file header." )
            .value( "open_start", tio::SinexEpochMeaning::openStart, "Zero epoch, meaning an interval that is open at its start." )
            .value( "open_end", tio::SinexEpochMeaning::openEnd, "Zero epoch, meaning an interval that is open at its end." );
    py::enum_< tio::SinexEccentricityBasis >( m, "SinexEccentricityBasis", R"doc(

Enumeration of the component orders of a SINEX eccentricity.

      )doc" )
            .value( "une", tio::SinexEccentricityBasis::une, "Up, north and east components." )
            .value( "xyz", tio::SinexEccentricityBasis::xyz, "Body-fixed Cartesian components." );
    py::enum_< tio::SinexEquipmentKind >( m, "SinexEquipmentKind", R"doc(

Enumeration of equipment types in a SINEX file.

      )doc" )
            .value( "receiver", tio::SinexEquipmentKind::receiver, "Receiver (``SITE/RECEIVER`` block)." )
            .value( "antenna", tio::SinexEquipmentKind::antenna, "Antenna (``SITE/ANTENNA`` block)." );

    py::class_< tio::SinexSiteKey >( m, "SinexSiteKey", R"doc(

Identification of a site (monument) in a SINEX file, by its four-character site code and its point code.

      )doc" )
            .def( py::init( []( std::string site, std::string point ) { return tio::SinexSiteKey{ site, point }; } ),
                  py::arg( "site_code" ),
                  py::arg( "point_code" ),
                  R"doc(

         Constructor.

         Parameters
         ----------
         site_code : str
             Four-character site code (for instance ``ABMF``).
         point_code : str
             Point code (for instance ``A``).

      )doc" )
            .def_readonly( "site_code", &tio::SinexSiteKey::siteCode, R"doc(

         **read-only**

         Four-character site code.

         :type: str
      )doc" )
            .def_readonly( "point_code", &tio::SinexSiteKey::pointCode, R"doc(

         **read-only**

         Point code.

         :type: str
      )doc" );
    py::class_< tio::SinexSolutionKey >( m, "SinexSolutionKey", R"doc(

Identification of a solution for a site in a SINEX file. A new solution number is used after each discontinuity of
the station position or velocity.

      )doc" )
            .def( py::init( []( std::string site, std::string point, std::string solution ) {
                      return tio::SinexSolutionKey{ { site, point }, solution };
                  } ),
                  py::arg( "site_code" ),
                  py::arg( "point_code" ),
                  py::arg( "solution_id" ),
                  R"doc(

         Constructor.

         Parameters
         ----------
         site_code : str
             Four-character site code.
         point_code : str
             Point code.
         solution_id : str
             Solution number, as given in the file.

      )doc" )
            .def_readonly( "site", &tio::SinexSolutionKey::site, R"doc(

         **read-only**

         Site of the solution.

         :type: SinexSiteKey
      )doc" )
            .def_readonly( "solution_id", &tio::SinexSolutionKey::solutionId, R"doc(

         **read-only**

         Solution number.

         :type: str
      )doc" );
    py::class_< tio::SinexCivilUtc >( m, "SinexCivilUtc", R"doc(

Calendar UTC epoch of a SINEX record (year, day of year and seconds of day). No leap-second table is used, so a seconds-of-day value of 86400 is retained as given.

      )doc" )
            .def_readonly( "year", &tio::SinexCivilUtc::year, R"doc(

         **read-only**

         Four-digit year (two-digit years 00-50 map to 2000-2050, and 51-99 to 1951-1999).

         :type: int
      )doc" )
            .def_readonly( "day_of_year", &tio::SinexCivilUtc::dayOfYear, R"doc(

         **read-only**

         Day of year.

         :type: int
      )doc" )
            .def_readonly( "seconds_of_day", &tio::SinexCivilUtc::secondsOfDay, R"doc(

         **read-only**

         Seconds of day.

         :type: int
      )doc" );
    py::class_< tio::SinexEpoch >( m, "SinexEpoch", R"doc(

Epoch of a SINEX record, as given in the file (``YY:DDD:SSSSS``), with its meaning.

      )doc" )
            .def_readonly( "meaning", &tio::SinexEpoch::meaning, R"doc(

         **read-only**

         Meaning of the epoch.

         :type: SinexEpochMeaning
      )doc" )
            .def_readonly( "civil_utc", &tio::SinexEpoch::civilUtc, R"doc(

         **read-only**

         Calendar UTC epoch, or ``None`` for a zero epoch.

         :type: SinexCivilUtc | None
      )doc" );
    py::class_< tio::SinexSiteRecord >( m, "SinexSiteRecord", R"doc(

Record of the ``SITE/ID`` block of a SINEX file.

      )doc" )
            .def_readonly( "key", &tio::SinexSiteRecord::key, R"doc(

         **read-only**

         Site of the record.

         :type: SinexSiteKey
      )doc" )
            .def_readonly( "domes_id", &tio::SinexSiteRecord::domesId, R"doc(

         **read-only**

         DOMES number of the monument.

         :type: str
      )doc" )
            .def_readonly( "observation_technique", &tio::SinexSiteRecord::observationTechnique, R"doc(

         **read-only**

         Observation technique (``C``, ``D``, ``L``, ``M``, ``P`` or ``R``).

         :type: str
      )doc" )
            .def_readonly( "description", &tio::SinexSiteRecord::description, R"doc(

         **read-only**

         Free-text description of the site.

         :type: str
      )doc" )
            .def_readonly( "approximate_longitude_radians", &tio::SinexSiteRecord::approximateLongitudeRadians, R"doc(

         **read-only**

         Approximate longitude in radians, or ``None``.

         :type: float | None
      )doc" )
            .def_readonly( "approximate_latitude_radians", &tio::SinexSiteRecord::approximateLatitudeRadians, R"doc(

         **read-only**

         Approximate latitude in radians, or ``None``.

         :type: float | None
      )doc" )
            .def_readonly( "approximate_height_metres", &tio::SinexSiteRecord::approximateHeightMetres, R"doc(

         **read-only**

         Approximate height in metres, or ``None``.

         :type: float | None
      )doc" );
    py::class_< tio::SinexParameterRecord >( m, "SinexParameterRecord", R"doc(

Estimated parameter of the ``SOLUTION/ESTIMATE`` block of a SINEX file (for instance a station coordinate, velocity, post-seismic deformation or seasonal term).

      )doc" )
            .def_readonly( "index", &tio::SinexParameterRecord::index, R"doc(

         **read-only**

         Parameter index in the file.

         :type: int
      )doc" )
            .def_readonly( "type", &tio::SinexParameterRecord::type, R"doc(

         **read-only**

         Parameter type (for instance ``STAX``, ``VELX``, ``AEXP_N``).

         :type: str
      )doc" )
            .def_readonly( "key", &tio::SinexParameterRecord::key, R"doc(

         **read-only**

         Solution to which the parameter belongs.

         :type: SinexSolutionKey
      )doc" )
            .def_readonly( "reference_epoch", &tio::SinexParameterRecord::referenceEpoch, R"doc(

         **read-only**

         Reference epoch of the parameter.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "original_unit", &tio::SinexParameterRecord::originalUnit, R"doc(

         **read-only**

         Unit of the parameter in the file.

         :type: str
      )doc" )
            .def_readonly( "constraint_code", &tio::SinexParameterRecord::constraintCode, R"doc(

         **read-only**

         Constraint code of the parameter.

         :type: str
      )doc" )
            .def_readonly( "original_value", &tio::SinexParameterRecord::originalValue, R"doc(

         **read-only**

         Value of the parameter in the file's unit.

         :type: float
      )doc" )
            .def_readonly( "original_standard_deviation", &tio::SinexParameterRecord::originalStandardDeviation, R"doc(

         **read-only**

         Standard deviation of the parameter in the file's unit.

         :type: float
      )doc" )
            .def_readonly( "value_si", &tio::SinexParameterRecord::valueSi, R"doc(

         **read-only**

         Value in SI units (velocities use a 365.25-day year), or ``None`` if not converted.

         :type: float | None
      )doc" )
            .def_readonly( "standard_deviation_si", &tio::SinexParameterRecord::standardDeviationSi, R"doc(

         **read-only**

         Standard deviation in SI units, or ``None`` if not converted.

         :type: float | None
      )doc" )
            .def_readonly( "si_unit", &tio::SinexParameterRecord::siUnit, R"doc(

         **read-only**

         SI unit of the converted value.

         :type: str
      )doc" );
    py::class_< tio::SinexObservationSpan >( m, "SinexObservationSpan", R"doc(

Record of the ``SOLUTION/EPOCHS`` block of a SINEX file: the span of the observations used for a solution.

      )doc" )
            .def_readonly( "key", &tio::SinexObservationSpan::key, R"doc(

         **read-only**

         Solution of the record.

         :type: SinexSolutionKey
      )doc" )
            .def_readonly( "observation_technique", &tio::SinexObservationSpan::observationTechnique, R"doc(

         **read-only**

         Observation technique.

         :type: str
      )doc" )
            .def_readonly( "start", &tio::SinexObservationSpan::start, R"doc(

         **read-only**

         First observation epoch.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "end", &tio::SinexObservationSpan::end, R"doc(

         **read-only**

         Last observation epoch.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "mean", &tio::SinexObservationSpan::mean, R"doc(

         **read-only**

         Mean observation epoch.

         :type: SinexEpoch
      )doc" );
    py::class_< tio::SinexDiscontinuityRecord >( m, "SinexDiscontinuityRecord", R"doc(

Record of the ``SOLUTION/DISCONTINUITY`` block (ITRF extension) of a SINEX file: the validity interval of a position or velocity solution.

      )doc" )
            .def_readonly( "key", &tio::SinexDiscontinuityRecord::key, R"doc(

         **read-only**

         Solution of the record.

         :type: SinexSolutionKey
      )doc" )
            .def_readonly( "observation_technique", &tio::SinexDiscontinuityRecord::observationTechnique, R"doc(

         **read-only**

         Observation technique.

         :type: str
      )doc" )
            .def_readonly( "kind", &tio::SinexDiscontinuityRecord::kind, R"doc(

         **read-only**

         Kind of discontinuity (``P`` for position, ``V`` for velocity).

         :type: str
      )doc" )
            .def_readonly( "start", &tio::SinexDiscontinuityRecord::start, R"doc(

         **read-only**

         Start of the validity interval.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "end", &tio::SinexDiscontinuityRecord::end, R"doc(

         **read-only**

         End of the validity interval.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "description", &tio::SinexDiscontinuityRecord::description, R"doc(

         **read-only**

         Description of the discontinuity (for instance an earthquake or an equipment change).

         :type: str
      )doc" );
    py::class_< tio::SinexEccentricityRecord >( m, "SinexEccentricityRecord", R"doc(

Record of the ``SITE/ECCENTRICITY`` block of a SINEX file: the offset of the antenna reference point from the monument marker.

      )doc" )
            .def_readonly( "key", &tio::SinexEccentricityRecord::key, R"doc(

         **read-only**

         Solution of the record.

         :type: SinexSolutionKey
      )doc" )
            .def_readonly( "observation_technique", &tio::SinexEccentricityRecord::observationTechnique, R"doc(

         **read-only**

         Observation technique.

         :type: str
      )doc" )
            .def_readonly( "start", &tio::SinexEccentricityRecord::start, R"doc(

         **read-only**

         Start of the validity interval.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "end", &tio::SinexEccentricityRecord::end, R"doc(

         **read-only**

         End of the validity interval.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "basis", &tio::SinexEccentricityRecord::basis, R"doc(

         **read-only**

         Component order of the eccentricity.

         :type: SinexEccentricityBasis
      )doc" )
            .def_readonly( "components_metres", &tio::SinexEccentricityRecord::componentsMetres, R"doc(

         **read-only**

         Eccentricity in metres, in the component order of :attr:`basis`.

         :type: list[float]
      )doc" );
    py::class_< tio::SinexEquipmentRecord >( m, "SinexEquipmentRecord", R"doc(

Record of the ``SITE/RECEIVER`` or ``SITE/ANTENNA`` block of a SINEX file.

      )doc" )
            .def_readonly( "kind", &tio::SinexEquipmentRecord::kind, R"doc(

         **read-only**

         Type of equipment.

         :type: SinexEquipmentKind
      )doc" )
            .def_readonly( "key", &tio::SinexEquipmentRecord::key, R"doc(

         **read-only**

         Solution of the record.

         :type: SinexSolutionKey
      )doc" )
            .def_readonly( "observation_technique", &tio::SinexEquipmentRecord::observationTechnique, R"doc(

         **read-only**

         Observation technique.

         :type: str
      )doc" )
            .def_readonly( "start", &tio::SinexEquipmentRecord::start, R"doc(

         **read-only**

         Start of the validity interval.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "end", &tio::SinexEquipmentRecord::end, R"doc(

         **read-only**

         End of the validity interval.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "equipment_type", &tio::SinexEquipmentRecord::equipmentType, R"doc(

         **read-only**

         Equipment model (for antennas, including the radome code).

         :type: str
      )doc" )
            .def_readonly( "serial_number", &tio::SinexEquipmentRecord::serialNumber, R"doc(

         **read-only**

         Serial number.

         :type: str
      )doc" )
            .def_readonly( "firmware_version", &tio::SinexEquipmentRecord::firmwareVersion, R"doc(

         **read-only**

         Firmware version (receivers only).

         :type: str
      )doc" );
    py::class_< tio::SinexCartesianSolution >( m, "SinexCartesianSolution", R"doc(

Cartesian position and (optional) velocity of a single solution of a SINEX file, as returned by :func:`get_station_coordinate_solution`.

      )doc" )
            .def_readonly( "key", &tio::SinexCartesianSolution::key, R"doc(

         **read-only**

         Solution.

         :type: SinexSolutionKey
      )doc" )
            .def_readonly( "site", &tio::SinexCartesianSolution::site, R"doc(

         **read-only**

         Site record of the solution.

         :type: SinexSiteRecord
      )doc" )
            .def_readonly( "position_reference_epoch", &tio::SinexCartesianSolution::positionReferenceEpoch, R"doc(

         **read-only**

         Reference epoch of the position.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "position_metres", &tio::SinexCartesianSolution::positionMetres, R"doc(

         **read-only**

         Cartesian position in metres.

         :type: list[float]
      )doc" )
            .def_readonly( "position_parameter_indices", &tio::SinexCartesianSolution::positionParameterIndices, R"doc(

         **read-only**

         Parameter indices of the position components.

         :type: list[int]
      )doc" )
            .def_readonly( "position_standard_deviations_metres", &tio::SinexCartesianSolution::positionStandardDeviationsMetres, R"doc(

         **read-only**

         Standard deviations of the position components in metres.

         :type: list[float]
      )doc" )
            .def_readonly( "velocity_metres_per_second", &tio::SinexCartesianSolution::velocityMetresPerSecond, R"doc(

         **read-only**

         Cartesian velocity components in metres per second (``None`` for components that are not given).

         :type: list[float | None]
      )doc" )
            .def_readonly( "velocity_parameter_indices", &tio::SinexCartesianSolution::velocityParameterIndices, R"doc(

         **read-only**

         Parameter indices of the velocity components (``None`` if not given).

         :type: list[int | None]
      )doc" )
            .def_readonly( "velocity_standard_deviations_metres_per_second",
                           &tio::SinexCartesianSolution::velocityStandardDeviationsMetresPerSecond,
                           R"doc(

         **read-only**

         Standard deviations of the velocity components in metres per second (``None`` if not given).

         :type: list[float | None]
      )doc" )
            .def_readonly( "velocity_reference_epochs", &tio::SinexCartesianSolution::velocityReferenceEpochs, R"doc(

         **read-only**

         Reference epochs of the velocity components (``None`` if not given).

         :type: list[SinexEpoch | None]
      )doc" );
    py::class_< tio::SinexStationHistory >( m, "SinexStationHistory", R"doc(

Contents of a SINEX station file (for instance an ITRF or IGS station coordinate solution, a post-seismic
deformation or seasonal-term file), as returned by :func:`read_sinex_station_history`.

The records are stored as read, without combining the different blocks or applying any model.

      )doc" )
            .def_readonly( "source", &tio::SinexStationHistory::source, R"doc(

         **read-only**

         Name of the file.

         :type: str
      )doc" )
            .def_readonly( "version", &tio::SinexStationHistory::version, R"doc(

         **read-only**

         SINEX format version.

         :type: str
      )doc" )
            .def_readonly( "agency", &tio::SinexStationHistory::agency, R"doc(

         **read-only**

         Agency that created the file.

         :type: str
      )doc" )
            .def_readonly( "creation_epoch", &tio::SinexStationHistory::creationEpoch, R"doc(

         **read-only**

         Creation epoch of the file.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "header_start", &tio::SinexStationHistory::headerStart, R"doc(

         **read-only**

         Start epoch of the data, from the header.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "header_end", &tio::SinexStationHistory::headerEnd, R"doc(

         **read-only**

         End epoch of the data, from the header.

         :type: SinexEpoch
      )doc" )
            .def_readonly( "observation_technique", &tio::SinexStationHistory::observationTechnique, R"doc(

         **read-only**

         Observation technique(s) of the file.

         :type: str
      )doc" )
            .def_readonly( "declared_parameter_count", &tio::SinexStationHistory::declaredParameterCount, R"doc(

         **read-only**

         Number of estimated parameters declared in the header.

         :type: int
      )doc" )
            .def_readwrite( "frame", &tio::SinexStationHistory::frame, R"doc(

         Terrestrial reference frame of the solutions. Not read from the file; may be set by the user.

         :type: str
      )doc" )
            .def_readwrite( "release", &tio::SinexStationHistory::release, R"doc(

         Release of the solutions. Not read from the file; may be set by the user.

         :type: str
      )doc" )
            .def_readwrite( "origin", &tio::SinexStationHistory::origin, R"doc(

         Reference origin of the solutions. Not read from the file; may be set by the user.

         :type: str
      )doc" )
            .def_readwrite( "tide_system", &tio::SinexStationHistory::tideSystem, R"doc(

         Permanent tide system of the solutions. Not read from the file; may be set by the user.

         :type: str
      )doc" )
            .def_readonly( "parameters", &tio::SinexStationHistory::parameters, R"doc(

         **read-only**

         Estimated parameters, with the parameter index as key.

         :type: dict[int, SinexParameterRecord]
      )doc" )
            .def_readonly( "observation_spans", &tio::SinexStationHistory::observationSpans, R"doc(

         **read-only**

         Observation spans of the solutions.

         :type: list[SinexObservationSpan]
      )doc" )
            .def_readonly( "discontinuities", &tio::SinexStationHistory::discontinuities, R"doc(

         **read-only**

         Validity intervals of the position and velocity solutions.

         :type: list[SinexDiscontinuityRecord]
      )doc" )
            .def_readonly( "eccentricities", &tio::SinexStationHistory::eccentricities, R"doc(

         **read-only**

         Antenna eccentricities.

         :type: list[SinexEccentricityRecord]
      )doc" )
            .def_readonly( "equipment_histories", &tio::SinexStationHistory::equipmentHistories, R"doc(

         **read-only**

         Receiver and antenna histories.

         :type: list[SinexEquipmentRecord]
      )doc" )
            .def_property_readonly(
                    "sites",
                    []( const tio::SinexStationHistory& history ) {
                        // A Python dict would silently erase duplicate site/point records.
                        std::vector< tio::SinexSiteRecord > records;
                        for( const auto& row : history.sites )
                        {
                            records.push_back( row.second );
                        }
                        return records;
                    },
                    R"doc(

         **read-only**

         Records of the ``SITE/ID`` block. In combined products, the same site and point code can be used for
         different monuments (techniques), so all records are retained.

         :type: list[SinexSiteRecord]
      )doc" );

    m.def( "read_sinex_station_history",
           py::overload_cast< const std::string& >( &tio::readSinexStationHistory ),
           py::arg( "file_name" ),
           R"doc(

Function to read a SINEX station file.

Function to read a SINEX (versions 2.00 to 2.02) station file, such as the ITRF and IGS station coordinate solutions,
including the IGN extensions (discontinuities, post-seismic deformation and seasonal-term parameters). All blocks are
read into records; unknown blocks are retained without interpretation.

Values are converted to SI units where the unit is defined (positions in metres, velocities in metres per second
using a 365.25-day year). Epochs are retained as calendar UTC, without a leap-second table. Malformed fields,
non-finite values, reversed intervals, duplicated components, inconsistent parameter counts and truncated blocks raise
an exception that gives the file and line.


Parameters
----------
file_name : str
    Path of the file that is to be read.

Returns
-------
SinexStationHistory
    Contents of the file.


Examples
--------
.. code-block:: python

   from tudatpy.data_input.environment_data import station_history

   history = station_history.read_sinex_station_history( "ITRF2020-IGS-TRF.SSC" )
   key = station_history.SinexSolutionKey( "AB07", "A", "1" )
   solution = station_history.get_station_coordinate_solution( history, key )

      )doc" );
    m.def( "get_station_coordinate_solution_keys",
           &tio::getSinexStationCoordinateSolutionKeys,
           py::arg( "history" ),
           R"doc(

Function to retrieve the identifiers of the station coordinate solutions in a SINEX station file: the solutions (one
per station and validity interval) that have a complete Cartesian position.

Parameters
----------
history : SinexStationHistory
    Contents of the SINEX file.

Returns
-------
list[SinexSolutionKey]
    Solutions that have a complete Cartesian position.

      )doc" );
    m.def( "get_station_coordinate_solution",
           &tio::getSinexStationCoordinateSolution,
           py::arg( "history" ),
           py::arg( "key" ),
           py::arg( "require_complete_velocity" ) = false,
           R"doc(

Function to retrieve a single station coordinate solution of a SINEX station file: the Cartesian position and velocity
of one station in one validity interval.

The solution must have all three position components. Velocity components that are not given are ``None``, unless a
complete velocity is required. A solution whose site is ambiguous (different monuments with the same site and point
code) raises an exception.


Parameters
----------
history : SinexStationHistory
    Contents of the SINEX file.
key : SinexSolutionKey
    Solution that is to be retrieved.
require_complete_velocity : bool, default = False
    Whether all three velocity components must be given.

Returns
-------
SinexCartesianSolution
    Position and velocity of the solution.

      )doc" );
}
}  // namespace station_history
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
