/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "expose_environment.h"
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/ground_stations/stationReferencePoint.h"

namespace tudatpy
{
namespace dynamics
{
namespace environment
{
namespace gs = tudat::ground_stations;
namespace io = tudat::input_output;

//! Register runtime station trajectories, their source coordinates and reference-point histories.
void expose_station_models( py::module& m )
{
    py::enum_< gs::StationReferencePointKind >( m, "StationReferencePointKind", R"doc(

Enumeration of the kinds of point on a station of which a station history gives the position.

      )doc" )
            .value( "marker", gs::StationReferencePointKind::marker, "Monument marker (the point of SINEX solutions)." )
            .value( "antenna_reference_point", gs::StationReferencePointKind::antennaReferencePoint, "Antenna reference point." )
            .value( "calibrated_instrument_point",
                    gs::StationReferencePointKind::calibratedInstrumentPoint,
                    "Instrument point defined by a user-supplied calibration (for instance a phase centre)." );

    py::class_< gs::StationHistoryInterval >( m, "StationHistoryInterval", R"doc(

Validity interval of a station history record, :math:`[t_{start}, t_{end})` in TAI. An end that is ``None`` is
unbounded.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "start_tai", &gs::StationHistoryInterval::startTai, R"doc(

         Start of the interval in TAI, or ``None`` if unbounded.

         :type: astro.time_representation.Time | None
      )doc" )
            .def_readwrite( "end_tai", &gs::StationHistoryInterval::endTai, R"doc(

         End of the interval in TAI (excluded), or ``None`` if unbounded.

         :type: astro.time_representation.Time | None
      )doc" );
    py::class_< gs::StationHistoryMetadata >( m, "StationHistoryMetadata", R"doc(

Description of a station history: the monument, the reference frame and conventions of its coordinates, and the point
of the station of which it gives the position.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "site", &gs::StationHistoryMetadata::site, R"doc(

         Site of the station.

         :type: ~tudatpy.data_input.environment_data.station_history.SinexSiteKey
      )doc" )
            .def_readwrite( "domes_id", &gs::StationHistoryMetadata::domesId, R"doc(

         DOMES number of the monument.

         :type: str
      )doc" )
            .def_readwrite( "frame", &gs::StationHistoryMetadata::frame, R"doc(

         Terrestrial reference frame realization of the coordinates (for instance ``ITRF2020``).

         :type: str
      )doc" )
            .def_readwrite( "release", &gs::StationHistoryMetadata::release, R"doc(

         Release of the solution.

         :type: str
      )doc" )
            .def_readwrite( "origin", &gs::StationHistoryMetadata::origin, R"doc(

         Reference origin of the coordinates.

         :type: str
      )doc" )
            .def_readwrite( "tide_system", &gs::StationHistoryMetadata::tideSystem, R"doc(

         Permanent tide system of the coordinates.

         :type: str
      )doc" )
            .def_readwrite( "reference_point_kind", &gs::StationHistoryMetadata::referencePointKind, R"doc(

         Kind of point of which the history gives the position.

         :type: StationReferencePointKind
      )doc" )
            .def_readwrite( "reference_point_name", &gs::StationHistoryMetadata::referencePointName, R"doc(

         Name of that point.

         :type: str
      )doc" )
            .def_readwrite( "sources", &gs::StationHistoryMetadata::sources, R"doc(

         Names of the source files.

         :type: list[str]
      )doc" )
            .def_readwrite( "signals", &gs::StationHistoryMetadata::signals, R"doc(

         Signals to which a calibrated instrument point applies.

         :type: list[str]
      )doc" );
    py::class_< gs::StationHistoryLinearSolution >( m, "StationHistoryLinearSolution", R"doc(

Linear (position and velocity) solution of a station history, valid on a single interval between discontinuities.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &gs::StationHistoryLinearSolution::source, R"doc(

         SINEX solution from which this solution was created.

         :type: ~tudatpy.data_input.environment_data.station_history.SinexCartesianSolution
      )doc" )
            .def_readwrite( "validity", &gs::StationHistoryLinearSolution::validity, R"doc(

         Validity interval of the position solution.

         :type: StationHistoryInterval
      )doc" )
            .def_readwrite( "velocity_solution_id", &gs::StationHistoryLinearSolution::velocitySolutionId, R"doc(

         Identifier of the velocity solution.

         :type: str
      )doc" )
            .def_readwrite( "velocity_validity", &gs::StationHistoryLinearSolution::velocityValidity, R"doc(

         Validity interval of the velocity solution.

         :type: StationHistoryInterval
      )doc" )
            .def_readwrite( "reference_time", &gs::StationHistoryLinearSolution::referenceTime, R"doc(

         Reference epoch of the position, in TAI.

         :type: astro.time_representation.Time
      )doc" )
            .def_readwrite( "position_metres", &gs::StationHistoryLinearSolution::positionMetres, R"doc(

         Cartesian position at the reference epoch, in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "velocity_metres_per_julian_year", &gs::StationHistoryLinearSolution::velocityMetresPerJulianYear, R"doc(

         Cartesian velocity in metres per Julian year.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" );
    py::enum_< gs::StationPostSeismicFunction >( m, "StationPostSeismicFunction", R"doc(

Enumeration of the functions of time of post-seismic deformation terms (ITRF2014/ITRF2020 models).

      )doc" )
            .value( "logarithmic",
                    gs::StationPostSeismicFunction::logarithmic,
                    R"doc(Logarithmic function :math:`A\log(1+\Delta t/\tau)`.)doc" )
            .value( "exponential",
                    gs::StationPostSeismicFunction::exponential,
                    R"doc(Exponential function :math:`A(1-e^{-\Delta t/\tau})`.)doc" );
    py::class_< gs::StationPostSeismicTerm >( m, "StationPostSeismicTerm", R"doc(

Post-seismic deformation term of a station history, applied from the epoch of an earthquake onwards.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "function", &gs::StationPostSeismicTerm::function, R"doc(

         Function of time of the term.

         :type: StationPostSeismicFunction
      )doc" )
            .def_readwrite( "event_time", &gs::StationPostSeismicTerm::eventTime, R"doc(

         Epoch of the earthquake, in TAI.

         :type: astro.time_representation.Time
      )doc" )
            .def_readwrite( "amplitude_metres", &gs::StationPostSeismicTerm::amplitudeMetres, R"doc(

         Amplitude of the east, north and up components, in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "relaxation_julian_years", &gs::StationPostSeismicTerm::relaxationJulianYears, R"doc(

         Relaxation time :math:`\tau` in Julian years.

         :type: float
      )doc" )
            .def_readwrite( "amplitude_parameter_index", &gs::StationPostSeismicTerm::amplitudeParameterIndex, R"doc(

         SINEX parameter index of the amplitude.

         :type: int
      )doc" )
            .def_readwrite( "relaxation_parameter_index", &gs::StationPostSeismicTerm::relaxationParameterIndex, R"doc(

         SINEX parameter index of the relaxation time.

         :type: int
      )doc" )
            .def_readonly( "provenance", &gs::StationPostSeismicTerm::provenance, R"doc(

         **read-only**

         Line of the SINEX file from which the term was read.

         :type: ~tudatpy.data_input.environment_data.station_history.SinexSourceRecord
      )doc" );
    py::class_< gs::StationSeasonalTerm >( m, "StationSeasonalTerm", R"doc(

Seasonal (annual or semi-annual) term of a station history, fitted for a single position solution. The term is
:math:`\mathbf{c}\cos(2\pi k y)+\mathbf{s}\sin(2\pi k y)`, with :math:`y` the calendar decimal year and :math:`k`
the number of cycles per year.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "solution_id", &gs::StationSeasonalTerm::solutionId, R"doc(

         Identifier of the position solution for which the term was fitted.

         :type: str
      )doc" )
            .def_readwrite( "cycles_per_calendar_year", &gs::StationSeasonalTerm::cyclesPerCalendarYear, R"doc(

         Number of cycles per calendar year :math:`k` (1 for annual, 2 for semi-annual).

         :type: int
      )doc" )
            .def_readwrite( "cosine_metres", &gs::StationSeasonalTerm::cosineMetres, R"doc(

         Cosine coefficients :math:`\mathbf{c}` in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "sine_metres", &gs::StationSeasonalTerm::sineMetres, R"doc(

         Sine coefficients :math:`\mathbf{s}` in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "parameter_indices", &gs::StationSeasonalTerm::parameterIndices, R"doc(

         SINEX parameter indices of the coefficients.

         :type: list[int]
      )doc" );

    // Station coordinate histories expose copied coefficient records so Python mutation cannot invalidate the native model.
    py::class_< gs::StationCoordinateHistory, std::shared_ptr< gs::StationCoordinateHistory > >( m, "StationCoordinateHistory", R"doc(

Station coordinate history: the catalogued position and velocity of a station point as a function of time, from linear
solutions between discontinuities, post-seismic deformation and seasonal terms.

At an epoch where a new solution starts, the new solution is used for both position and velocity. Epochs outside the
validity of the solutions, or in a gap between them, raise an exception. Objects of this class are normally created
with :func:`~tudatpy.dynamics.environment_setup.ground_station.ign_station_history` (and transformed to another frame
with :func:`~tudatpy.dynamics.environment_setup.ground_station.transform_station_history`), and added to a ground
station with :func:`~tudatpy.dynamics.environment_setup.ground_station.station_history_motion`.

      )doc" )
            .def( py::init< std::vector< gs::StationHistoryLinearSolution >,
                            std::vector< gs::StationPostSeismicTerm >,
                            std::vector< gs::StationSeasonalTerm >,
                            gs::StationHistoryMetadata >( ),
                  py::arg( "solutions" ),
                  py::arg( "post_seismic_terms" ),
                  py::arg( "seasonal_terms" ),
                  py::arg( "metadata" ),
                  R"doc(

         Constructor, from user-defined records.

         Parameters
         ----------
         solutions : list[StationHistoryLinearSolution]
             Linear solutions, with consecutive validity intervals.
         post_seismic_terms : list[StationPostSeismicTerm]
             Post-seismic deformation terms.
         seasonal_terms : list[StationSeasonalTerm]
             Seasonal terms.
         metadata : StationHistoryMetadata
             Description of the history.

      )doc" )
            .def( "position",
                  &gs::StationCoordinateHistory::getPosition,
                  py::arg( "tai" ),
                  R"doc(

         Function to compute the position of the station.

         Parameters
         ----------
         tai : astro.time_representation.Time
             Epoch (TAI) at which the position is computed.

         Returns
         -------
         numpy.ndarray[numpy.float64[3, 1]]
             Cartesian position in metres, in the frame given by the metadata.

      )doc" )
            .def( "state",
                  &gs::StationCoordinateHistory::getState,
                  py::arg( "tai" ),
                  R"doc(

         Function to compute the position and velocity of the station.

         Parameters
         ----------
         tai : astro.time_representation.Time
             Epoch (TAI) at which the state is computed.

         Returns
         -------
         numpy.ndarray[numpy.float64[6, 1]]
             Cartesian position and velocity, in metres and metres per second.

      )doc" )
            .def_property_readonly(
                    "metadata", []( const gs::StationCoordinateHistory& history ) { return history.getMetadata( ); }, R"doc(

         **read-only**

         Description of the history.

         :type: StationHistoryMetadata
      )doc" )
            .def_property_readonly( "solutions", &gs::StationCoordinateHistory::getSolutions, py::return_value_policy::copy, R"doc(

         **read-only**

         Linear solutions (a copy).

         :type: list[StationHistoryLinearSolution]
      )doc" )
            .def_property_readonly(
                    "post_seismic_terms", &gs::StationCoordinateHistory::getPostSeismicTerms, py::return_value_policy::copy, R"doc(

         **read-only**

         Post-seismic deformation terms (a copy).

         :type: list[StationPostSeismicTerm]
      )doc" )
            .def_property_readonly( "seasonal_terms", &gs::StationCoordinateHistory::getSeasonalTerms, py::return_value_policy::copy, R"doc(

         **read-only**

         Seasonal terms (a copy).

         :type: list[StationSeasonalTerm]
      )doc" )
            .def_property_readonly( "has_fitted_seasonal_terms", &gs::StationCoordinateHistory::hasFittedSeasonalTerms, R"doc(

         **read-only**

         Whether the history contains seasonal terms (which should then not also be modelled by a loading model).

         :type: bool
      )doc" );

    m.def(
            "sinex_civil_utc_to_tai",
            []( int year, int day, int seconds ) { return gs::convertSinexCivilUtcToTai( io::SinexCivilUtc{ year, day, seconds } ); },
            py::arg( "year" ),
            py::arg( "day_of_year" ),
            py::arg( "seconds_of_day" ),
            R"doc(

Function to convert a calendar UTC epoch (as used in SINEX files) to TAI. An epoch during an inserted leap second
(seconds of day 86400) is supported.


Parameters
----------
year : int
    Four-digit year.
day_of_year : int
    Day of year.
seconds_of_day : int
    Seconds of day.

Returns
-------
astro.time_representation.Time
    Epoch in TAI.

      )doc" );

    py::class_< gs::StationReferencePointIdentity >( m, "StationReferencePointIdentity", R"doc(

Identification of a point on a station: its monument and the kind and name of the point.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "site", &gs::StationReferencePointIdentity::site, R"doc(

         Site of the station.

         :type: ~tudatpy.data_input.environment_data.station_history.SinexSiteKey
      )doc" )
            .def_readwrite( "domes_id", &gs::StationReferencePointIdentity::domesId, R"doc(

         DOMES number of the monument.

         :type: str
      )doc" )
            .def_readwrite( "kind", &gs::StationReferencePointIdentity::kind, R"doc(

         Kind of point.

         :type: StationReferencePointKind
      )doc" )
            .def_readwrite( "name", &gs::StationReferencePointIdentity::name, R"doc(

         Name of the point.

         :type: str
      )doc" );
    py::enum_< gs::StationLocalVertical >( m, "StationLocalVertical", R"doc(

Enumeration of definitions of the local vertical, used to convert local (up, north, east) offsets to Cartesian
coordinates.

      )doc" )
            .value( "geocentric_radial", gs::StationLocalVertical::geocentricRadial, "Geocentric radial direction." )
            .value( "geodetic_ellipsoid_normal", gs::StationLocalVertical::geodeticEllipsoidNormal, "Normal to the reference ellipsoid." )
            .value( "explicit_fixed_basis", gs::StationLocalVertical::explicitFixedBasis, "User-defined fixed basis." );
    py::class_< gs::StationReferencePointGeometry >( m, "StationReferencePointGeometry", R"doc(

Fixed local (east, north, up) basis of a station, used to convert local offsets of its reference points to Cartesian
coordinates. Objects of this class are normally created with
:func:`~tudatpy.dynamics.environment_setup.ground_station.geocentric_reference_point_geometry` or
:func:`~tudatpy.dynamics.environment_setup.ground_station.geodetic_reference_point_geometry`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "frame", &gs::StationReferencePointGeometry::frame, R"doc(

         Terrestrial reference frame of the basis.

         :type: str
      )doc" )
            .def_readwrite( "vertical", &gs::StationReferencePointGeometry::vertical, R"doc(

         Definition of the local vertical.

         :type: StationLocalVertical
      )doc" )
            .def_readwrite( "enu_to_terrestrial", &gs::StationReferencePointGeometry::enuToTerrestrial, R"doc(

         Rotation matrix from the local (east, north, up) basis to body-fixed Cartesian coordinates.

         :type: numpy.ndarray[numpy.float64[3, 3]]
      )doc" )
            .def_readwrite( "anchor_metres", &gs::StationReferencePointGeometry::anchorMetres, R"doc(

         Cartesian position at which the basis is defined, in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "equatorial_radius_metres", &gs::StationReferencePointGeometry::equatorialRadiusMetres, R"doc(

         Equatorial radius of the reference ellipsoid in metres (geodetic vertical only).

         :type: float | None
      )doc" )
            .def_readwrite( "flattening", &gs::StationReferencePointGeometry::flattening, R"doc(

         Flattening of the reference ellipsoid (geodetic vertical only).

         :type: float | None
      )doc" );
    py::class_< gs::SinexEquipmentHeaderTimeContext >( m, "SinexEquipmentHeaderTimeContext", R"doc(

Validity interval to use for SINEX equipment records whose start or end refers to the file header, when the header does
not give that epoch.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "start_tai", &gs::SinexEquipmentHeaderTimeContext::startTai, R"doc(

         Epoch to use for a start that refers to the header, in TAI.

         :type: astro.time_representation.Time | None
      )doc" )
            .def_readwrite( "end_tai", &gs::SinexEquipmentHeaderTimeContext::endTai, R"doc(

         Epoch to use for an end that refers to the header, in TAI.

         :type: astro.time_representation.Time | None
      )doc" )
            .def_readwrite( "source", &gs::SinexEquipmentHeaderTimeContext::source, R"doc(

         Description of where these epochs come from.

         :type: str
      )doc" );
    py::enum_< gs::StationReferencePointOffsetBasis >( m, "StationReferencePointOffsetBasis", R"doc(

Enumeration of the component orders of a reference point offset.

      )doc" )
            .value( "une", gs::StationReferencePointOffsetBasis::une, "Up, north and east components." )
            .value( "xyz", gs::StationReferencePointOffsetBasis::xyz, "Body-fixed Cartesian components." );
    py::class_< gs::StationReferencePointOffset >( m, "StationReferencePointOffset", R"doc(

Offset from one point of a station to another (for instance from the marker to the antenna reference point), valid on
a single interval.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "from_point", &gs::StationReferencePointOffset::from, R"doc(

         Point from which the offset is defined.

         :type: StationReferencePointIdentity
      )doc" )
            .def_readwrite( "to_point", &gs::StationReferencePointOffset::to, R"doc(

         Point to which the offset is defined.

         :type: StationReferencePointIdentity
      )doc" )
            .def_readwrite( "validity", &gs::StationReferencePointOffset::validity, R"doc(

         Validity interval of the offset.

         :type: StationHistoryInterval
      )doc" )
            .def_readwrite( "frame", &gs::StationReferencePointOffset::frame, R"doc(

         Terrestrial reference frame of the offset.

         :type: str
      )doc" )
            .def_readwrite( "basis", &gs::StationReferencePointOffset::basis, R"doc(

         Component order of the offset.

         :type: StationReferencePointOffsetBasis
      )doc" )
            .def_readwrite( "components_metres", &gs::StationReferencePointOffset::componentsMetres, R"doc(

         Offset in metres, in the component order of :attr:`basis`.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "equipment_type", &gs::StationReferencePointOffset::equipmentType, R"doc(

         Type of the equipment to which the offset applies.

         :type: str
      )doc" )
            .def_readwrite( "equipment_serial_number", &gs::StationReferencePointOffset::equipmentSerialNumber, R"doc(

         Serial number of the equipment.

         :type: str
      )doc" )
            .def_readwrite( "calibration_id", &gs::StationReferencePointOffset::calibrationId, R"doc(

         Identifier of the calibration that defines the offset.

         :type: str
      )doc" )
            .def_readwrite( "signal", &gs::StationReferencePointOffset::signal, R"doc(

         Signal to which the offset applies.

         :type: str
      )doc" )
            .def_readwrite( "source", &gs::StationReferencePointOffset::source, R"doc(

         Description of the source of the offset.

         :type: str
      )doc" )
            .def_readwrite( "position_function", &gs::StationReferencePointOffset::positionFunction, R"doc(

         Time-dependent offset (in the component order of :attr:`basis`) as a function of TAI, used instead of
         :attr:`components_metres` if set.

         :type: callable[[astro.time_representation.Time], numpy.ndarray[numpy.float64[3, 1]]]
      )doc" )
            .def_readwrite( "state_function", &gs::StationReferencePointOffset::stateFunction, R"doc(

         Time-dependent offset and its rate (per second), required for state evaluation if :attr:`position_function` is set.

         :type: callable[[astro.time_representation.Time], numpy.ndarray[numpy.float64[6, 1]]]
      )doc" )
            .def_readwrite( "direction_dependent", &gs::StationReferencePointOffset::directionDependent, R"doc(

         Whether the offset depends on the observation direction (which is not modelled).

         :type: bool
      )doc" )
            .def_readonly( "provenance", &gs::StationReferencePointOffset::provenance, R"doc(

         **read-only**

         Line of the SINEX file that gives the offset.

         :type: ~tudatpy.data_input.environment_data.station_history.SinexSourceRecord
      )doc" )
            .def_readonly( "equipment_provenance", &gs::StationReferencePointOffset::equipmentProvenance, R"doc(

         **read-only**

         Line of the SINEX file that gives the equipment.

         :type: ~tudatpy.data_input.environment_data.station_history.SinexSourceRecord
      )doc" )
            .def_readonly( "header_time_context", &gs::StationReferencePointOffset::headerTimeContext, R"doc(

         **read-only**

         Validity interval used for header-referenced epochs, if any.

         :type: SinexEquipmentHeaderTimeContext | None
      )doc" );
    py::class_< gs::StationReferencePointOffsetHistory, std::shared_ptr< gs::StationReferencePointOffsetHistory > >(
            m, "StationReferencePointOffsetHistory", R"doc(

Sequence of offsets from one point of a station to another, with consecutive validity intervals. At an epoch where a
new offset starts, the new offset is used for both position and velocity.

      )doc" )
            .def( py::init< std::vector< gs::StationReferencePointOffset >, std::optional< gs::StationReferencePointGeometry > >( ),
                  py::arg( "offsets" ),
                  py::arg( "geometry" ) = py::none( ),
                  R"doc(

         Constructor.

         Parameters
         ----------
         offsets : list[StationReferencePointOffset]
             Offsets, with consecutive validity intervals and the same start and end points.
         geometry : StationReferencePointGeometry, default = None
             Local basis used to convert offsets in up, north, east components; required for such offsets.

      )doc" )
            .def( "position_correction",
                  &gs::StationReferencePointOffsetHistory::getPositionCorrection,
                  py::arg( "tai" ),
                  R"doc(

         Function to compute the Cartesian offset at a given epoch.

         Parameters
         ----------
         tai : astro.time_representation.Time
             Epoch (TAI) at which the offset is computed.

         Returns
         -------
         numpy.ndarray[numpy.float64[3, 1]]
             Cartesian offset in metres.

      )doc" )
            .def( "state_correction",
                  &gs::StationReferencePointOffsetHistory::getStateCorrection,
                  py::arg( "tai" ),
                  R"doc(

         Function to compute the Cartesian offset and its rate at a given epoch.

         Parameters
         ----------
         tai : astro.time_representation.Time
             Epoch (TAI) at which the offset is computed.

         Returns
         -------
         numpy.ndarray[numpy.float64[6, 1]]
             Cartesian offset and its rate, in metres and metres per second.

      )doc" )
            .def_property_readonly( "offsets", &gs::StationReferencePointOffsetHistory::getOffsets, py::return_value_policy::copy, R"doc(

         **read-only**

         Offsets (a copy).

         :type: list[StationReferencePointOffset]
      )doc" )
            .def_property_readonly( "geometry", []( const gs::StationReferencePointOffsetHistory& h ) { return h.getGeometry( ); }, R"doc(

         **read-only**

         Local basis of the offsets, if any.

         :type: StationReferencePointGeometry | None
      )doc" );
    m.def( "station_reference_point_identity",
           &gs::getStationReferencePointIdentity,
           py::arg( "metadata" ),
           R"doc(

Function to retrieve the identification of the point whose position a station coordinate history describes.

Parameters
----------
metadata : StationHistoryMetadata
    Description of the station coordinate history.

Returns
-------
StationReferencePointIdentity
    Identification of the point.

      )doc" );
}
}  // namespace environment
}  // namespace dynamics
}  // namespace tudatpy
