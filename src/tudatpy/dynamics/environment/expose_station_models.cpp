/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "expose_environment.h"
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/ground_stations/stationReferencePoint.h"
#include "tudat/astro/ground_stations/transformedStationHistoryTrajectory.h"

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
    // Source coordinates and records support advanced custom trajectories and retain provider provenance.
    py::enum_< gs::StationReferencePointKind >( m, "StationReferencePointKind" )
            .value( "marker", gs::StationReferencePointKind::marker )
            .value( "antenna_reference_point", gs::StationReferencePointKind::antennaReferencePoint )
            .value( "calibrated_instrument_point", gs::StationReferencePointKind::calibratedInstrumentPoint );
    py::class_< gs::StationHistoryTime >(
            m, "StationHistoryTime", "Split source coordinates and optional derivatives with respect to TDB seconds." )
            .def( py::init<>( ) )
            .def_readwrite( "tai", &gs::StationHistoryTime::tai )
            .def_readwrite( "mjd_day", &gs::StationHistoryTime::mjdDay )
            .def_readwrite( "mjd_fraction", &gs::StationHistoryTime::mjdFraction )
            .def_readwrite( "calendar_year", &gs::StationHistoryTime::calendarYear )
            .def_readwrite( "year_fraction", &gs::StationHistoryTime::yearFraction )
            .def_readwrite( "mjd_rate", &gs::StationHistoryTime::mjdRate )
            .def_readwrite( "decimal_year_rate", &gs::StationHistoryTime::decimalYearRate )
            .def_readwrite( "convention", &gs::StationHistoryTime::convention );
    py::class_< gs::StationHistoryInterval >( m, "StationHistoryInterval" )
            .def( py::init<>( ) )
            .def_readwrite( "start_tai", &gs::StationHistoryInterval::startTai )
            .def_readwrite( "end_tai", &gs::StationHistoryInterval::endTai );
    py::class_< gs::StationHistoryMetadata >( m, "StationHistoryMetadata" )
            .def( py::init<>( ) )
            .def_readwrite( "site", &gs::StationHistoryMetadata::site )
            .def_readwrite( "domes_id", &gs::StationHistoryMetadata::domesId )
            .def_readwrite( "frame", &gs::StationHistoryMetadata::frame )
            .def_readwrite( "release", &gs::StationHistoryMetadata::release )
            .def_readwrite( "origin", &gs::StationHistoryMetadata::origin )
            .def_readwrite( "tide_system", &gs::StationHistoryMetadata::tideSystem )
            .def_readwrite( "time_convention", &gs::StationHistoryMetadata::timeConvention )
            .def_readwrite( "reference_point_kind", &gs::StationHistoryMetadata::referencePointKind )
            .def_readwrite( "reference_point_name", &gs::StationHistoryMetadata::referencePointName )
            .def_readwrite( "sources", &gs::StationHistoryMetadata::sources )
            .def_readwrite( "signals", &gs::StationHistoryMetadata::signals );
    py::class_< gs::StationHistoryLinearSolution >( m, "StationHistoryLinearSolution" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &gs::StationHistoryLinearSolution::source )
            .def_readwrite( "validity", &gs::StationHistoryLinearSolution::validity )
            .def_readwrite( "velocity_solution_id", &gs::StationHistoryLinearSolution::velocitySolutionId )
            .def_readwrite( "velocity_validity", &gs::StationHistoryLinearSolution::velocityValidity )
            .def_readwrite( "reference_time", &gs::StationHistoryLinearSolution::referenceTime )
            .def_readwrite( "position_metres", &gs::StationHistoryLinearSolution::positionMetres )
            .def_readwrite( "velocity_metres_per_julian_year", &gs::StationHistoryLinearSolution::velocityMetresPerJulianYear );
    py::enum_< gs::StationPostSeismicFunction >( m, "StationPostSeismicFunction" )
            .value( "logarithmic", gs::StationPostSeismicFunction::logarithmic )
            .value( "exponential", gs::StationPostSeismicFunction::exponential );
    py::class_< gs::StationPostSeismicTerm >( m, "StationPostSeismicTerm" )
            .def( py::init<>( ) )
            .def_readwrite( "function", &gs::StationPostSeismicTerm::function )
            .def_readwrite( "event_time", &gs::StationPostSeismicTerm::eventTime )
            .def_readwrite( "amplitude_metres", &gs::StationPostSeismicTerm::amplitudeMetres )
            .def_readwrite( "relaxation_julian_years", &gs::StationPostSeismicTerm::relaxationJulianYears )
            .def_readwrite( "amplitude_parameter_index", &gs::StationPostSeismicTerm::amplitudeParameterIndex )
            .def_readwrite( "relaxation_parameter_index", &gs::StationPostSeismicTerm::relaxationParameterIndex )
            .def_readonly( "provenance", &gs::StationPostSeismicTerm::provenance );
    py::class_< gs::StationSeasonalTerm >( m, "StationSeasonalTerm" )
            .def( py::init<>( ) )
            .def_readwrite( "solution_id", &gs::StationSeasonalTerm::solutionId )
            .def_readwrite( "cycles_per_calendar_year", &gs::StationSeasonalTerm::cyclesPerCalendarYear )
            .def_readwrite( "cosine_metres", &gs::StationSeasonalTerm::cosineMetres )
            .def_readwrite( "sine_metres", &gs::StationSeasonalTerm::sineMetres )
            .def_readwrite( "parameter_indices", &gs::StationSeasonalTerm::parameterIndices );
    // Trajectories expose copied coefficient records so Python mutation cannot invalidate the native model.
    py::class_< gs::StationHistoryTrajectory, std::shared_ptr< gs::StationHistoryTrajectory > >(
            m, "StationHistoryTrajectory", "Absolute terrestrial trajectory with explicit source time and analytic TDB velocity." )
            .def( py::init< std::vector< gs::StationHistoryLinearSolution >,
                            std::vector< gs::StationPostSeismicTerm >,
                            std::vector< gs::StationSeasonalTerm >,
                            gs::StationHistoryMetadata >( ),
                  py::arg( "solutions" ),
                  py::arg( "post_seismic_terms" ),
                  py::arg( "seasonal_terms" ),
                  py::arg( "metadata" ) )
            .def( "position",
                  &gs::StationHistoryTrajectory::getPosition,
                  py::arg( "source_time" ),
                  "Evaluate absolute Cartesian metres in the frame and origin declared by metadata." )
            .def( "state",
                  &gs::StationHistoryTrajectory::getState,
                  py::arg( "source_time" ),
                  "Evaluate Cartesian metres and metres per TDB second; source-coordinate rates must be available." )
            .def_property_readonly( "metadata", []( const gs::StationHistoryTrajectory& trajectory ) { return trajectory.getMetadata( ); } )
            .def_property_readonly( "solutions", &gs::StationHistoryTrajectory::getSolutions, py::return_value_policy::copy )
            .def_property_readonly(
                    "post_seismic_terms", &gs::StationHistoryTrajectory::getPostSeismicTerms, py::return_value_policy::copy )
            .def_property_readonly( "seasonal_terms", &gs::StationHistoryTrajectory::getSeasonalTerms, py::return_value_policy::copy )
            .def_property_readonly( "has_fitted_seasonal_terms", &gs::StationHistoryTrajectory::hasFittedSeasonalTerms );
    // These conversions implement the named source profile, including optional derivatives with respect to TDB.
    m.def( "ign_nominal_utc_station_time_from_tai",
           &gs::getIgnNominalUtcStationTimeFromTai,
           py::arg( "tai" ),
           py::arg( "tai_rate_with_respect_to_tdb" ) = py::none( ) );
    m.def( "ign_nominal_utc_station_time_from_tdb",
           &gs::getIgnNominalUtcStationTimeFromTdb,
           py::arg( "tdb" ),
           py::arg( "converter" ) = nullptr,
           py::arg( "require_rates" ) = false );
    m.def(
            "ign_nominal_utc_station_time",
            []( int year, int day, int seconds ) { return gs::getIgnNominalUtcStationTime( io::SinexCivilUtc{ year, day, seconds } ); },
            py::arg( "year" ),
            py::arg( "day_of_year" ),
            py::arg( "seconds_of_day" ) );
    m.def(
            "sinex_civil_utc_to_tai",
            []( int year, int day, int seconds ) { return gs::convertSinexCivilUtcToTai( io::SinexCivilUtc{ year, day, seconds } ); },
            py::arg( "year" ),
            py::arg( "day_of_year" ),
            py::arg( "seconds_of_day" ) );
    // Frame wrappers retain input identity, catalogue provenance and the explicit release-alignment decision.
    py::class_< gs::TransformedStationHistoryTrajectory,
                std::shared_ptr< gs::TransformedStationHistoryTrajectory >,
                gs::StationHistoryTrajectory >( m, "TransformedStationHistoryTrajectory" )
            .def_property_readonly( "source_trajectory", &gs::TransformedStationHistoryTrajectory::getSourceTrajectory )
            .def_property_readonly( "catalogue",
                                    []( const gs::TransformedStationHistoryTrajectory& trajectory ) { return trajectory.getCatalogue( ); } )
            .def_property_readonly( "acknowledged_release_alignment",
                                    &gs::TransformedStationHistoryTrajectory::getAcknowledgedReleaseAlignment );

    // Offset histories retain equipment identity and explicit local geometry for ordered point composition.
    py::class_< gs::StationReferencePointIdentity >( m, "StationReferencePointIdentity" )
            .def( py::init<>( ) )
            .def_readwrite( "site", &gs::StationReferencePointIdentity::site )
            .def_readwrite( "domes_id", &gs::StationReferencePointIdentity::domesId )
            .def_readwrite( "kind", &gs::StationReferencePointIdentity::kind )
            .def_readwrite( "name", &gs::StationReferencePointIdentity::name );
    py::enum_< gs::StationLocalVertical >( m, "StationLocalVertical" )
            .value( "geocentric_radial", gs::StationLocalVertical::geocentricRadial )
            .value( "geodetic_ellipsoid_normal", gs::StationLocalVertical::geodeticEllipsoidNormal )
            .value( "explicit_fixed_basis", gs::StationLocalVertical::explicitFixedBasis );
    py::class_< gs::StationReferencePointGeometry >( m, "StationReferencePointGeometry" )
            .def( py::init<>( ) )
            .def_readwrite( "frame", &gs::StationReferencePointGeometry::frame )
            .def_readwrite( "vertical", &gs::StationReferencePointGeometry::vertical )
            .def_readwrite( "enu_to_terrestrial", &gs::StationReferencePointGeometry::enuToTerrestrial )
            .def_readwrite( "anchor_metres", &gs::StationReferencePointGeometry::anchorMetres )
            .def_readwrite( "equatorial_radius_metres", &gs::StationReferencePointGeometry::equatorialRadiusMetres )
            .def_readwrite( "flattening", &gs::StationReferencePointGeometry::flattening );
    py::class_< gs::SinexEquipmentHeaderTimeContext >( m, "SinexEquipmentHeaderTimeContext" )
            .def( py::init<>( ) )
            .def_readwrite( "start_tai", &gs::SinexEquipmentHeaderTimeContext::startTai )
            .def_readwrite( "end_tai", &gs::SinexEquipmentHeaderTimeContext::endTai )
            .def_readwrite( "source", &gs::SinexEquipmentHeaderTimeContext::source );
    py::enum_< gs::StationReferencePointOffsetBasis >( m, "StationReferencePointOffsetBasis" )
            .value( "une", gs::StationReferencePointOffsetBasis::une )
            .value( "xyz", gs::StationReferencePointOffsetBasis::xyz );
    py::class_< gs::StationReferencePointOffset >( m, "StationReferencePointOffset" )
            .def( py::init<>( ) )
            .def_readwrite( "from_point", &gs::StationReferencePointOffset::from )
            .def_readwrite( "to_point", &gs::StationReferencePointOffset::to )
            .def_readwrite( "validity", &gs::StationReferencePointOffset::validity )
            .def_readwrite( "frame", &gs::StationReferencePointOffset::frame )
            .def_readwrite( "basis", &gs::StationReferencePointOffset::basis )
            .def_readwrite( "components_metres", &gs::StationReferencePointOffset::componentsMetres )
            .def_readwrite( "equipment_type", &gs::StationReferencePointOffset::equipmentType )
            .def_readwrite( "equipment_serial_number", &gs::StationReferencePointOffset::equipmentSerialNumber )
            .def_readwrite( "calibration_id", &gs::StationReferencePointOffset::calibrationId )
            .def_readwrite( "signal", &gs::StationReferencePointOffset::signal )
            .def_readwrite( "source", &gs::StationReferencePointOffset::source )
            .def_readwrite( "position_function", &gs::StationReferencePointOffset::positionFunction )
            .def_readwrite( "state_function", &gs::StationReferencePointOffset::stateFunction )
            .def_readwrite( "direction_dependent", &gs::StationReferencePointOffset::directionDependent )
            .def_readonly( "provenance", &gs::StationReferencePointOffset::provenance )
            .def_readonly( "equipment_provenance", &gs::StationReferencePointOffset::equipmentProvenance )
            .def_readonly( "header_time_context", &gs::StationReferencePointOffset::headerTimeContext );
    py::class_< gs::StationReferencePointOffsetHistory, std::shared_ptr< gs::StationReferencePointOffsetHistory > >(
            m, "StationReferencePointOffsetHistory" )
            .def( py::init< std::vector< gs::StationReferencePointOffset >, std::optional< gs::StationReferencePointGeometry > >( ),
                  py::arg( "offsets" ),
                  py::arg( "geometry" ) = py::none( ) )
            .def( "position_correction", &gs::StationReferencePointOffsetHistory::getPositionCorrection, py::arg( "source_time" ) )
            .def( "state_correction", &gs::StationReferencePointOffsetHistory::getStateCorrection, py::arg( "source_time" ) )
            .def_property_readonly( "offsets", &gs::StationReferencePointOffsetHistory::getOffsets, py::return_value_policy::copy )
            .def_property_readonly( "geometry", []( const gs::StationReferencePointOffsetHistory& h ) { return h.getGeometry( ); } );
    py::class_< gs::StationReferencePointTrajectory, std::shared_ptr< gs::StationReferencePointTrajectory >, gs::StationHistoryTrajectory >(
            m, "StationReferencePointTrajectory" )
            .def_property_readonly( "reference_point_identity",
                                    []( const gs::StationReferencePointTrajectory& p ) { return p.getReferencePointIdentity( ); } )
            .def_property_readonly( "source_trajectory", &gs::StationReferencePointTrajectory::getSourceTrajectory )
            .def_property_readonly( "offset_history", &gs::StationReferencePointTrajectory::getOffsetHistory );
    m.def( "station_reference_point_identity", &gs::getStationReferencePointIdentity, py::arg( "trajectory" ) );
}
}  // namespace environment
}  // namespace dynamics
}  // namespace tudatpy
