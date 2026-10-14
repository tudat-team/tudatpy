/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "expose_ground_station.h"
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/ground_stations/transformedStationHistoryTrajectory.h"

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace ground_station
{
namespace gs = tudat::ground_stations;
namespace ss = tudat::simulation_setup;
namespace io = tudat::input_output;

void expose_station_history_models( py::module& m )
{
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
    py::class_< gs::StationHistoryTrajectory, std::shared_ptr< gs::StationHistoryTrajectory > >( m, "StationHistoryTrajectory" )
            .def( py::init< std::vector< gs::StationHistoryLinearSolution >,
                            std::vector< gs::StationPostSeismicTerm >,
                            std::vector< gs::StationSeasonalTerm >,
                            gs::StationHistoryMetadata >( ),
                  py::arg( "solutions" ),
                  py::arg( "post_seismic_terms" ),
                  py::arg( "seasonal_terms" ),
                  py::arg( "metadata" ) )
            .def( "position", &gs::StationHistoryTrajectory::getPosition, py::arg( "source_time" ) )
            .def( "state", &gs::StationHistoryTrajectory::getState, py::arg( "source_time" ) )
            .def_property_readonly( "metadata", []( const gs::StationHistoryTrajectory& trajectory ) { return trajectory.getMetadata( ); } )
            .def_property_readonly( "solutions", &gs::StationHistoryTrajectory::getSolutions, py::return_value_policy::copy )
            .def_property_readonly(
                    "post_seismic_terms", &gs::StationHistoryTrajectory::getPostSeismicTerms, py::return_value_policy::copy )
            .def_property_readonly( "seasonal_terms", &gs::StationHistoryTrajectory::getSeasonalTerms, py::return_value_policy::copy )
            .def_property_readonly( "has_fitted_seasonal_terms", &gs::StationHistoryTrajectory::hasFittedSeasonalTerms );
    m.def( "ign_nominal_utc_station_time_from_tai",
           &gs::getIgnNominalUtcStationTimeFromTai,
           py::arg( "tai" ),
           py::arg( "tai_rate_with_respect_to_tdb" ) = py::none( ) );
    m.def( "ign_nominal_utc_station_time_from_tdb",
           py::overload_cast< const tudat::Time&, std::shared_ptr< tudat::earth_orientation::TerrestrialTimeScaleConverter >, bool >(
                   &gs::getIgnNominalUtcStationTimeFromTdb ),
           py::arg( "tdb" ),
           py::arg( "converter" ) = nullptr,
           py::arg( "require_rates" ) = false );
    m.def( "ign_nominal_utc_station_time_from_tdb",
           py::overload_cast< double, std::shared_ptr< tudat::earth_orientation::TerrestrialTimeScaleConverter >, bool >(
                   &gs::getIgnNominalUtcStationTimeFromTdb ),
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
    m.def( "ign_station_history",
           &gs::createIgnStationHistoryTrajectory,
           py::arg( "solutions" ),
           py::arg( "discontinuities" ),
           py::arg( "site" ),
           py::arg( "post_seismic" ) = nullptr,
           py::arg( "seasonal" ) = nullptr,
           py::arg( "enu_to_terrestrial" ) = py::none( ),
           py::arg( "domes_selector" ) = py::none( ),
           py::arg( "velocity_consistency_tolerance_metres_per_julian_year" ) = 1.0e-8,
           "Construct the explicit IGN source profile from preloaded station products." );
    py::class_< gs::TransformedStationHistoryTrajectory,
                std::shared_ptr< gs::TransformedStationHistoryTrajectory >,
                gs::StationHistoryTrajectory >( m, "TransformedStationHistoryTrajectory" )
            .def_property_readonly( "source_trajectory", &gs::TransformedStationHistoryTrajectory::getSourceTrajectory )
            .def_property_readonly( "catalogue",
                                    []( const gs::TransformedStationHistoryTrajectory& trajectory ) { return trajectory.getCatalogue( ); } )
            .def_property_readonly( "acknowledged_release_alignment",
                                    &gs::TransformedStationHistoryTrajectory::getAcknowledgedReleaseAlignment );
    m.def( "transform_station_history",
           &gs::createTransformedStationHistoryTrajectory,
           py::arg( "source" ),
           py::arg( "catalogue" ),
           py::arg( "target_frame" ),
           py::arg( "target_origin" ),
           py::arg( "target_tide_system" ),
           py::arg( "acknowledge_aligned_releases" ) = false,
           "Transform the complete absolute trajectory using preloaded affine parameters and shared source time."
           "Origin/tide changes reject; release alignment requires an explicit declaration." );
    py::class_< ss::StationHistoryGroundStationMotionSettings,
                std::shared_ptr< ss::StationHistoryGroundStationMotionSettings >,
                ss::GroundStationMotionSettings >( m, "StationHistoryGroundStationMotionSettings" );
    m.def( "station_history_motion",
           &ss::stationHistoryGroundStationMotionSettings,
           py::arg( "trajectory" ),
           py::arg( "source_time_function" ) = py::none( ),
           "Absolute trajectory motion. Optional callback receives (TDB Time, require_rates)." );
}
}  // namespace ground_station
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
