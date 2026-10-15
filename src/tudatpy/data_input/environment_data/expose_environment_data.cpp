/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_environment_data.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "tudat/io/readSp3File.h"

#include "coma/expose_coma.h"
#include "earth_orientation/expose_earth_orientation.h"
#include "gravity_field/expose_gravity_field.h"
#include "atmospheric_loading/expose_atmospheric_loading.h"
#include "station_history/expose_station_history.h"
#include "station_loading/expose_station_loading.h"
#include "ocean_loading/expose_ocean_loading.h"
#include "terrestrial_frames/expose_terrestrial_frames.h"
#include "ilrs/expose_ilrs.h"
#include "space_weather/expose_space_weather.h"
#include "interface/spice/expose_spice.h"
#include "missions/expose_missions.h"

namespace py = pybind11;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{

//! Register neutral environmental products and their existing provider modules.
void expose_environment_data( py::module& m )
{

    // SP3 is a raw ephemeris product; its epochs remain in the declared file time system.
    py::class_< tudat::input_output::Sp3Record >(
            m,
            "Sp3Record",
            "One raw SP3 sample. Position and velocity use metres and metres/second; epoch is a Time in the file time coordinate." )
            .def( py::init<>( ) )
            .def_readwrite( "epoch", &tudat::input_output::Sp3Record::epoch )
            .def_readwrite( "source_line", &tudat::input_output::Sp3Record::sourceLine )
            .def_readwrite( "position", &tudat::input_output::Sp3Record::position )
            .def_readwrite( "velocity", &tudat::input_output::Sp3Record::velocity )
            .def_readwrite( "position_available", &tudat::input_output::Sp3Record::positionAvailable )
            .def_readwrite( "velocity_available", &tudat::input_output::Sp3Record::velocityAvailable )
            .def_readwrite( "clock_correction", &tudat::input_output::Sp3Record::clockCorrection )
            .def_readwrite( "clock_rate", &tudat::input_output::Sp3Record::clockRate )
            .def_readwrite( "clock_event", &tudat::input_output::Sp3Record::clockEvent )
            .def_readwrite( "clock_predicted", &tudat::input_output::Sp3Record::clockPredicted )
            .def_readwrite( "maneuver", &tudat::input_output::Sp3Record::maneuver )
            .def_readwrite( "orbit_predicted", &tudat::input_output::Sp3Record::orbitPredicted );

    py::class_< tudat::input_output::Sp3Data >( m,
                                                "Sp3Data",
                                                "SP3c/d source header and ordered satellite records. Frame/time labels and provider "
                                                "comments are retained; no physical conversion or interpolation is inferred." )
            .def( py::init<>( ) )
            .def_readwrite( "source", &tudat::input_output::Sp3Data::source )
            .def_readwrite( "version", &tudat::input_output::Sp3Data::version )
            .def_readwrite( "has_velocity_records", &tudat::input_output::Sp3Data::hasVelocityRecords )
            .def_readwrite( "declared_epoch_count", &tudat::input_output::Sp3Data::declaredEpochCount )
            .def_readwrite( "epoch_interval", &tudat::input_output::Sp3Data::epochInterval )
            .def_readwrite( "data_used", &tudat::input_output::Sp3Data::dataUsed )
            .def_readwrite( "coordinate_system", &tudat::input_output::Sp3Data::coordinateSystem )
            .def_readwrite( "orbit_type", &tudat::input_output::Sp3Data::orbitType )
            .def_readwrite( "agency", &tudat::input_output::Sp3Data::agency )
            .def_readwrite( "time_system", &tudat::input_output::Sp3Data::timeSystem )
            .def_readwrite( "comments", &tudat::input_output::Sp3Data::comments )
            .def_readwrite( "satellite_ids", &tudat::input_output::Sp3Data::satelliteIds )
            .def_readwrite( "records", &tudat::input_output::Sp3Data::records );

    m.def( "read_sp3_file",
           py::overload_cast< const std::string& >( &tudat::input_output::readSp3File ),
           py::arg( "file_name" ),
           "Read an uncompressed SP3c/d ephemeris, preserving frame/time metadata and missing/predicted flags. Units are converted to SI "
           "once; no orbit interpolation or time conversion is performed." );

    // Register provider gravity conventions before consuming dynamics settings.
    auto gravity_field_submodule = m.def_submodule( "gravity_field" );

    gravity_field::expose_gravity_field( gravity_field_submodule );

    auto station_loading_submodule = m.def_submodule( "station_loading" );
    station_loading::expose_station_loading( station_loading_submodule );
    auto terrestrial_frames_submodule = m.def_submodule( "terrestrial_frames" );
    terrestrial_frames::expose_terrestrial_frames( terrestrial_frames_submodule );

    auto atmospheric_loading_submodule = m.def_submodule( "atmospheric_loading" );
    atmospheric_loading::expose_atmospheric_loading( atmospheric_loading_submodule );

    auto station_history_submodule = m.def_submodule( "station_history" );
    station_history::expose_station_history( station_history_submodule );

    auto ocean_loading_submodule = m.def_submodule( "ocean_loading" );
    ocean_loading::expose_ocean_loading( ocean_loading_submodule );

    auto earth_orientation_submodule = m.def_submodule( "earth_orientation" );
    earth_orientation::expose_earth_orientation( earth_orientation_submodule );

    auto coma_submodule = m.def_submodule( "coma" );
    coma::expose_coma( coma_submodule );

    auto missions_submodule = m.def_submodule( "missions" );
    missions::expose_missions( missions_submodule );

    auto ilrs_submodule = m.def_submodule( "ilrs" );
    ilrs::expose_ilrs( ilrs_submodule );

    auto space_weather_submodule = m.def_submodule( "space_weather" );
    space_weather::expose_space_weather( space_weather_submodule );

    auto spice_submodule = m.def_submodule( "spice" );
    tudatpy::interface::spice::expose_spice( spice_submodule );
}

}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
