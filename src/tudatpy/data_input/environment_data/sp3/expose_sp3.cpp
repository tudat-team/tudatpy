/*    Copyright (c) 2010-2026, Delft University of Technology
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

#include "expose_sp3.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "tudat/io/readSp3File.h"

namespace py = pybind11;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace sp3
{

void expose_sp3( py::module& m )
{
    py::class_< tudat::input_output::Sp3FileContents, std::shared_ptr< tudat::input_output::Sp3FileContents > >(
            m, "Sp3FileContents", R"doc(Parsed contents and metadata of an SP3 precise-orbit file.)doc" )
            .def_property_readonly(
                    "format_version",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return std::string( 1, contents.formatVersion ); } )
            .def_property_readonly( "has_velocity_records",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.hasVelocityRecords; } )
            .def_property_readonly( "velocities_were_derived",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.velocitiesWereDerived; } )
            .def_property_readonly( "start_epoch",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.startEpoch; } )
            .def_property_readonly( "reference_julian_day",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.referenceJulianDay; } )
            .def_property_readonly( "declared_number_of_epochs",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredNumberOfEpochs; } )
            .def_property_readonly(
                    "declared_number_of_satellites",
                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredNumberOfSatellites; } )
            .def_property_readonly( "declared_epoch_interval",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredEpochInterval; } )
            .def_property_readonly( "frame_name",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.frameName; } )
            .def_property_readonly( "analysis_center",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.analysisCenter; } )
            .def_property_readonly( "time_scale",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.timeScale; } )
            .def_property_readonly( "satellite_states",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.satelliteStates; } )
            .def( "get_satellite_state_history",
                  &tudat::input_output::Sp3FileContents::getSatelliteStateHistory,
                  py::arg( "satellite_id" ),
                  R"doc(
Return the state history for one satellite.

States are returned in the reference frame and time system declared by the SP3
file. No frame or time conversion is performed by the data-input module.
)doc" );

    m.def( "read_sp3_file",
           &tudat::input_output::readSp3File,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read an SP3 precise-orbit file.

Epochs are expressed in seconds since ``reference_julian_day`` in the time
scale declared by the file. Positions and velocities are converted to metres
and metres per second, respectively. For position-only files, velocities are
reconstructed with second-order three-point finite differences. SP3-a through
SP3-d and the GPS, GLO, GAL, BDT, TAI, UTC, IRN, and QZS time-system tags are
supported. Reference-frame and time-system tags are retained as metadata; the
loader does not transform either one.
)doc" );
}

}  // namespace sp3
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
