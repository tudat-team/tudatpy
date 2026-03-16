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

#include "coma/expose_coma.h"
#include "ilrs/expose_ilrs.h"
#include "space_weather/expose_space_weather.h"
#include "interface/spice/expose_spice.h"
#include "missions/expose_missions.h"
#include "tudat/io/readSp3File.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{

void expose_environment_data( py::module& m )
{
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

    py::class_< tudat::input_output::Sp3FileContents, std::shared_ptr< tudat::input_output::Sp3FileContents > >(
            m, "Sp3FileContents", R"doc(Parsed contents and metadata of an SP3 precise-orbit file.)doc" )
            .def_property_readonly( "start_epoch",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.startEpoch; } )
            .def_property_readonly( "declared_number_of_epochs",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredNumberOfEpochs; } )
            .def_property_readonly( "declared_epoch_interval",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.declaredEpochInterval; } )
            .def_property_readonly( "frame_name",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.frameName; } )
            .def_property_readonly( "analysis_center",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.analysisCenter; } )
            .def_property_readonly( "time_scale",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.timeScale; } )
            .def_property_readonly( "satellite_states",
                                    []( const tudat::input_output::Sp3FileContents& contents ) { return contents.satelliteStates; } );

    m.def( "read_sp3_file",
           &tudat::input_output::readSp3File,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read an SP3 precise-orbit file.

Epochs are expressed in seconds since ``reference_julian_day`` in the time
scale declared by the file. Positions and velocities are converted to metres
and metres per second, respectively.
)doc" );

    m.def( "read_sp3c_file",
           &tudat::input_output::readSp3cFile,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(Deprecated spelling retained as an alias of :func:`read_sp3_file`.)doc" );
}

}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
