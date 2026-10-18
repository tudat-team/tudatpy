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

#include <pybind11/stl.h>

#include "tudat/io/readDsnAmcFile.h"

#include "coma/expose_coma.h"
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

void expose_environment_data( py::module& m )
{
    using tudat::input_output::DsnAmcData;
    using tudat::input_output::DsnAmcRecord;

    py::class_< DsnAmcRecord >( m, "DsnAmcRecord", R"doc(
        One unfiltered DSN AMC record. Delays are zenith path lengths in metres;
        time is UTC seconds since J2000. Pointing angles are the AWVR values in
        degrees, cloud liquid is in micrometres and wind speed is in mph.
        Unknown quality flags and retrieval identifiers are preserved.
    )doc" )
            .def_readonly( "utc_time", &DsnAmcRecord::utcTime )
            .def_readonly( "azimuth", &DsnAmcRecord::azimuth )
            .def_readonly( "elevation", &DsnAmcRecord::elevation )
            .def_readonly( "zenith_wet_delay", &DsnAmcRecord::zenithWetDelay )
            .def_readonly( "zenith_dry_delay", &DsnAmcRecord::zenithDryDelay )
            .def_readonly( "cloud_liquid", &DsnAmcRecord::cloudLiquid )
            .def_readonly( "algorithm_id", &DsnAmcRecord::algorithmId )
            .def_readonly( "cloud_flag", &DsnAmcRecord::cloudFlag )
            .def_readonly( "wvr_flag", &DsnAmcRecord::wvrFlag )
            .def_readonly( "mtp_flag", &DsnAmcRecord::mtpFlag )
            .def_readonly( "surface_flag", &DsnAmcRecord::surfaceFlag )
            .def_readonly( "wind_flag", &DsnAmcRecord::windFlag )
            .def_readonly( "path_delay_flag", &DsnAmcRecord::pathDelayFlag )
            .def_readonly( "wind_speed", &DsnAmcRecord::windSpeed )
            .def_readonly( "wind_direction", &DsnAmcRecord::windDirection );

    py::class_< DsnAmcData, std::shared_ptr< DsnAmcData > >( m, "DsnAmcData", R"doc(
        Records from one AMC unit and their source filename. The reader does
        not assign the unit to a tracking station or screen its quality flags.
    )doc" )
            .def_property_readonly( "records", &DsnAmcData::getRecords )
            .def_property_readonly( "file_name", &DsnAmcData::getFileName );

    m.def( "read_dsn_amc_file", &tudat::input_output::readDsnAmcFile, py::arg( "file_name" ), R"doc(
        Read a DSN Advanced Media Calibration path-delay file.

        Convert 1970 UT epochs to UTC seconds since J2000 and centimetres of
        zenith wet/dry path length to metres. Preserve all records, quality
        flags, algorithm identifiers and gaps. Malformed rows and duplicate
        or decreasing epochs raise RuntimeError. Station association, quality
        screening and interpolation are left to the correction settings.
    )doc" );

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
