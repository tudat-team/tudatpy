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
    py::class_< tudat::input_output::DsnAmcRecord >( m, "DsnAmcRecord", R"doc(
        Single record of a DSN Advanced Media Calibration (AMC) path delay file. Time is in UTC seconds since J2000,
        zenith path delays in m, azimuth/elevation in deg, cloud liquid in micrometer, wind speed in mph. Quality flags are 0 for
        nominal data.
    )doc" )
            .def_readonly( "utc_time", &tudat::input_output::DsnAmcRecord::utcTime )
            .def_readonly( "azimuth", &tudat::input_output::DsnAmcRecord::azimuth )
            .def_readonly( "elevation", &tudat::input_output::DsnAmcRecord::elevation )
            .def_readonly( "zenith_wet_delay", &tudat::input_output::DsnAmcRecord::zenithWetDelay )
            .def_readonly( "zenith_dry_delay", &tudat::input_output::DsnAmcRecord::zenithDryDelay )
            .def_readonly( "cloud_liquid", &tudat::input_output::DsnAmcRecord::cloudLiquid )
            .def_readonly( "algorithm_id", &tudat::input_output::DsnAmcRecord::algorithmId )
            .def_readonly( "cloud_flag", &tudat::input_output::DsnAmcRecord::cloudFlag )
            .def_readonly( "wvr_flag", &tudat::input_output::DsnAmcRecord::wvrFlag )
            .def_readonly( "mtp_flag", &tudat::input_output::DsnAmcRecord::mtpFlag )
            .def_readonly( "surface_flag", &tudat::input_output::DsnAmcRecord::surfaceFlag )
            .def_readonly( "wind_flag", &tudat::input_output::DsnAmcRecord::windFlag )
            .def_readonly( "path_delay_flag", &tudat::input_output::DsnAmcRecord::pathDelayFlag )
            .def_readonly( "wind_speed", &tudat::input_output::DsnAmcRecord::windSpeed )
            .def_readonly( "wind_direction", &tudat::input_output::DsnAmcRecord::windDirection );

    m.def( "read_dsn_amc_file", &tudat::input_output::readDsnAmcFile, py::arg( "file_name" ), R"doc(
        Read all records of a DSN Advanced Media Calibration (AMC) path delay file, without any quality screening.

        Parameters
        ----------
        file_name : str
            Name of AMC file

        Returns
        -------
        list[DsnAmcRecord]
            Records in file
    )doc" );

    m.def( "get_valid_dsn_amc_zenith_delays", &tudat::input_output::getValidDsnAmcZenithDelays, py::arg( "records" ), R"doc(
        Retrieve the valid zenith path delays from DSN AMC records: those for which all quality flags are 0, and the
        path delays are valid (non-negative, and not equal to the 999 cm placeholder for a missing retrieval).

        Parameters
        ----------
        records : list[DsnAmcRecord]
            AMC records, from :func:`read_dsn_amc_file`

        Returns
        -------
        dict[float, numpy.ndarray]
            Zenith [dry, wet] path delays [m], with UTC time (seconds since J2000) as key, which can be used as input for
            :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.measured_tropospheric_light_time_correction`
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
