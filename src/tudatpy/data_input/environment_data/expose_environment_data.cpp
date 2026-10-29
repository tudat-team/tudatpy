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
        Single record of a DSN Advanced Media Calibration (AMC) path delay file (see DSN document DORS_002_020501).
    )doc" )
            .def_readonly( "utc_time", &tudat::input_output::DsnAmcRecord::utcTime, R"doc(UTC time, in seconds since J2000.)doc" )
            .def_readonly( "azimuth",
                           &tudat::input_output::DsnAmcRecord::azimuth,
                           R"doc(Azimuth of the Advanced Water Vapor Radiometer (AWVR) [deg].)doc" )
            .def_readonly( "elevation", &tudat::input_output::DsnAmcRecord::elevation, R"doc(Elevation of the AWVR [deg].)doc" )
            .def_readonly( "zenith_wet_delay",
                           &tudat::input_output::DsnAmcRecord::zenithWetDelay,
                           R"doc(Zenith wet path delay [m] (9.99 m denotes a missing retrieval).)doc" )
            .def_readonly( "zenith_dry_delay",
                           &tudat::input_output::DsnAmcRecord::zenithDryDelay,
                           R"doc(Zenith dry (static) path delay [m] (9.99 m denotes a missing retrieval).)doc" )
            .def_readonly( "cloud_liquid", &tudat::input_output::DsnAmcRecord::cloudLiquid, R"doc(Cloud liquid [micrometer].)doc" )
            .def_readonly(
                    "algorithm_id",
                    &tudat::input_output::DsnAmcRecord::algorithmId,
                    R"doc(Identifier of retrieval algorithm (e.g. 005/006: statistical, 200: Bayesian), with leading zeros retained.)doc" )
            .def_readonly(
                    "cloud_flag",
                    &tudat::input_output::DsnAmcRecord::cloudFlag,
                    R"doc(Cloud flag: 0: liquid water < 20 micrometer, 1: 20-100 micrometer, 2: 100-200 micrometer, 3: > 200 micrometer.)doc" )
            .def_readonly(
                    "wvr_flag",
                    &tudat::input_output::DsnAmcRecord::wvrFlag,
                    R"doc(AWVR flag: 0: all AWVR observables within acceptable range, 1: AWVR observable(s) outside acceptable range.)doc" )
            .def_readonly(
                    "mtp_flag",
                    &tudat::input_output::DsnAmcRecord::mtpFlag,
                    R"doc(Microwave Temperature Profiler (MTP) flag: 0: all MTP observables within acceptable range, 1: MTP observable(s) outside acceptable range.)doc" )
            .def_readonly(
                    "surface_flag",
                    &tudat::input_output::DsnAmcRecord::surfaceFlag,
                    R"doc(Surface meteorology flag: 0: all surface meteorology observables within acceptable range, 1: observable(s) outside acceptable range.)doc" )
            .def_readonly(
                    "wind_flag",
                    &tudat::input_output::DsnAmcRecord::windFlag,
                    R"doc(Wind flag: 0: wind observables within acceptable range, 1: wind observable(s) outside acceptable range.)doc" )
            .def_readonly(
                    "path_delay_flag",
                    &tudat::input_output::DsnAmcRecord::pathDelayFlag,
                    R"doc(Path delay flag: 0: retrieved zenith path delay within acceptable range, 1: retrieved zenith path delay outside acceptable range, 2: elevation outside prescribed range (no retrieval for Bayesian algorithm; retrieval possibly erroneous due to ground pickup for statistical algorithm).)doc" )
            .def_readonly( "wind_speed", &tudat::input_output::DsnAmcRecord::windSpeed, R"doc(Wind speed [mph] (NaN if not in file).)doc" )
            .def_readonly( "wind_direction",
                           &tudat::input_output::DsnAmcRecord::windDirection,
                           R"doc(Wind direction [deg] (NaN if not in file).)doc" );

    m.def( "read_dsn_amc_file", &tudat::input_output::readDsnAmcFile, py::arg( "file_name" ), R"doc(
        Read all records of a DSN Advanced Media Calibration (AMC) path delay file, without any quality screening. Both 16-column
        files (e.g. Juno) and 14-column files without wind speed and direction (e.g. Cassini) are accepted, where the columns
        include the leading 'P'. For 14-column files, the wind speed and direction are set to NaN.

        Parameters
        ----------
        file_name : str
            Name of AMC file

        Returns
        -------
        list[DsnAmcRecord]
            Records in file
    )doc" );

    m.def( "get_valid_dsn_amc_zenith_delays",
           &tudat::input_output::getValidDsnAmcZenithDelays,
           py::arg( "records" ),
           py::arg( "minimum_rejected_cloud_flag" ) = 1,
           py::arg( "reject_on_elevation_flag" ) = true,
           py::arg( "reject_on_mtp_flag" ) = true,
           py::arg( "reject_on_surface_flag" ) = true,
           py::arg( "reject_on_wind_flag" ) = true,
           R"doc(
        Retrieve the valid zenith path delays from DSN AMC records. Records with an invalid wet path delay or bad quality flags
        are rejected, where the input arguments define part of the rejection criteria. An invalid dry path delay is stored as NaN.

        Parameters
        ----------
        records : list[DsnAmcRecord]
            AMC records, from :func:`read_dsn_amc_file`
        minimum_rejected_cloud_flag : int, default = 1
            Cloud flag value at and above which records are rejected (cloud flag 0: liquid water < 20 micrometer, 1: 20-100
            micrometer, 2: 100-200 micrometer, 3: > 200 micrometer)
        reject_on_elevation_flag : bool, default = True
            Boolean denoting whether records with elevation outside the prescribed range are rejected
        reject_on_mtp_flag : bool, default = True
            Boolean denoting whether records with non-zero Microwave Temperature Profiler (MTP) flag are rejected
        reject_on_surface_flag : bool, default = True
            Boolean denoting whether records with non-zero surface meteorology flag are rejected
        reject_on_wind_flag : bool, default = True
            Boolean denoting whether records with non-zero wind flag are rejected

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
