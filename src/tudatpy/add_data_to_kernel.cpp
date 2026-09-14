#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/io/readCrdFile.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/simulation/environment_setup/defaultGroundStationSettings.h"

namespace py = pybind11;
namespace tio = tudat::input_output;
namespace tss = tudat::simulation_setup;

void add_data_to_kernel( py::module_& m )
{
    m.def( "read_crd_file", &tio::readCrdFile, py::arg( "file_name" ) );

    m.def( "read_tracking_txt_file",
           &tio::createTrackingTxtFileContents,
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           py::arg( "comment_symbol" ),
           py::arg( "value_separators" ),
           py::arg( "ignore_omitted_columns" ),
           py::arg( "data_filter_method" ) );

    m.def( "read_ifms_file",
           &tio::readIfmsFile,
           py::arg( "file_name" ),
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true );

    const std::vector< std::string > defaultFdetsColumnTypes = {
        "utc_datetime_string", "signal_to_noise_ratio", "normalised_spectral_max", "doppler_measured_frequency_hz", "doppler_noise_hz"
    };

    m.def( "read_fdets_file",
           py::overload_cast< const std::string&, const std::vector< std::string >& >( &tio::readFdetsFile ),
           py::arg( "file_name" ),
           py::arg( "column_types" ) = defaultFdetsColumnTypes );

    m.def( "set_dsn_weather_data_in_ground_stations",
           py::overload_cast< tss::SystemOfBodies&,
                              const std::vector< std::string >&,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::map< int, std::vector< std::string > >&,
                              const std::string& >( &tio::setDsnWeatherDataInGroundStations ),
           py::arg( "bodies" ),
           py::arg( "weather_file_names" ),
           py::arg_v( "interpolator_settings", tudat::interpolators::linearInterpolation( ), "..." ),
           py::arg_v( "ground_stations_per_complex",
                      tss::getDefaultDsnStationNamesPerComplex( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_default_dsn_station_names_per_complex()" ),
           py::arg( "body_with_ground_stations_name" ) = "Earth" );

    m.def( "set_estrack_weather_data_in_ground_stations",
           py::overload_cast< tss::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::string,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::string& >( &tio::setEstrackWeatherDataInGroundStation ),
           py::arg( "bodies" ),
           py::arg( "weather_file_names" ),
           py::arg( "ground_station_name" ),
           py::arg_v( "interpolator_settings",
                      tudat::interpolators::cubicSplineInterpolation(
                              tudat::interpolators::AvailableLookupScheme::huntingAlgorithm,
                              tudat::interpolators::BoundaryInterpolationType::use_boundary_value_with_warning ),
                      "..." ),
           py::arg( "body_with_ground_stations_name" ) = "Earth" );
}
