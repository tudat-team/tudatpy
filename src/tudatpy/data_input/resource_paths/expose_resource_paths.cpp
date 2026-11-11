#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_resource_paths.h"

#include <pybind11/pybind11.h>

#include "tudat/io/basicInputOutput.h"

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{

namespace resource_paths
{

void expose_resource_paths( py::module& m )
{
    m.def( "get_resource_path",
           &tudat::paths::get_resource_path,
           R"doc(
         Return the root directory of the Tudat resources.

         Returns
         -------
         str
             Root directory of the Tudat resources.
      )doc" );

    m.def( "get_tudat_data_path",
           &tudat::paths::get_tudat_data_path,
           R"doc(
         Return the root directory of the Tudat data resources.

         Returns
         -------
         str
             Root directory of the Tudat data resources.
      )doc" );

    m.def( "get_tudat_path",
           &tudat::paths::get_tudat_path,
           R"doc(
         Return the Tudat resource root directory.

         Returns
         -------
         str
             Tudat resource root directory.
      )doc" );

    m.def( "get_default_output_path",
           &tudat::paths::get_default_output_path,
           R"doc(
         Return the default Tudat output directory.

         Returns
         -------
         str
             Default Tudat output directory.
      )doc" );

    m.def( "get_test_data_path",
           &tudat::paths::getTudatTestDataPath,
           R"doc(
         Return the Tudat test-data directory.

         Returns
         -------
         str
             Tudat test-data directory.
      )doc" );

    m.def( "get_ephemeris_path",
           &tudat::paths::getEphemerisDataFilesPath,
           R"doc(
         Return the Tudat ephemeris resource directory.

         Returns
         -------
         str
             Tudat ephemeris resource directory.
      )doc" );

    m.def( "get_earth_deformation_path",
           &tudat::paths::getEarthDeformationDataFilesPath,
           R"doc(
         Return the Tudat Earth-deformation resource directory.

         Returns
         -------
         str
             Tudat Earth-deformation resource directory.
      )doc" );

    m.def( "get_earth_orientation_path",
           &tudat::paths::getEarthOrientationDataFilesPath,
           R"doc(
         Return the Tudat Earth-orientation resource directory.

         Returns
         -------
         str
             Tudat Earth-orientation resource directory.
      )doc" );

    m.def( "get_quadrature_path",
           &tudat::paths::getQuadratureDataPath,
           R"doc(
         Return the Tudat Gaussian-quadrature resource directory.

         Returns
         -------
         str
             Tudat Gaussian-quadrature resource directory.
      )doc" );

    m.def( "get_spice_kernel_path",
           &tudat::paths::getSpiceKernelPath,
           R"doc(
         Return the Tudat SPICE-kernel resource directory.

         Returns
         -------
         str
             Tudat SPICE-kernel resource directory.
      )doc" );

    m.def( "get_atmosphere_tables_path",
           &tudat::paths::getAtmosphereTablesPath,
           R"doc(
         Return the Tudat atmosphere-table resource directory.

         Returns
         -------
         str
             Tudat atmosphere-table resource directory.
      )doc" );

    m.def( "get_gravity_models_path",
           &tudat::paths::getGravityModelsPath,
           R"doc(
         Return the Tudat gravity-model resource directory.

         Returns
         -------
         str
             Tudat gravity-model resource directory.
      )doc" );

    m.def( "get_space_weather_path",
           &tudat::paths::getSpaceWeatherDataPath,
           R"doc(
         Return the Tudat space-weather resource directory.

         Returns
         -------
         str
             Tudat space-weather resource directory.
      )doc" );

    m.def( "get_station_location_path",
           &tudat::paths::getStationLocationDataPath,
           R"doc(
         Return the Tudat station-location resource directory.

         Returns
         -------
         str
             Tudat station-location resource directory.
      )doc" );

    m.def( "get_nequick2_path",
           &tudat::paths::getNeQuick2DataPath,
           R"doc(
         Return the Tudat NeQuick2 resource directory.

         Returns
         -------
         str
             Tudat NeQuick2 resource directory.
      )doc" );
}

}  // namespace resource_paths

}  // namespace data_input

}  // namespace tudatpy
