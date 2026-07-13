#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_paths.h"

#include <pybind11/pybind11.h>

#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readHistoryFromFile.h"

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

namespace paths
{

void expose_paths( py::module& m )
{
    m.def( "get_resource_path", &tudat::paths::get_resource_path, R"doc(Return the root directory of the Tudat resources.)doc" );

    m.def( "get_tudat_data_path", &tudat::paths::get_tudat_data_path, R"doc(Return the root directory of the Tudat resources.)doc" );

    m.def( "get_tudat_path", &tudat::paths::get_tudat_path, R"doc(Return the Tudat resource root directory.)doc" );

    m.def( "get_default_output_path", &tudat::paths::get_default_output_path, R"doc(Return the default Tudat output directory.)doc" );

    m.def( "get_test_data_path", &tudat::paths::getTudatTestDataPath, R"doc(Return the Tudat test-data directory.)doc" );

    m.def( "get_ephemeris_path", &tudat::paths::getEphemerisDataFilesPath, R"doc(Return the Tudat ephemeris resource directory.)doc" );

    m.def( "get_earth_deformation_path",
           &tudat::paths::getEarthDeformationDataFilesPath,
           R"doc(Return the Tudat Earth-deformation resource directory.)doc" );

    m.def( "get_earth_orientation_path",
           &tudat::paths::getEarthOrientationDataFilesPath,
           R"doc(Return the Tudat Earth-orientation resource directory.)doc" );

    m.def( "get_quadrature_path",
           &tudat::paths::getQuadratureDataPath,
           R"doc(Return the Tudat Gaussian-quadrature resource directory.)doc" );

    m.def( "get_spice_kernel_path", &tudat::paths::getSpiceKernelPath, R"doc(Return the Tudat SPICE-kernel resource directory.)doc" );

    m.def( "get_atmosphere_tables_path",
           &tudat::paths::getAtmosphereTablesPath,
           R"doc(Return the Tudat atmosphere-table resource directory.)doc" );

    m.def( "get_gravity_models_path", &tudat::paths::getGravityModelsPath, R"doc(Return the Tudat gravity-model resource directory.)doc" );

    m.def( "get_space_weather_path",
           &tudat::paths::getSpaceWeatherDataPath,
           R"doc(Return the Tudat space-weather resource directory.)doc" );

    m.def( "get_station_location_path",
           &tudat::paths::getStationLocationDataPath,
           R"doc(Return the Tudat station-location resource directory.)doc" );

    m.def( "get_nequick2_path", &tudat::paths::getNeQuick2DataPath, R"doc(Return the Tudat NeQuick2 resource directory.)doc" );

    m.def( "read_vector_history_from_file",
           &tudat::input_output::readVectorHistoryFromFile< double, double >,
           py::arg( "vector_size" ),
           py::arg( "file_name" ),
           R"doc(

 Read a vector history from a file.


 Parameters
 ----------
 vector_size : int
     Size of the vector at each epoch.
 file_name : str
     Name of the file containing the vector history.
 Returns
 -------
 Dict[float, numpy.ndarray]
     Dictionary mapping epochs to the vector at the given epoch.






     )doc" );
}

}  // namespace paths

}  // namespace data_access

}  // namespace tudatpy
