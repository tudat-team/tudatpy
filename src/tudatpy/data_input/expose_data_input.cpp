#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_data_input.h"

#include "environment_data/expose_environment_data.h"
#include "resource_paths/expose_resource_paths.h"
#include "tracking_data/expose_tracking_data.h"

namespace tudatpy
{

namespace data_input
{

void expose_data_input( py::module& m )
{
    auto resource_paths = m.def_submodule( "resource_paths", R"doc(Resource-path helpers for locating Tudat data files.)doc" );
    tudatpy::data_input::resource_paths::expose_resource_paths( resource_paths );

    auto tracking_data =
            m.def_submodule( "tracking_data", R"doc(Tracking-data readers, containers, and low-level parsed file interfaces.)doc" );
    tudatpy::data_input::tracking_data::expose_tracking_data( tracking_data );

    auto environment_data = m.def_submodule( "environment_data", R"doc(Environment-data readers and query interfaces.)doc" );
    tudatpy::data_input::environment_data::expose_environment_data( environment_data );
}

}  // namespace data_input

}  // namespace tudatpy
