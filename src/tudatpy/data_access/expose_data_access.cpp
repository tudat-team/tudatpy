#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_data_access.h"

#include "paths/expose_paths.h"
#include "tracking/expose_tracking.h"

namespace tudatpy
{

namespace data_access
{

void expose_data_access( py::module& m )
{
    auto paths = m.def_submodule( "paths" );
    tudatpy::data_access::paths::expose_paths( paths );

    auto tracking = m.def_submodule( "tracking" );
    tudatpy::data_access::tracking::expose_tracking( tracking );
}

}  // namespace data_access

}  // namespace tudatpy
