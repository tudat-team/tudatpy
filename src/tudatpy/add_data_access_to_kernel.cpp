#include <pybind11/pybind11.h>

#include "data_access/expose_data_access.h"

namespace py = pybind11;

void add_data_access_to_kernel( py::module_& m )
{
    tudatpy::data_access::expose_data_access( m );
}
