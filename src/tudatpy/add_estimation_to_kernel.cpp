#include <pybind11/pybind11.h>
#include "estimation/expose_estimation.h"

namespace py = pybind11;

void add_estimation_types_to_kernel( py::module_& m )
{
    tudatpy::estimation::expose_estimation_types( m );
}

void add_estimation_to_kernel( py::module_& m )
{
    tudatpy::estimation::expose_estimation( m );
}
