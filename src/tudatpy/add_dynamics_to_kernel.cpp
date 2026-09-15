#include <pybind11/pybind11.h>
#include "dynamics/expose_dynamics.h"

namespace py = pybind11;

void add_dynamics_types_to_kernel( py::module_& m )
{
    tudatpy::dynamics::expose_dynamics_types( m );
}

void add_dynamics_to_kernel( py::module_& m )
{
    tudatpy::dynamics::expose_dynamics( m );
}
