#include <pybind11/pybind11.h>
#include "trajectory_design/expose_trajectory_design.h"

namespace py = pybind11;

void add_trajectory_design_to_kernel(py::module_& m)
{
    tudatpy::trajectory_design::expose_trajectory_design(m);
}