//
// Created by Filippo Oggionni on 12/07/21.
//

#ifndef TUDATPY_EXPOSE_TORQUE_SETUP_H
#define TUDATPY_EXPOSE_TORQUE_SETUP_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace propagation_setup
{
namespace torque
{

void expose_torque_setup( py::module& m );

}  // namespace torque
}  // namespace propagation_setup
}  // namespace dynamics
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_TORQUE_SETUP_H
