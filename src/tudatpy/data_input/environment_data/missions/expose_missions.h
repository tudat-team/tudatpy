//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATPY_EXPOSE_MISSIONS_H
#define TUDATPY_EXPOSE_MISSIONS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{
namespace environment_data
{
namespace missions
{

void expose_missions( py::module& m );

}  // namespace missions
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_MISSIONS_H
