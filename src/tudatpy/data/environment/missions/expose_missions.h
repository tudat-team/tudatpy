//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATPY_EXPOSE_MISSIONS_H
#define TUDATPY_EXPOSE_MISSIONS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data
{
namespace environment
{
namespace missions
{

void expose_missions( py::module& m );

}  // namespace missions
}  // namespace environment
}  // namespace data
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_MISSIONS_H
