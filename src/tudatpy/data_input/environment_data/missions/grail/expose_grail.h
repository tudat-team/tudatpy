//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATPY_EXPOSE_GRAIL_H
#define TUDATPY_EXPOSE_GRAIL_H

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
namespace grail
{

void expose_grail( py::module& m );

}  // namespace grail
}  // namespace missions
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_GRAIL_H
