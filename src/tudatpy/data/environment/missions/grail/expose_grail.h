//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATPY_EXPOSE_GRAIL_H
#define TUDATPY_EXPOSE_GRAIL_H

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
namespace grail
{

void expose_grail( py::module& m );

}  // namespace grail
}  // namespace missions
}  // namespace environment
}  // namespace data
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_GRAIL_H
