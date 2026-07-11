#ifndef TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_GENERIC_TEXT_FILE_H
#define TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_GENERIC_TEXT_FILE_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace generic_text_file
{

void expose_generic_text_file( py::module& m );

}  // namespace generic_text_file

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_GENERIC_TEXT_FILE_H
