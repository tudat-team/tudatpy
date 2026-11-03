#ifndef TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_GENERIC_TEXT_FILE_H
#define TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_GENERIC_TEXT_FILE_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

namespace generic_text_file
{

void expose_generic_text_file( py::module& m );

}  // namespace generic_text_file

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_GENERIC_TEXT_FILE_H
