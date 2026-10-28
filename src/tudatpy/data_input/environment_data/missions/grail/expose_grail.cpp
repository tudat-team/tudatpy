/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_grail.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
// #include <pybind11/native_enum.h>
#include <tudat/io/basicInputOutput.h>
#include <tudat/io/readVariousPdsFiles.h>
#include <string>
#include <vector>

namespace py = pybind11;
namespace tio = tudat::input_output;

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

void expose_grail( py::module& m )
{
    m.def( "grail_antenna_file_reader",
           &tio::grailAntennaFileReader,
           py::arg( "file_name" ),
           R"doc(
         Read a GRAIL antenna file.

         Parameters
         ----------
         file_name : str
             Path to the GRAIL antenna file.

         Returns
         -------
         dict
             Parsed GRAIL antenna-file contents.
      )doc" );
    m.def( "grail_mass_level_0_file_reader",
           &tio::grailMassLevel0FileReader,
           py::arg( "file_name" ),
           R"doc(
         Read a GRAIL level-0 mass file.

         Parameters
         ----------
         file_name : str
             Path to the GRAIL mass file.

         Returns
         -------
         dict
             Parsed GRAIL mass-file contents.
      )doc" );
    m.def( "grail_mass_level_1_file_reader",
           &tio::grailMassLevel1FileReader,
           py::arg( "file_name" ),
           py::arg( "data_level" ) = "1b",
           R"doc(
         Read a GRAIL level-1 mass file.

         Parameters
         ----------
         file_name : str
             Path to the GRAIL mass file.
         data_level : str, default="1b"
             Data level identifier used by the reader.

         Returns
         -------
         dict
             Parsed GRAIL mass-file contents.
      )doc" );
};

}  // namespace grail
}  // namespace missions
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
