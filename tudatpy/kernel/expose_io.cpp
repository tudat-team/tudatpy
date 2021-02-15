/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_io.h"

#include "tudat/io/readHistoryFromFile.h"

namespace py = pybind11;

using namespace tudat::input_output;

namespace tudatpy {

void expose_io(py::module &m) {
  m.def("get_resource_path", &tudat::paths::get_resource_path);
  m.def("get_ephemeris_path", &tudat::paths::getEphemerisDataFilesPath);
  m.def("get_earth_orientation_path", &tudat::paths::getEarthOrientationDataFilesPath);
  m.def("get_quadrature_path", &tudat::paths::getQuadratureDataPath);
  m.def("get_spice_kernel_path", &tudat::paths::getSpiceKernelPath);
  m.def("get_atmosphere_tables_path", &tudat::paths::getAtmosphereTablesPath);
  m.def("get_gravity_models_path", &tudat::paths::getGravityModelsPath);
  m.def("get_space_weather_path", &tudat::paths::getSpaceWeatherDataPath);

  m.def("read_vector_history_from_file", &tudat::input_output::readVectorHistoryFromFile< double, double >,
        py::arg("vector_size"),
        py::arg("file_name") );

  m.def("read_matrix_history_from_file", &tudat::input_output::readMatrixHistoryFromFile< double, double >,
        py::arg("matrix_rows"),
        py::arg("matrix_columns"),
        py::arg("file_name") );
};

}// namespace tudatpy
