/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "expose_util.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "tudat/io/readHistoryFromFile.h"

namespace py = pybind11;

namespace tudatpy
{
namespace util
{

void expose_util( py::module& m )
{
    m.def( "read_vector_history_from_file",
           &tudat::input_output::readVectorHistoryFromFile< double, double >,
           py::arg( "vector_size" ),
           py::arg( "file_name" ),
           R"doc(
         Read a vector history from a file.

         Parameters
         ----------
         vector_size : int
             Size of the vector at each epoch.
         file_name : str
             Name of the file containing the vector history.

         Returns
         -------
         dict[float, numpy.ndarray]
             Dictionary mapping epochs to the vector at the given epoch.
      )doc" );

    m.def( "read_matrix_history_from_file",
           &tudat::input_output::readMatrixHistoryFromFile< double, double >,
           py::arg( "matrix_rows" ),
           py::arg( "matrix_columns" ),
           py::arg( "file_name" ),
           R"doc(

 Read a matrix history from a file.


 Parameters
 ----------
 matrix_rows : int
     Number of rows in the matrix at each epoch.
 matrix_columns : int
     Number of columns in the matrix at each epoch.
 file_name : str
     Name of the file containing the matrix history.
 Returns
 -------
 Dict[float, numpy.ndarray]
     Dictionary mapping epochs to the matrix at the given epoch.


     )doc" );
}

}  // namespace util
}  // namespace tudatpy
