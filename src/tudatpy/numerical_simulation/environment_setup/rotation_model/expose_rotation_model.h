/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_ROTATION_MODEL_SETUP_H
#define TUDATPY_EXPOSE_ROTATION_MODEL_SETUP_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{
namespace rotation_model
{

void expose_rotation_model_setup( py::module &m );

}  // namespace rotation_model
}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_ROTATION_MODEL_SETUP_H
