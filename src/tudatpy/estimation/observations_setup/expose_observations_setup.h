/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_OBSERVATIONS_SETUP_H
#define TUDATPY_EXPOSE_OBSERVATIONS_SETUP_H

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ancillary_settings/expose_ancillary_settings.h"
#include "observations_dependent_variables/expose_observations_dependent_variables.h"
#include "observations_simulation_settings/expose_observations_simulation_settings.h"
#include "observations_wrapper/expose_observations_wrapper.h"
#include "random_noise/expose_random_noise.h"
#include "viability/expose_viability.h"


namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{

void expose_observations_setup( py::module &m );

}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_OBSERVATIONS_SETUP_H
