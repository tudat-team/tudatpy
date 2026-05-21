/*    Copyright (c) 2010-2021, Delft University of Technology
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
#include "expose_observations_wrapper.h"

#include "expose_observations_wrapper_bindings.h"

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_wrapper
{

void expose_observations_wrapper( py::module& m )
{
    py::module_::import( "tudatpy.estimation.observations" ).attr( "ObservationCollection" );

    expose_observations_wrapper_io_bindings( m );
    expose_observations_wrapper_simulation_bindings( m );
}

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
