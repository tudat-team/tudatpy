/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include <pybind11/pybind11.h>

#include "scalarTypes.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace ephemeris
{
namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
{

void expose_ephemeris_setup_state_scalar( py::module& m )
{
    m.def( "create_ephemeris",
           &tss::createBodyEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ephemeris_settings" ),
           py::arg( "body_name" ) );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace ephemeris
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
