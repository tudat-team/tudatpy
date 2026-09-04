/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_environment_state_scalar.h"

#include "scalarTypes.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace environment
{
namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
{

void expose_environment_state_scalar( BodyPythonBinding& bodyBinding, SystemOfBodiesPythonBinding& systemOfBodiesBinding )
{
    bodyBinding.def( "state_in_base_frame_from_ephemeris",
                     &tss::Body::getStateInBaseFrameFromEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
                     py::arg( "time" ),
                     R"doc(
Return the body's ephemeris state, translated to the system's global frame origin when needed.

Parameters
----------
time : astro.time_representation.Time
    Epoch at which the Cartesian state is evaluated.
      )doc" );

    systemOfBodiesBinding
            .def( "create_empty_body",
                  &tss::SystemOfBodies::createEmptyBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = true,
                  R"doc(
Create an empty body and add it to this system of bodies.
      )doc" )
            .def( "add_body",
                  &tss::SystemOfBodies::addBody< STATE_SCALAR_TYPE, TIME_TYPE >,
                  py::arg( "body_to_add" ),
                  py::arg( "body_name" ),
                  py::arg( "process_body" ) = true,
                  R"doc(
Add an existing body to this system of bodies.
      )doc" );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace environment
}  // namespace dynamics
}  // namespace tudatpy
