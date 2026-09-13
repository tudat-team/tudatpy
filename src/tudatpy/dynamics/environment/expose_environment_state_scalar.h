#ifndef TUDATPY_EXPOSE_ENVIRONMENT_STATE_SCALAR_H
#define TUDATPY_EXPOSE_ENVIRONMENT_STATE_SCALAR_H

#include <pybind11/pybind11.h>

#include "tudat/simulation/environment_setup/body.h"

namespace py = pybind11;

namespace tudatpy
{
namespace dynamics
{
namespace environment
{

using BodyPythonBinding = py::class_< tudat::simulation_setup::Body, std::shared_ptr< tudat::simulation_setup::Body > >;
using SystemOfBodiesPythonBinding =
        py::class_< tudat::simulation_setup::SystemOfBodies, std::shared_ptr< tudat::simulation_setup::SystemOfBodies > >;

void expose_environment_state_scalar( BodyPythonBinding& bodyBinding, SystemOfBodiesPythonBinding& systemOfBodiesBinding );

}  // namespace environment
}  // namespace dynamics
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_ENVIRONMENT_STATE_SCALAR_H
