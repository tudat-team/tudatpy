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
#include "expose_propagation_bindings.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <tudat/astro/propulsion/thrustMagnitudeWrapper.h>

namespace py = pybind11;
namespace tpr = tudat::propulsion;

namespace tudatpy
{
namespace dynamics
{
namespace propagation
{

void expose_propagation_thrust_types( py::module& m )
{
    py::class_< tpr::ThrustMagnitudeWrapper, std::shared_ptr< tpr::ThrustMagnitudeWrapper > >( m, "ThrustMagnitudeWrapper" , R"doc(

         Base class for computing the thrust magnitude and mass flow rate of an engine.

         Object supplying the thrust magnitude and specific impulse used to evaluate the force and propellant consumption.
         Instances of this class are associated with an :class:`~tudatpy.dynamics.environment.EngineModel`.

      )doc" );
}

void expose_propagation_thrust_bindings( py::module& m )
{
    py::class_< tpr::ConstantThrustMagnitudeWrapper, std::shared_ptr< tpr::ConstantThrustMagnitudeWrapper >, tpr::ThrustMagnitudeWrapper >(
            m, "ConstantThrustMagnitudeWrapper" , R"doc(

         Object defining a constant engine thrust magnitude and specific impulse.

         :class:`~ThrustMagnitudeWrapper` derived class for an engine whose thrust magnitude can be reset through the
         :attr:`~ConstantThrustMagnitudeWrapper.constant_thrust_magnitude` property.

      )doc" )
            .def_property( "constant_thrust_magnitude",
                           &tpr::ConstantThrustMagnitudeWrapper::getConstantThrustForceMagnitude,
                           &tpr::ConstantThrustMagnitudeWrapper::resetConstantThrustForceMagnitude , R"doc(

         Constant thrust force magnitude produced by the engine, in N.

         :type: float

      )doc" );

    py::class_< tpr::CustomThrustMagnitudeWrapper, std::shared_ptr< tpr::CustomThrustMagnitudeWrapper >, tpr::ThrustMagnitudeWrapper >(
            m, "CustomThrustMagnitudeWrapper" , R"doc(

         Object defining an engine thrust magnitude from a user-provided function.

         :class:`~ThrustMagnitudeWrapper` derived class evaluating the thrust magnitude and specific impulse as functions of time.
         The thrust magnitude function can be replaced using :attr:`~CustomThrustMagnitudeWrapper.custom_thrust_magnitude`.

      )doc" )
            .def_property( "custom_thrust_magnitude", nullptr, &tpr::CustomThrustMagnitudeWrapper::resetThrustMagnitudeFunction , R"doc(

         **write-only**

         Function returning the thrust magnitude, in N, at a time in seconds since J2000. Assigning this property replaces the current thrust magnitude function.

         :type: Callable[[float], float]

      )doc" );
}

}  // namespace propagation
}  // namespace dynamics
}  // namespace tudatpy
