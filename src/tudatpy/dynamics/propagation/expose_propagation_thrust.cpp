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

void expose_propagation_thrust_bindings( py::module& m )
{
    py::class_< tpr::ThrustMagnitudeWrapper, std::shared_ptr< tpr::ThrustMagnitudeWrapper > >(
            m, "ThrustMagnitudeWrapper" );

    py::class_< tpr::ConstantThrustMagnitudeWrapper,
                std::shared_ptr< tpr::ConstantThrustMagnitudeWrapper >,
                tpr::ThrustMagnitudeWrapper >( m, "ConstantThrustMagnitudeWrapper" )
            .def_property(
                    "constant_thrust_magnitude",
                    &tpr::ConstantThrustMagnitudeWrapper::getConstantThrustForceMagnitude,
                    &tpr::ConstantThrustMagnitudeWrapper::resetConstantThrustForceMagnitude );

    py::class_< tpr::CustomThrustMagnitudeWrapper,
                std::shared_ptr< tpr::CustomThrustMagnitudeWrapper >,
                tpr::ThrustMagnitudeWrapper >( m, "CustomThrustMagnitudeWrapper" )
            .def_property( "custom_thrust_magnitude",
                           nullptr,
                           &tpr::CustomThrustMagnitudeWrapper::resetThrustMagnitudeFunction );
}

}  // namespace propagation
}  // namespace dynamics
}  // namespace tudatpy
