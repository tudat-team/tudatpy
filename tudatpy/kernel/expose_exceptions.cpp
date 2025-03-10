/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "expose_exceptions.h"
#include <tudat/basics/tudatExceptions.h>

namespace py = pybind11;
namespace te = tudat::exceptions;

namespace tudatpy
{
namespace exceptions
{

void expose_exceptions( py::module &m )
{
    auto tudatErrorExp = py::register_exception< te::TudatError >( m, "TudatError", PyExc_RuntimeError );
    tudatErrorExp.doc( ) = R"(Base Error thrown by Tudat.)";

    auto interpolationOutOfBoundsErrorExp =
            py::register_exception< te::InterpolationOutOfBoundsError >( m, "InterpolationOutOfBoundsError", tudatErrorExp.ptr( ) );
    interpolationOutOfBoundsErrorExp.doc( ) = R"(Error thrown when the interpolation is out of bounds.)";

    auto lagrangeInterpolationOutOfBoundsErrorExp = py::register_exception< te::LagrangeInterpolationOutOfBoundsError >(
            m, "LagrangeInterpolationOutOfBoundsError", interpolationOutOfBoundsErrorExp.ptr( ) );
    lagrangeInterpolationOutOfBoundsErrorExp.doc( ) = R"(Error thrown when the Lagrange interpolation is out of bounds.)";

    auto maximumIterationsExceededErrorExp =
            py::register_exception< te::MaximumIterationsExceededError >( m, "MaximumIterationsExceededError", tudatErrorExp.ptr( ) );
    maximumIterationsExceededErrorExp.doc( ) = R"(Error thrown when the maximum number of iterations is exceeded.)";

    auto stepSizeViolatedErrorExp =
            py::register_exception< te::StepSizeViolationError >( m, "StepSizeViolationError", tudatErrorExp.ptr( ) );
    stepSizeViolatedErrorExp.doc( ) = R"(Error thrown when the step size in a numerical integration is not valid.)";

    auto minimumStepSizeViolatedErrorExp = py::register_exception< te::MinimumStepSizeViolatedError >(
            m, "MinimumStepSizeViolatedError", stepSizeViolatedErrorExp.ptr( ) );
    minimumStepSizeViolatedErrorExp.doc( ) = R"(Error thrown when the defined minimum step size is violated in a numerical integration.)";
}

}  // namespace exceptions
}  // namespace tudatpy
