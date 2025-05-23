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
#include <pybind11/eval.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "expose_exceptions.h"
#include <tudat/basics/tudatExceptions.h>

namespace py = pybind11;
namespace te = tudat::exceptions;

namespace tudatpy
{
namespace exceptions
{

void expose_exceptions( py::module& m )
{
    /*
    This module is based on the solution to exposing C++ exceptions to Python described in the following GitHub issue:
    https://github.com/pybind/pybind11/issues/1281

    At the time of implementation, it seems this is the only way to expose the
    attributes of the C++ exceptions to Python.
    See https://github.com/tudat-team/tudat/pull/297 for more details.
     */

    py::exec( R"pybind(
class TudatError(RuntimeError):
    """
    Base Error thrown by Tudat.

    """
    def __init__(self, message: str):
        super().__init__(message)

class InterpolationOutOfBoundsError(TudatError):
    """
    Error thrown when the independent variable data point is out of the bounds of the data to be interpolated.

    Attributes
    ----------
    message : str
            The error message.
    requested_value : float
            The requested value that is out of bounds.
    lower_bound : float
            The lower bound of the data to be interpolated.
    upper_bound : float
            The upper bound of the data to be interpolated.
    """

    def __init__(self, message: str, requested_value: float, lower_bound: float, upper_bound: float):
        super().__init__(message)
        self.requested_value = requested_value
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

class LagrangeInterpolationOutOfBoundsError(InterpolationOutOfBoundsError):
    """
    Error thrown when the independent variable data point of a Lagrange interpolation is outside the reliable bounds of the data to be interpolated.
    For more information, see :func:`~tudatpy.math.interpolators.lagrange_interpolation`.

    Attributes
    ----------
    message : str
            The error message.
    requested_value : float
            The requested value that is outside the reliable bounds.
    lower_bound : float
            The lower bound of the data that can be reliably interpolated.
    upper_bound : float
            The upper bound of the data that can be reliably interpolated.
    """

    def __init__(self, message: str, requested_value: float, lower_bound: float, upper_bound: float):
        super().__init__(message, requested_value, lower_bound, upper_bound)

class MaximumIterationsExceededError(TudatError):
    """
    Error thrown when the maximum number of iterations of an iterative operation is exceeded.

    Attributes
    ----------
    message : str
            The error message.
    number_of_iterations : int
            The number of iterations that have been completed.
    maximum_number_of_iterations : int
            The maximum number of iterations that was specified.

    """

    def __init__(self, message: str, number_of_iterations: int, maximum_number_of_iterations: int):
        super().__init__(message)
        self.number_of_iterations = number_of_iterations
        self.maximum_number_of_iterations = maximum_number_of_iterations

class StepSizeViolationError(TudatError):
    """
    Error thrown when the step size in a numerical integration is not valid.
    """
    def __init__(self, message: str):
        super().__init__(message)

class MinimumStepSizeViolatedError(StepSizeViolationError):
    """
    Error thrown when the step size requested by the step size controller is smaller than the defined minimum step size in a numerical integration.
    
    Attributes
    ----------
    message : str
            The error message.
    minimum_step_size : float
            The minimum step size that is allowed.
    recommended_step_size : float
            The step size recommended by the stepsize controller, which is smaller than the minimum step size.
    
    """

    def __init__(self, message: str, minimum_step_size: float, recommended_step_size: float):
        super().__init__(message)
        self.minimum_step_size = minimum_step_size
        self.recommended_step_size = recommended_step_size

        )pybind",
              m.attr( "__dict__" ),
              m.attr( "__dict__" ) );

    // Retrieve a handle the the exception type just created by executing the string literal above.
    const py::object tudatErrorExp = m.attr( "TudatError" );
    const py::object interpolationOutOfBoundsErrorExp = m.attr( "InterpolationOutOfBoundsError" );
    const py::object lagrangeInterpolationOutOfBoundsErrorExp = m.attr( "LagrangeInterpolationOutOfBoundsError" );
    const py::object maximumIterationsExceededErrorExp = m.attr( "MaximumIterationsExceededError" );
    const py::object stepSizeViolatedErrorExp = m.attr( "StepSizeViolationError" );
    const py::object minimumStepSizeViolatedErrorExp = m.attr( "MinimumStepSizeViolatedError" );

    // Register a function that will translate our C++ exceptions to the
    // appropriate Python exception type.
    //
    // Note that capturing lambdas are not allowed here, so we must
    // `import` the exception type in the body of the function.
    py::register_exception_translator( []( std::exception_ptr p ) {
        // Accept variable number of arguments after exc to populate the exception instance
        const auto setPyException = []( const char* pyTypeName, const auto& exc, auto... args ) {
            const py::object pyClass = py::module_::import( "tudatpy.exceptions" ).attr( pyTypeName );

            // Pass the exception message and any additional arguments to the exception constructor
            const py::object pyInstance = pyClass( exc.what( ), args... );
            PyErr_SetObject( pyClass.ptr( ), pyInstance.ptr( ) );
        };

        // Handle the different possible C++ exceptions, creating the
        // corresponding Python exception and setting it as the active
        // exception in this thread.
        try
        {
            if( p ) std::rethrow_exception( p );
        }
        catch( const te::MaximumIterationsExceededError& exc )
        {
            setPyException( "MaximumIterationsExceededError", exc, exc.numberOfIterations, exc.maximumNumberOfIterations );
        }
        catch( const te::MinimumStepSizeViolatedError< TIME_TYPE >& exc )
        {
            setPyException( "MinimumStepSizeViolatedError", exc, exc.minimumStepSize, exc.recommendedStepSize );
        }
        catch( const te::StepSizeViolationError& exc )
        {
            setPyException( "StepSizeViolationError", exc );
        }
        catch( const te::LagrangeInterpolationOutOfBoundsError< TIME_TYPE >& exc )
        {
            setPyException( "LagrangeInterpolationOutOfBoundsError", exc, exc.requestedValue, exc.lowerBound, exc.upperBound );
        }
        catch( const te::InterpolationOutOfBoundsError< TIME_TYPE >& exc )
        {
            setPyException( "InterpolationOutOfBoundsError", exc, exc.requestedValue, exc.lowerBound, exc.upperBound );
        }
        catch( const te::TudatError& exc )
        {
            setPyException( "TudatError", exc );
        }
    } );
}

}  // namespace exceptions
}  // namespace tudatpy
