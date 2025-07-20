import typing
from . import spice_exceptions
__all__ = ['InterpolationOutOfBoundsError', 'LagrangeInterpolationOutOfBoundsError', 'MaximumIterationsExceededError', 'MinimumStepSizeViolatedError', 'StepSizeViolationError', 'TudatError', 'spice_exceptions']

class InterpolationOutOfBoundsError(TudatError):
    """Error thrown when the independent variable data point is out of the bounds of the data to be interpolated.
    
    Attributes
    ----------
    message : str
            The error message.
    requested_value : float
            The requested value that is out of bounds.
    lower_bound : float
            The lower bound of the data to be interpolated.
    upper_bound : float
            The upper bound of the data to be interpolated."""

    def __init__(self, message: str, requested_value: float, lower_bound: float, upper_bound: float):
        ...

class LagrangeInterpolationOutOfBoundsError(InterpolationOutOfBoundsError):
    """Error thrown when the independent variable data point of a Lagrange interpolation is outside the reliable bounds of the data to be interpolated.
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
            The upper bound of the data that can be reliably interpolated."""

    def __init__(self, message: str, requested_value: float, lower_bound: float, upper_bound: float):
        ...

class MaximumIterationsExceededError(TudatError):
    """Error thrown when the maximum number of iterations of an iterative operation is exceeded.
    
    Attributes
    ----------
    message : str
            The error message.
    number_of_iterations : int
            The number of iterations that have been completed.
    maximum_number_of_iterations : int
            The maximum number of iterations that was specified."""

    def __init__(self, message: str, number_of_iterations: int, maximum_number_of_iterations: int):
        ...

class MinimumStepSizeViolatedError(StepSizeViolationError):
    """Error thrown when the step size requested by the step size controller is smaller than the defined minimum step size in a numerical integration.
    
    Attributes
    ----------
    message : str
            The error message.
    minimum_step_size : float
            The minimum step size that is allowed.
    recommended_step_size : float
            The step size recommended by the stepsize controller, which is smaller than the minimum step size."""

    def __init__(self, message: str, minimum_step_size: float, recommended_step_size: float):
        ...

class StepSizeViolationError(TudatError):
    """Error thrown when the step size in a numerical integration is not valid."""

    def __init__(self, message: str):
        ...

class TudatError(RuntimeError):
    """Base Error thrown by Tudat."""

    def __init__(self, message: str):
        ...