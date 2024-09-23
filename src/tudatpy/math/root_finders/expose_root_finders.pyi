import typing
import typing
__all__ = ['MaximumIterationHandling', 'NewtonRaphsonCore', 'RootFinderCore', 'RootFinderSettings', 'accept_result', 'accept_result_with_warning', 'bisection', 'halley', 'newton_raphson', 'secant', 'throw_exception']

class MaximumIterationHandling:
    """Members:
	
	accept_result :
	
	accept_result_with_warning :
	
	throw_exception :
	"""
    __members__: typing.ClassVar[dict[str, MaximumIterationHandling]]
    accept_result: typing.ClassVar[MaximumIterationHandling]
    accept_result_with_warning: typing.ClassVar[MaximumIterationHandling]
    throw_exception: typing.ClassVar[MaximumIterationHandling]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class NewtonRaphsonCore(RootFinderCore):

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, x_tol: float, max_iter: int) -> None:
        ...

class RootFinderCore:

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RootFinderSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

def bisection(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    ...

def halley(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    ...

def newton_raphson(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    ...

def secant(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    ...
accept_result: MaximumIterationHandling
accept_result_with_warning: MaximumIterationHandling
throw_exception: MaximumIterationHandling