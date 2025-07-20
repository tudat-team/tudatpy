import typing
__all__ = ['MaximumIterationHandling', 'NewtonRaphsonCore', 'RootFinderCore', 'RootFinderSettings', 'accept_result', 'accept_result_with_warning', 'bisection', 'halley', 'newton_raphson', 'secant', 'throw_exception']

class MaximumIterationHandling:
    """Enumeration of types of behaviour to be used when the convergence criterion on maximum number of iterations is reached.
    
    
    
    
    
          
    
    Members:
    
      accept_result : 
    
    The program will accept the root at the final iteration, without any additional output.
    
    
    
      accept_result_with_warning : 
    
    The program will accept the root at the final iteration, but will print a warning to the terminal that the root finder may not have converged.
    
          
    
      throw_exception : 
    
    The program will not accept the root at the final iteration, and will throw a :class:`~tudatpy.exceptions.MaximumIterationsExceededError` error."""
    __members__: typing.ClassVar[dict[str, MaximumIterationHandling]]
    accept_result: typing.ClassVar[MaximumIterationHandling]
    accept_result_with_warning: typing.ClassVar[MaximumIterationHandling]
    throw_exception: typing.ClassVar[MaximumIterationHandling]

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

    def __init__(self, x_tol: float, max_iter: int) -> None:
        ...

class RootFinderCore:
    pass

class RootFinderSettings:
    """Class to define settings for a root finder."""

def bisection(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    """Function to create settings for a bisection root-finder.
    
    Function to create settings for a bisection root finder. This root finder approximates the root by initializing with
    two initial guesses :math:`x_{\\downarrow,0}` and :math:`x_{\\uparrow,0}`, for which it is required that
    :math:`f(x_{\\downarrow,0}) < 0` and :math:`f(x_{\\uparrow,0}) > 0`. At each iteration :math:`i`, the current guess of
    the root :math:`x_{i}` is:
    
    .. math::
       x_{i}=\\begin{cases}
       x_{\\downarrow,i}, & |f(x_{\\downarrow,i})|<|f(x_{\\uparrow,i})|\\\\
          x_{\\uparrow,i}, & |f(x_{\\downarrow,i})|\\ge|f(x_{\\uparrow,i})|
               \\end{cases}
    
    The midpoint :math:`x_{m,i}` of :math:`x_{\\downarrow,i}` and :math:`x_{\\uparrow,i}` is then computed from :math:`x_{m,i}=(x_{\\downarrow,i}-x_{\\uparrow,i})/2`.
    Depending on the sign of :math:`f(x_{m,i})`, it then replaces either :math:`x_{\\downarrow,i}` or :math:`x_{\\uparrow,i}` (depending on whether
    its sign matches :math:`f(x_{\\downarrow,i})` for iteration :math:`i+1` and :math:`f(x_{\\uparrow,i})`), while the other point from iteration :math:`i` is retained.
    
    Although slow, the algorithm is ensured to converge to a root, if the two initial guesses indeed have opposite signs (if not, an exception is thrown).
    
    
    Parameters
    ----------
    relative_variable_tolerance : float, default = nan
        Relative tolerance :math:`\\epsilon_{r}` (setting not used if nan)
    absolute_variable_tolerance : float, default = nan
        Relative absolute :math:`\\epsilon_{a}` (setting not used if nan)
    root_function_tolerance : float, default = nan
        Root function tolerance :math:`\\epsilon_{f}` (setting not used if nan)
    maximum_iteration : int, default = 1000
        Maximum number of iterations :math:`N`
    maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
        Algorithm behaviour if maximum number of iterations :math:`N` is reached
    Returns
    -------
    RootFinderSettings
        Bisection root-finding settings object"""

def halley(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    """Function to create settings for a Halley root-finder.
    
    Function to create settings for a Halley root finder. This root finder approximates the root by initializing with
    a single initial guesses :math:`x_{0}` and requires an analytical formulation for :math:`f(x)`, :math:`f'(x)=\\frac{d}{dx}f(x)` and :math:`f''(x)=\\frac{d^{2}}{dx^{2}}f(x)`.
    The algorithm uses the following equation to iterate:
    
    .. math::
       x_{i+1}=x_{i}-\\frac{2f(x_{i})f'(x_{i})}{2(f'(x_{i}))^{2}-f(x_{i})f''(x_{i})}
    
    
    
    
    Parameters
    ----------
    relative_variable_tolerance : float, default = nan
        Relative tolerance :math:`\\epsilon_{r}` (setting not used if nan)
    absolute_variable_tolerance : float, default = nan
        Relative absolute :math:`\\epsilon_{a}` (setting not used if nan)
    root_function_tolerance : float, default = nan
        Root function tolerance :math:`\\epsilon_{f}` (setting not used if nan)
    maximum_iteration : int, default = 1000
        Maximum number of iterations :math:`N`
    maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
        Algorithm behaviour if maximum number of iterations :math:`N` is reached
    Returns
    -------
    RootFinderSettings
        Halley root-finding settings object"""

def newton_raphson(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    """Function to create settings for a Newton-Raphson root-finder.
    
    Function to create settings for a bisection root finder. This root finder approximates the root by initializing with
    a single initial guesses :math:`x_{0}` and requires an analytical formulation for :math:`f(x)` and :math:`f'(x)=\\frac{d}{dx}f(x)`.
    The algorithm uses the following equation to iterate:
    
    .. math::
       x_{i+1}=x_{i}-\\frac{f(x_{i})}{f'(x_{i})}
    
    
    Parameters
    ----------
    relative_variable_tolerance : float, default = nan
        Relative tolerance :math:`\\epsilon_{r}` (setting not used if nan)
    absolute_variable_tolerance : float, default = nan
        Relative absolute :math:`\\epsilon_{a}` (setting not used if nan)
    root_function_tolerance : float, default = nan
        Root function tolerance :math:`\\epsilon_{f}` (setting not used if nan)
    maximum_iteration : int, default = 1000
        Maximum number of iterations :math:`N`
    maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
        Algorithm behaviour if maximum number of iterations :math:`N` is reached
    Returns
    -------
    RootFinderSettings
        Newton-Raphson root-finding settings object"""

def secant(relative_variable_tolerance: float=..., absolute_variable_tolerance: float=..., root_function_tolerance: float=..., maximum_iteration: int=1000, maximum_iteration_handling: MaximumIterationHandling=...) -> RootFinderSettings:
    """Function to create settings for a secant method root-finder.
    
    Function to create settings for a root finder using the secant method. This root finder approximates the root by initializing with
    two initial guesses :math:`x_{0}` and :math:`x_{1}`. The algorithm uses the following equation to iterate:
    
    .. math::
       x_{i+1}=x_{i}-f(x_{i})\\frac{x_{i}-x_{i-1}}{f(x_{i})-f(x_{i-1})}
    
    
    
    Parameters
    ----------
    relative_variable_tolerance : float, default = nan
        Relative tolerance :math:`\\epsilon_{r}` (setting not used if nan)
    absolute_variable_tolerance : float, default = nan
        Relative absolute :math:`\\epsilon_{a}` (setting not used if nan)
    root_function_tolerance : float, default = nan
        Root function tolerance :math:`\\epsilon_{f}` (setting not used if nan)
    maximum_iteration : int, default = 1000
        Maximum number of iterations :math:`N`
    maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
        Algorithm behaviour if maximum number of iterations :math:`N` is reached
    Returns
    -------
    RootFinderSettings
        Secant root-finding settings object"""
accept_result: MaximumIterationHandling
accept_result_with_warning: MaximumIterationHandling
throw_exception: MaximumIterationHandling