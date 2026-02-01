import functools
import warnings
from types import ModuleType
from typing import Callable


def deprecation_warning(old_name: str, new_name: str) -> None:
    warnings.warn(
        f"{old_name} is deprecated. Use {new_name} instead.",
        stacklevel=3,
    )

def object_deprecation(old_name: str, new_name: str) -> Callable:
    """
    Decorator to throw a deprecation warning when calling the object (function or class).

    Parameters
    ----------
    old_name : str
        The name of the deprecated object.
    new_name : str
        The name of the new object to use instead of the deprecated one.

    Returns
    -------
    Callable
        A decorator that wraps the original object and issues a deprecation warning when called.

    Example
    -------

    ```python
    from tudatpy._deprecation import object_deprecation

    @object_deprecation("old_function", "new_function")
    def old_function():
        pass
        
    old_function()  # This will issue a deprecation warning.
    ```
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            deprecation_warning(old_name, new_name)
            return func(*args, **kwargs)
        return wrapper
    return decorator


def register_deprecated_object(destination_module: ModuleType, old_name: str, superseding_object) -> None:
    """
    Register an object (function or class) under a deprecated name.

    This function makes an object (function or class) available under a different, deprecated name with the same functionality. This is useful e.g. when fixing a typo, where only the corrected function name remains implemented, while they deprecated function should stay available for backwards compatibility purposes.

    Parameters
    ----------
    destination_module : ModuleType
        Module where the deprecated function should be registered.
    old_name : str
        Name under which the deprecated function should be registered.
    superseding_object
        Object that super-seeds the deprecated object. This can be a function or a class.
    
    Example
    -------
    This example registers a new function `new_function` under the deprecated name `old_function` in the module where the `register_deprecated_object` function is called.
    ```python
    import sys
    from tudatpy._deprecation import register_deprecated_object
    
    def new_function():
        pass
        
    register_deprecated_object(sys.modules[__name__], "old_function", new_function)
    ```
    """

    interface = object_deprecation(old_name, superseding_object.__name__)(superseding_object)
    setattr(destination_module, old_name, interface)

def keep_typo_backwards_compatible(destination_module: ModuleType, typo: str, correct_spelling: str) -> None:
    """
    .. warning::

        This function makes rather strong assumptions and should be used with care.
    
    This is a convenience function to register functions _with_ typos to maintain backwards compatibility.

    This function assumes that:
    - the developer has fixed all typos in the function/class signatures of a module
    - backwards compatibility should be maintained, i.e. the functions should be callable with typos
    - all functions/classes that have the correct spelling had the wrong spelling previously.
    
    The function will search for all functions/classes in the given module that contain the correct spelling and register them under the name with the typo.

    Parameters
    ----------
    destination_module : ModuleType
        Module where the deprecated functions should be registered.
    typo : str
        Typo in the function/class name that should be registered as deprecated.
    correct_spelling : str
        Correct spelling of the function/class name.

    Example
    -------
    This example registers all functions/classes in the `tudatpy.kernel.estimation.observations_setup.ancillary_settings` module that contain the correct spelling "ancillary" under the deprecated name with the typo "ancilliary".

    ```python
    from tudatpy.kernel.estimation.observations_setup.ancillary_settings import *
    from tudatpy._deprecation import keep_typo_backwards_compatible

    import sys
    keep_typo_backwards_compatible(sys.modules[__name__], "ancilliary", "ancillary")
    ```

    """


    for argument in dir(destination_module):
        if correct_spelling in argument:
            old_name = argument.replace(correct_spelling, typo)
            superseding_func = getattr(destination_module, argument)
            register_deprecated_object(destination_module, old_name, superseding_func)



def property_deprecation(old_name: str, new_name: str) -> Callable:
    """
    Wrap a property to issue a deprecation warning when accessed.

    Parameters
    ----------
    old_name : str
        Deprecated name of the property.
    new_name : str
        New name of the property.

    Returns
    -------
    Callable
        A decorator that wraps the property and issues a deprecation warning when accessed.

    Example
    -------
    In this example, a `ancilliary_settings` property is added to the `SingleObservationSet` class to maintain backwards compatibility, which wraps the `ancillary_settings` property and issues a deprecation warning when accessed.

    ```python
    from tudatpy._deprecation import property_deprecation_decorator
    from tudatpy.kernel.estimation.observations import SingleObservationSet

    SingleObservationSet.ancilliary_settings = property_deprecation_decorator("SingleObservationSet.ancilliary_settings", "SingleObservationSet.ancillary_settings")(SingleObservationSet.ancillary_settings)
    ```
    """
    def decorator(prop):
        def property_deprecation(old_name: str, new_name: str, prop: property)-> property:

            def wrap_get(fget):
                if fget is None:
                    return None
                @functools.wraps(fget)
                def getter(*args, **kwargs):
                    deprecation_warning(old_name, new_name)
                    return fget(*args, **kwargs)
                return getter

            def wrap_set(fset):
                if fset is None:
                    return None
                @functools.wraps(fset)
                def setter(*args, **kwargs):
                    deprecation_warning(old_name, new_name)
                    return fset(*args, **kwargs)
                return setter

            def wrap_del(fdel):
                if fdel is None:
                    return None
                @functools.wraps(fdel)
                def deleter(*args, **kwargs):
                    deprecation_warning(old_name, new_name)
                    return fdel(*args, **kwargs)
                return deleter

            return property(
                wrap_get(prop.fget),
                wrap_set(prop.fset),
                wrap_del(prop.fdel),
                prop.__doc__,
            )
        return property_deprecation(old_name, new_name, prop)
    return decorator