import functools
import warnings
from types import ModuleType


def object_deprecation_decorator(old_name: str, new_name: str):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn(
                f"{old_name} is deprecated. Use {new_name} instead.",
                stacklevel=2,
            )
            return func(*args, **kwargs)
        return wrapper
    return decorator


def deprecate_object(destination_module: ModuleType, old_name: str, superseding_object):
    interface = object_deprecation_decorator(old_name, superseding_object.__name__)(superseding_object)
    setattr(destination_module, old_name, interface)

def keep_typo_backwards_compatible(destination_module: ModuleType, typo: str, correct: str):
    for argument in dir(destination_module):
        if correct in argument:
            old_name = argument.replace(correct, typo)
            superseding_func = getattr(destination_module, argument)
            deprecate_object(destination_module, old_name, superseding_func)

def deprecated_property(old_name: str, new_name: str, prop: property):
    if not isinstance(prop, property):
        raise TypeError("deprecated_property expects a property")

    def warn():
        warnings.warn(
            f"{old_name} is deprecated. Use {new_name} instead.",
            stacklevel=2,
        )

    def wrap_get(fget):
        if fget is None:
            return None
        @functools.wraps(fget)
        def getter(*args, **kwargs):
            warn()
            return fget(*args, **kwargs)
        return getter

    def wrap_set(fset):
        if fset is None:
            return None
        @functools.wraps(fset)
        def setter(*args, **kwargs):
            warn()
            return fset(*args, **kwargs)
        return setter

    def wrap_del(fdel):
        if fdel is None:
            return None
        @functools.wraps(fdel)
        def deleter(*args, **kwargs):
            warn()
            return fdel(*args, **kwargs)
        return deleter

    return property(
        wrap_get(prop.fget),
        wrap_set(prop.fset),
        wrap_del(prop.fdel),
        prop.__doc__,
    )

def property_deprecation_decorator(old_name: str, new_name: str):
    def decorator(prop):
        return deprecated_property(old_name, new_name, prop)
    return decorator