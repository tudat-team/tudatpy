import numpy
import typing
__all__ = ['RigidBodyPropertiesSettings', 'constant_rigid_body_properties', 'custom_mass_dependent_rigid_body_properties', 'custom_time_dependent_rigid_body_properties']

class RigidBodyPropertiesSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def body_mass_property_type(self) -> ...:
        ...

def constant_rigid_body_properties(mass: float, center_of_mass: numpy.ndarray=..., inertia_tensor: numpy.ndarray=...) -> RigidBodyPropertiesSettings:
    ...

def custom_mass_dependent_rigid_body_properties(mass: float, center_of_mass_function: typing.Callable[[float], numpy.ndarray]=None, inertia_tensor_function: typing.Callable[[float], numpy.ndarray]=None) -> RigidBodyPropertiesSettings:
    ...

def custom_time_dependent_rigid_body_properties(mass_function: typing.Callable[[float], float], center_of_mass_function: typing.Callable[[float], numpy.ndarray]=None, inertia_tensor_function: typing.Callable[[float], numpy.ndarray]=None) -> RigidBodyPropertiesSettings:
    ...