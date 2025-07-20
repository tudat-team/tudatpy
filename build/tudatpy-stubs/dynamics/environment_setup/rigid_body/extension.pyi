import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['RigidBodyPropertiesSettings', 'constant_rigid_body_properties', 'custom_mass_dependent_rigid_body_properties', 'custom_time_dependent_rigid_body_properties']

class RigidBodyPropertiesSettings:
    """Base class for providing settings for rigid body model creation.
    
    This class is a functional base class for settings of gravity field models that require no information in addition to their type.
    Gravity field model classes requiring additional information must be created using an object derived from this class."""

    @property
    def body_mass_property_type(self) -> ...:
        """
                 **read-only**
        
                 Type of rigid body model that is to be created.
        
                 :type: RigidBodyPropertiesType
        """

def constant_rigid_body_properties(mass: float, center_of_mass: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]=..., inertia_tensor: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]=...) -> RigidBodyPropertiesSettings:
    """Function for creating constant rigid body properties.
    
    Function for creating constant rigid body properties (mass, center of mass, inertia tensor). The center of mass and/or inertia tensor can be left empty by setting them
    to NaN (default), in which case no center of mass or inertia tensor are defined
    
    
    Parameters
    ----------
    mass : float
        Constant mass of the body
    center_of_mass : np.array, default = np.full([3, 1], np.nan)
        Constant center of mass of the body (in a body-fixed frame)
    inertia_tensor : np.array, default = np.full([3, 3], np.nan)
        Constant inertia tensor of the body (in a body-fixed frame)
    Returns
    -------
    RigidBodyPropertiesSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings"""

def custom_mass_dependent_rigid_body_properties(mass: float, center_of_mass_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=None, inertia_tensor_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]]=None) -> RigidBodyPropertiesSettings:
    """Function for creating custom (time-dependent) rigid body properties.
    
    Function for creating custom rigid body properties, center of mass and inertia tensor are defined by user-defined functions as a function of mass.
    This functionality is typically used for a body under thrust, where the center of mass and inertia tensor are defined as a function of expended mass.
    
    
    Parameters
    ----------
    mass : Callable[[float], float]
        Mass of the body (to be overridden during propagation if mass is propagated)
    center_of_mass_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]] = None
        Function returning the center of mass as a function of mass (to be used during the propagation)
    inertia_tensor_function : Callable[[float], numpy.ndarray[numpy.float64[3, 3]]] = None
        Function returning the inertia tensor as a function of mass (to be used during the propagation)
    Returns
    -------
    RigidBodyPropertiesSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings"""

def custom_time_dependent_rigid_body_properties(mass_function: typing.Callable[[float], float], center_of_mass_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=None, inertia_tensor_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]]=None) -> RigidBodyPropertiesSettings:
    """Function for creating custom (time-dependent) rigid body properties.
    
    Function for creating custom rigid body properties, where the mass, center of mass and inertia tensor are defined by user-defined functions (as a function of time).
    The center of mass and/or inertia tensor functions can be left empty by setting them
    to None (default), in which case no center of mass or inertia tensor are defined
    
    
    Parameters
    ----------
    mass_function : Callable[[float], float]
        Function returning the mass as a function of time (to be used during the propagation)
    center_of_mass_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]] = None
        Function returning the center of mass as a function of time (to be used during the propagation)
    inertia_tensor_function : Callable[[float], numpy.ndarray[numpy.float64[3, 3]]] = None
        Function returning the inertia tensor as a function of time (to be used during the propagation)
    Returns
    -------
    RigidBodyPropertiesSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings"""