import numpy
import pybind11_stubgen.typing_ext
from ...dynamics import parameters_setup
import typing
__all__ = ['EstimatableParameterSet', 'print_parameter_names']

class EstimatableParameterSet:
    """Class containing a consolidated set of estimatable parameters.
    
    Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
    The user typically creates instances of this class via the :func:`~tudatpy.dynamics.parameters_setup.create_parameter_set` function."""

    def indices_for_parameter_type(self, parameter_type: tuple[parameters_setup.EstimatableParameterTypes, tuple[str, str]]) -> list[tuple[int, int]]:
        """
                 Function to retrieve the indices of a given type of parameter.
        
                 Function to retrieve the index of all parameters of a given type from the parameter set.
                 This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.
        
        
                 Parameters
                 ----------
                 parameter_type : Tuple[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, Tuple[str, str] ]
                     help
                 Returns
                 -------
                 List[ Tuple[int, int] ]
                     help
        """

    @property
    def constraints_size(self) -> int:
        """
                 **read-only**
        
                 Total size of linear constraint that is to be applied during estimation.
        
                 :type: int
        """

    @property
    def initial_multi_arc_states_size(self) -> int:
        """
                 **read-only**
        
                 Amount of initial state parameters in the set, which are treated in a multi-arc fashion.
        
                 :type: int
        """

    @property
    def initial_single_arc_states_size(self) -> int:
        """
                 **read-only**
        
                 Amount of initial state parameters in the set, which are treated in a single-arc fashion.
        
                 :type: int
        """

    @property
    def initial_states_size(self) -> int:
        """
                 **read-only**
        
                 Amount of initial state parameters contained in the set.
        
                 :type: int
        """

    @property
    def parameter_set_size(self) -> int:
        """
                 **read-only**
        
                 Size of the parameter set, i.e. amount of estimatable parameters contained in the set.
        
                 :type: int
        """

    @property
    def parameter_vector(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 Vector containing the parameter values of all parameters in the set.
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]
        """

    @parameter_vector.setter
    def parameter_vector(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        ...

def print_parameter_names(parameter_set: EstimatableParameterSet) -> None:
    """Function for printing a list of estimatable parameter names.
    
    Function that allows you to print a verbose list of all parameters that shall be estimated. Consider parameters are listed separately.
    
    
    Parameters
    ----------
    parameter_set : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`.
        Instance of :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.
    
    Returns
    -------
    List[]
        Verbose List of all parameters that shall be estimated. Consider parameters are listed separately.
    
    Examples
    --------
    .. code-block:: python
    
        import numpy as np
        from tudatpy.interface import spice
        from tudatpy.dynamics import environment_setup, propagation_setup, parameters_setup, parameters
        from tudatpy
        from tudatpy.astro.time_representation import DateTime
    
        # Load SPICE kernels
        spice.load_standard_kernels()
    
        # Set simulation epochs
        simulation_start_epoch = DateTime(2000, 1, 1).epoch()
        simulation_end_epoch = DateTime(2000, 1, 4).epoch()
    
        # Create bodies
        bodies_to_create = ["Sun", "Earth"]
        global_frame_origin = "Earth"
        global_frame_orientation = "J2000"
        body_settings = environment_setup.get_default_body_settings(
          bodies_to_create, global_frame_origin, global_frame_orientation)
    
        # Create vehicle
        body_settings.add_empty_settings("Delfi-C3")
        body_settings.get("Delfi-C3").constant_mass = 2.2
        bodies = environment_setup.create_system_of_bodies(body_settings)
    
        # Define propagation settings
        bodies_to_propagate = ["Delfi-C3"]
        central_bodies = ["Earth"]
        accelerations_settings_delfi_c3 = dict(
          Sun=[propagation_setup.acceleration.point_mass_gravity()],
          Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(8, 8)]
        )
        acceleration_settings = {"Delfi-C3": accelerations_settings_delfi_c3}
        acceleration_models = propagation_setup.create_acceleration_models(
          bodies, acceleration_settings, bodies_to_propagate, central_bodies)
        initial_state = np.zeros(6)  # Use a real initial state if needed
    
        # Integrator settings
        integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(60.0,
                                                                                       coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87)
    
        # Create propagator
        termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)
        propagator_settings = propagation_setup.propagator.translational(
          central_bodies, acceleration_models, bodies_to_propagate, initial_state,
          simulation_start_epoch, integrator_settings, termination_condition)
    
        # Define parameters to estimate
        parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
        parameter_settings.append(parameters_setup.gravitational_parameter("Earth"))
        parameters_to_estimate = parameters_setup.create_parameter_set(parameter_settings, bodies)
    
        # Print parameter names
        print(parameters.print_parameter_names(parameters_to_estimate))"""