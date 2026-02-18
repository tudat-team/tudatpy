.. _random_noise:

``random_noise``
=====================

When simulating observations, it is often useful to add random noise to the observations to emulate the behaviour of real data. The functions in this module are used to settings for this noise to the simulation of observations. At present, simple Python functions are used directly as input. The noise functions take the time as input (typically as a ``float``) and return the associated randomly generated noise as a ``np.ndarray`` (rather than a ``float``, to permit consistency with multi-valued observations; for single-valued observations such as range or Doppler the noise function should return an array of size 1).

The main interface with Tudat is that these noise functions are can be added to :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects, which are in turn used as input to observation simulation (see :ref:`observations_simulation_settings`). The noise function can be set in these objects directly using the :attr:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings.noise_function` attribute. For convenience there are options to set a custom-defined noise function to a list of simulation settings (using the :func:`~tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_all`, :func:`~tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_observable` and :func:`~tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_observable_for_link_ends` functions). Similarly, there are functions to create a Gaussian noise function (with zero mean) and add this to the simulation settings (sing the :func:`~tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_all`, :func:`~tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_observable` and :func:`~tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_observable_for_link_ends` functions).

As an example, adding Gaussian noise to all observations of type ``one_way_range_type`` is done by:

.. code-block:: python

    from tudatpy.estimation.observations_setup import random_noise
    from tudatpy.estimation.observable_models_setup import model_settings
    
    # Create a list[ObservationSimulationSettings] 
    observation_simulation_settings_list = ....  
    
    # Add random noise settings
    noise_level = 0.1
    random_noise.add_gaussian_noise_to_observable(
        observation_simulation_settings_list,
        noise_level,
        model_settings.one_way_range_type
    )

which will add 0.1 m standard deviation Gaussian random noise to each one-way range setting in the ``observation_simulation_settings_list`` list. In this case (the :func:`~tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_observable` function).

An example of doing the exact same procedure, but using the custom noise function interface:

.. code-block:: python

    import numpy as np
    from tudatpy.estimation.observations_setup import random_noise
    from tudatpy.estimation.observable_models_setup import model_settings
    
    def custom_noise_function(current_time):
        return np.array([np.random.lognormal(0.0, 1.0)])
    
    random_noise.add_noise_function_to_observable(
        observation_simulation_settings_list,
        custom_noise_function,
        model_settings.one_way_range_type
    )

where it is important to realize that the noise function *must* have a single float representing time as input, and returns a vector (of the size of a single observation) as output (for range, this size is equal to 1)

Functions
---------
.. currentmodule:: tudatpy.estimation.observations_setup.random_noise

.. autosummary::

   add_noise_function_to_all

   add_noise_function_to_observable

   add_noise_function_to_observable_for_link_ends

   add_gaussian_noise_to_all

   add_gaussian_noise_to_observable

   add_gaussian_noise_to_observable_for_link_ends


.. autofunction:: tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_all

.. autofunction:: tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_observable

.. autofunction:: tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_observable_for_link_ends

.. autofunction:: tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_all

.. autofunction:: tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_observable

.. autofunction:: tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_observable_for_link_ends