import typing
__all__ = ['BaseFunctionHodographicShaping', 'hodograph_constant', 'hodograph_cosine', 'hodograph_exponential', 'hodograph_exponential_cosine', 'hodograph_exponential_sine', 'hodograph_power', 'hodograph_power_cosine', 'hodograph_power_sine', 'hodograph_scaled_exponential', 'hodograph_scaled_exponential_cosine', 'hodograph_scaled_exponential_sine', 'hodograph_scaled_power', 'hodograph_scaled_power_cosine', 'hodograph_scaled_power_sine', 'hodograph_sine', 'recommended_axial_hodograph_functions', 'recommended_normal_hodograph_functions', 'recommended_radial_hodograph_functions']

class BaseFunctionHodographicShaping:
    """Base class for defining settings of the shape functions for hodographic shaping method.
    
    Base class for defining settings of the shape functions for Hodograph shaping method. Objects derived
    from this class are created by calling the dedicated functions in this module"""

def hodograph_constant() -> BaseFunctionHodographicShaping:
    """Function for creating a constant contribution to hodographic trajectory shaping.
    
    Function for creating a constant contribution to hodographic trajectory shaping. This adds a contribution
    :math:`K` to the selected velocity component, with :math:`K` a free parameter.
    
    Returns
    -------
    BaseFunctionHodographicShaping
        Settings object for a constant contribution to hodographic shaping."""

def hodograph_cosine(frequency: float) -> BaseFunctionHodographicShaping:
    """Function for creating a cosine contribution to hodographic trajectory shaping.
    
    Function for creating a cosine contribution to hodographic trajectory shaping. For a
    provided frequency :math:`f`, this adds a contribution :math:`K\\cos(f\\cdot T)` to the selected
    velocity component, with :math:`T` the time since departure, and :math:`K` a free parameter.
    
    
    Parameters
    ----------
    frequency : float
        Frequency of the cosine contribution to the shape function.
    Returns
    -------
    BaseFunctionHodographicShaping
        Settings object for a cosine contribution to hodographic shaping."""

def hodograph_exponential(exponent: float) -> BaseFunctionHodographicShaping:
    ...

def hodograph_exponential_cosine(exponent: float, frequency: float) -> BaseFunctionHodographicShaping:
    ...

def hodograph_exponential_sine(exponent: float, frequency: float) -> BaseFunctionHodographicShaping:
    ...

def hodograph_power(exponent: float) -> BaseFunctionHodographicShaping:
    ...

def hodograph_power_cosine(exponent: float, frequency: float, scale_factor: float=1.0) -> BaseFunctionHodographicShaping:
    """Function for creating a power cosine function contribution to hodographic trajectory shaping.
    
    Function for creating a power cosine function contribution to hodographic trajectory shaping. For a
    provided exponent :math:`r`, (optional) scale factor :math:`c` and frequency :math:`f`, this adds a contribution :math:`K\\cdot c\\cos(f\\cdot t)\\cdot t^{r}` to the selected
    velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.
    
    
    Parameters
    ----------
    frequency : float
        Frequency of the cosine contribution to the shape function.
    exponent : float
        Exponent of the power function contribution to the shape function.
    scale_factor : float, default = 1.0
        Optional scale factor, which can be used to scale the physical meaning of the free parameter :math:`K`.
    Returns
    -------
    BaseFunctionHodographicShaping
        Settings object for a power cosine function contribution to hodographic shaping."""

def hodograph_power_sine(exponent: float, frequency: float, scale_factor: float=1.0) -> BaseFunctionHodographicShaping:
    """Function for creating a power sine function contribution to hodographic trajectory shaping.
    
    Function for creating a power sine function contribution to hodographic trajectory shaping. For a
    provided exponent :math:`r`, (optional) scale factor :math:`c` and frequency :math:`f`, this adds a contribution :math:`K\\cdot c\\sin(f\\cdot t)\\cdot t^{r}` to the selected
    velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.
    
    
    Parameters
    ----------
    frequency : float
        Frequency of the sine contribution to the shape function.
    exponent : float
        Exponent of the power function contribution to the shape function.
    scale_factor : float, default = 1.0
        Optional scale factor, which can be used to scale the physical meaning of the free parameter :math:`K`.
    Returns
    -------
    BaseFunctionHodographicShaping
        Settings object for a power sine function contribution to hodographic shaping."""

def hodograph_scaled_exponential(exponent: float, scale_factor: float=1.0) -> BaseFunctionHodographicShaping:
    """No documentation found."""

def hodograph_scaled_exponential_cosine(exponent: float, frequency: float, scale_factor: float=1.0) -> BaseFunctionHodographicShaping:
    """No documentation found."""

def hodograph_scaled_exponential_sine(exponent: float, frequency: float, scale_factor: float=1.0) -> BaseFunctionHodographicShaping:
    """No documentation found."""

def hodograph_scaled_power(exponent: float, scale_factor: float=1.0) -> BaseFunctionHodographicShaping:
    """No documentation found."""

def hodograph_scaled_power_cosine(exponent: float, frequency: float, scale_factor: float) -> BaseFunctionHodographicShaping:
    """No documentation found."""

def hodograph_scaled_power_sine(exponent: float, frequency: float, scale_factor: float) -> BaseFunctionHodographicShaping:
    """No documentation found."""

def hodograph_sine(frequency: float) -> BaseFunctionHodographicShaping:
    """Function for creating a sine contribution to hodographic trajectory shaping.
    
    Function for creating a sine contribution to hodographic trajectory shaping. For a
    provided frequency :math:`f`, this adds a contribution :math:`K\\sin(f\\cdot t)` to the selected
    velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.
    
    
    Parameters
    ----------
    frequency : float
        Frequency of the sine contribution to the shape function.
    Returns
    -------
    BaseFunctionHodographicShaping
        Settings object for a cosine contribution to hodographic shaping."""

def recommended_axial_hodograph_functions(time_of_flight: float, number_of_revolutions: int) -> list[BaseFunctionHodographicShaping]:
    """Function for creating the default axial hodograph  ic trajectory shaping functions.
    
    Function for creating the default axial hodographic trajectory shaping functions. This function
    (and its counterparts radial and normal components) provided three shaping functions that have been found in
    literature to work well for this method. For a given time-of-flight :math:`T` and number of revolutions :math:`N`, this function returns a list of
    three shaping functions:
    
    * Cosine term, see :func:`hodograph_cosine` with frequency = :math:`\\frac{2\\pi(N+1/2)}{T}`
    * Power cosine function term, see :func:`hodograph_power_cosine` with  exponent = 3.0, frequency = :math:`\\frac{2\\pi(N+1/2)}{T}`, scale_factor = :math:`1/T`
    * Power sine function term, see :func:`hodograph_power_sine` with  exponent = 3.0, frequency = :math:`\\frac{2\\pi(N+1/2)}{T}`, scale_factor = :math:`1/T`
    
    
    Parameters
    ----------
    time_of_flight : float
        Total time of flight (in seconds) of the trajectory that is to be generated.
    number_of_revolutions : int
        Number of full revolutions around the central body that are to be used.
    Returns
    -------
    list[BaseFunctionHodographicShaping]
        List of default settings object for axial hodographic shaping"""

def recommended_normal_hodograph_functions(time_of_flight: float) -> list[BaseFunctionHodographicShaping]:
    """Function for creating the default normal hodographic trajectory shaping functions.
    
    Function for creating the default normal hodographic trajectory shaping functions. This function
    (and its counterparts radial and axial components) provided three shaping functions that have been found in
    literature to work well for this method. For a given time-of-flight :math:`T`, this function returns a list of
    three shaping functions:
    
    * Constant term, see :func:`hodograph_constant`
    * Power function, see :func:`hodograph_power`, with exponent = 1.0, scale_factor = :math:`1/T`
    * Power function, see :func:`hodograph_power`, with exponent = 2.0, scale_factor = :math:`1/T`
    
    
    Parameters
    ----------
    time_of_flight : float
        Total time of flight (in seconds) of the trajectory that is to be generated.
    Returns
    -------
    list[BaseFunctionHodographicShaping]
        List of default settings object for axial hodographic shaping"""

def recommended_radial_hodograph_functions(time_of_flight: float) -> list[BaseFunctionHodographicShaping]:
    """Function for creating the default radial hodographic trajectory shaping functions.
    
    Function for creating the default radial hodographic trajectory shaping functions. This function
    (and its counterparts normal and axial components) provided three shaping functions that have been found in
    literature to work well for this method. For a given time-of-flight :math:`T`, this function returns a list of
    three shaping functions:
    
    * Constant term, see :func:`hodograph_constant`
    * Power function, see :func:`hodograph_power`, with exponent = 1.0, scale_factor = :math:`1/T`
    * Power function, see :func:`hodograph_power`, with exponent = 2.0, scale_factor = :math:`1/T`
    
    
    Parameters
    ----------
    time_of_flight : float
        Total time of flight (in seconds) of the trajectory that is to be generated.
    Returns
    -------
    list[BaseFunctionHodographicShaping]
        List of default settings object for radial hodographic shaping"""