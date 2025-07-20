import typing
__all__ = ['BasicSolidBodyDeformationSettings', 'BodyDeformationSettings', 'basic_solid_body_tidal', 'degree_two_basic_solid_body_tidal', 'iers_2010_solid_body_tidal', 'ocean_tidal', 'pole_tidal']

class BasicSolidBodyDeformationSettings(BodyDeformationSettings):
    """Class for defining model settings for simple tidal solid-body shape deformation.
    
    `BodyDeformationSettings` derived class for simple tidal solid-body shape deformation."""

class BodyDeformationSettings:
    """Base class for providing settings for body shape deformation model.
    
    Functional (base) class for settings of body shape deformation models that require no information in addition to their type.
    Body shape deformation model settings requiring additional information must be defined using an object derived from this class."""

def basic_solid_body_tidal(tide_raising_bodies: list[str], displacement_love_numbers: dict[int, tuple[float, float]], reference_radius: float=...) -> BasicSolidBodyDeformationSettings:
    """Function for creating basic tidal solid-body shape deformation
    
    Function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{m}` and :math:`l_{m}` (with only :math:`m=2,3` currently supported). This function implements equations (7.5) and (7.6) of the `IERS Conventions 2010 <https://iers-conventions.obspm.fr/conventions_material.php>`_.
    
    
    Parameters
    ----------
    tide_raising_bodies : list[ string ]
        List of bodies that raise a tide that induces the shape variation.
    displacement_love_numbers : dict[ int, [float,float] ]
        Dictionary of pairs. The dictionary key the spherical harmonic degree :math:`l` of the tidal deformation (2 or 3 are currently supported). The dictionary value is comprised of a pair of floats representing the :math:`h_{2}` and :math:`l_{2}` deformation Love numbers
    reference_radius : float, default = NaN
        Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
    Returns
    -------
    BasicSolidBodyDeformationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create a settings for degree 2 tidal deformation of the Earth due to the Sun and Moon:
    
    .. code-block:: python
    
      # Create Love numbers
      love_numbers = dict()
      love_numbers[2] = (0.6, 0.08)
    
      # Create tide raising bodies
      tide_raising_bodies = ["Sun", "Moon"]
    
      # Append shape variation settings to existing (default is empty) list
      body_settings.get( "Earth" ).shape_deformation_settings.append(
          environment_setup.shape_deformation.basic_solid_body_tidal(
              tide_raising_bodies, love_numbers ) )"""

def degree_two_basic_solid_body_tidal(tide_raising_bodies: list[str], love_number: float, shida_number: float, reference_radius: float=...) -> BasicSolidBodyDeformationSettings:
    """Function for creating degree 2 basic tidal solid-body shape deformation
    
    Function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{2}` and :math:`l_{2}`. This function implements equations (7.5) of the IERS Conventions 2010, and provides a simplified interface (for degree 2 only) of :func:`basic_solid_body_tidal`.
    
    
    Parameters
    ----------
    tide_raising_bodies : list[ string ]
        List of bodies that raise a tide that induces the shape variation.
    love_number : float
        Value of :math:`h_{2}` deformation Love number`
    shida_number : float
        Value of :math:`l_{2}` deformation Shida number`
    reference_radius : float, default = NaN
        Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
    Returns
    -------
    BasicSolidBodyDeformationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create a settings for degree 2 tidal deformation of the Earth due to the Sun and Moon:
    
    .. code-block:: python
    
      # Create Love numbers
      h2_love_number = 0.6
      l2_shida_number = 0.08
    
      # Create tide raising bodies
      tide_raising_bodies = ["Sun", "Moon"]
    
      # Append shape variation settings to existing (default is empty) list
      body_settings.get( "Earth" ).shape_deformation_settings.append(
          environment_setup.shape_deformation.degree_two_basic_solid_body_tidal(
              tide_raising_bodies, h2_love_number, l2_shida_number ) )"""

def iers_2010_solid_body_tidal() -> BodyDeformationSettings:
    """Function for creating full IERS 2010 shape deformation model
    
    Function for creating full IERS 2010 shape deformation model, computing the tidal shape variation due to the full model defined in
    Section 7.1.1 of the IERS Conventions 2010, implementing Eqs. (7.5)-(7.13), including all terms from Tables 7.3a and 7.3b. At present, none of the input parameters of the model can be varied.
    
    Returns
    -------
    BodyDeformationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` defining the IERS 2010 settings"""

def ocean_tidal(blq_files: list[str]) -> BodyDeformationSettings:
    """No documentation found."""

def pole_tidal() -> BodyDeformationSettings:
    """No documentation found."""