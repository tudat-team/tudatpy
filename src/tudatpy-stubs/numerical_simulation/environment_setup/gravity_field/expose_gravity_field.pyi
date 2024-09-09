import numpy
import typing
__all__ = ['CentralGravityFieldSettings', 'FromFileSphericalHarmonicsGravityFieldSettings', 'GravityFieldSettings', 'GravityFieldType', 'PolyhedronGravityFieldSettings', 'PredefinedSphericalHarmonicsModel', 'SphericalHarmonicsGravityFieldSettings', 'central', 'central_gravity', 'central_spice', 'central_spice_gravity', 'egm96', 'from_file_spherical_harmonic', 'gggrx1200', 'ggm02c', 'ggm02s', 'glgm3150', 'goco05c', 'jgmess160a', 'jgmro120d', 'lpe200', 'polyhedron_from_density', 'polyhedron_from_mu', 'polyhedron_gravity', 'predefined_spherical_harmonic', 'ring_gravity', 'ring_model', 'sh_triaxial_ellipsoid_from_density', 'sh_triaxial_ellipsoid_from_gravitational_parameter', 'shgj180u', 'spherical_harmonic', 'spherical_harmonic_gravity', 'spherical_harmonic_triaxial_body']

class CentralGravityFieldSettings(GravityFieldSettings):
    """`GravityFieldSettings` derived class defining settings of point mass gravity field.
	
	Derived class of `GravityFieldSettings` for central gravity fields, which are defined by a single gravitational parameter.
	"""

    @property
    def gravitational_parameter(self) -> float:
        """
        Gravitational parameter of central gravity field.
        	
        """

    @gravitational_parameter.setter
    def gravitational_parameter(self, arg1: float) -> None:
        ...

class FromFileSphericalHarmonicsGravityFieldSettings(SphericalHarmonicsGravityFieldSettings):
    """
		"""

class GravityFieldSettings:
    """Base class for providing settings for automatic gravity field model creation.
	
	This class is a functional base class for settings of gravity field models that require no information in addition to their type.
	Gravity field model classes requiring additional information must be created using an object derived from this class.
	"""

    @property
    def gravity_field_type(self) -> GravityFieldType:
        """
        Type of gravity field model that is to be created.
        	
        """

class GravityFieldType:
    """Enumeration of gravity field types.
	
	Enumeration of gravity field types supported by tudat.
	
	
	:member polyhedron:
	:member central_gravity:
	:member central_spice_gravity:
	:member spherical_harmonic_gravity:
	:member polyhedron_gravity:
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, GravityFieldType]]
    central_gravity: typing.ClassVar[GravityFieldType]
    central_spice_gravity: typing.ClassVar[GravityFieldType]
    polyhedron_gravity: typing.ClassVar[GravityFieldType]
    ring_gravity: typing.ClassVar[GravityFieldType]
    spherical_harmonic_gravity: typing.ClassVar[GravityFieldType]

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

class PolyhedronGravityFieldSettings(GravityFieldSettings):
    """`GravityFieldSettings` derived class defining settings of a polyhedron gravity field representation.
	
	Derived class of `GravityFieldSettings` for gravity fields, which are defined by a polyhedron gravity field representation.
	"""

    @property
    def associated_reference_frame(self) -> str:
        """
        Identifier for body-fixed reference frame with which the vertices coordinates are associated.
        	
        """

    @associated_reference_frame.setter
    def associated_reference_frame(self, arg1: str) -> None:
        ...

    @property
    def density(self) -> float:
        """
        Density of the polyhedron.
        	
        """

    @density.setter
    def density(self, arg1: float) -> None:
        ...

    @property
    def gravitational_parameter(self) -> float:
        """
        Gravitational parameter of gravity field.
        	
        """

    @gravitational_parameter.setter
    def gravitational_parameter(self, arg1: float) -> None:
        ...

    @property
    def vertices_coordinates(self) -> numpy.ndarray:
        """
        Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
        row per vertex, 3 columns).
        """

    @property
    def vertices_defining_each_facet(self) -> numpy.ndarray:
        """
        Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
        the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
        when seen from the outside of the polyhedron.
        """

class PredefinedSphericalHarmonicsModel:
    """Enumeration of predefined spherical harmonics models.
	
	Enumeration of predefined spherical harmonics models supported by tudat, for which thee coefficient files are automatically available (downloaded from
	`here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_). The directory where these files are stored can be
	extracted using the :func:`~tudatpy.io.get_gravity_models_path` function.
	
	
	:member egm96:
	:member ggm02c:
	:member ggm02s:
	:member goco05c:
	:member glgm3150:
	:member lpe200:
	:member gggrx1200:
	:member jgmro120d:
	:member jgmess160a:
	:member shgj180u:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, PredefinedSphericalHarmonicsModel]]
    egm96: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    gggrx1200: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    ggm02c: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    ggm02s: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    glgm3150: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    goco05c: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    jgmess160a: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    jgmro120d: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    lpe200: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    shgj180u: typing.ClassVar[PredefinedSphericalHarmonicsModel]

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

class SphericalHarmonicsGravityFieldSettings(GravityFieldSettings):
    """`GravityFieldSettings` derived class defining settings of spherical harmonic gravity field representation.
	
	Derived class of `GravityFieldSettings` for gravity fields, which are defined by a spherical harmonic gravity field representation.
	"""
    scaled_mean_moment_of_inertia: float

    @property
    def associated_reference_frame(self) -> str:
        """
        Identifier for body-fixed reference frame with which the coefficients are associated.
        	
        """

    @associated_reference_frame.setter
    def associated_reference_frame(self, arg1: str) -> None:
        ...

    @property
    def create_time_dependent_field(self) -> bool:
        """
        Boolean that denotes whether the field should be created as time-dependent (even if no variations are imposed initially).
        	
        """

    @create_time_dependent_field.setter
    def create_time_dependent_field(self, arg1: bool) -> None:
        ...

    @property
    def gravitational_parameter(self) -> float:
        """
        Gravitational parameter of gravity field.
        	
        """

    @gravitational_parameter.setter
    def gravitational_parameter(self, arg1: float) -> None:
        ...

    @property
    def normalized_cosine_coefficients(self) -> numpy.ndarray:
        """
        Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
        	
        """

    @normalized_cosine_coefficients.setter
    def normalized_cosine_coefficients(self, arg1: numpy.ndarray) -> None:
        ...

    @property
    def normalized_sine_coefficients(self) -> numpy.ndarray:
        """
        Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
        	
        """

    @normalized_sine_coefficients.setter
    def normalized_sine_coefficients(self, arg1: numpy.ndarray) -> None:
        ...

    @property
    def reference_radius(self) -> float:
        """
        Reference radius of spherical harmonic field expansion.
        	
        """

def central(gravitational_parameter: float) -> GravityFieldSettings:
    """Factory function for central gravity field settings object.
	
	Factory function for settings object, defining a point-mass gravity field model with user-defined gravitational parameter :math:`\\mu`. The gravitational potential is the defined as:
	
	.. math::
	   U(\\mathbf{r})=
	
	
	with :math:`\\mathbf{r}` the position vector measured from the body's center of mass.
	
	
	:param gravitational_parameter:
			Gravitational parameter defining the point-mass gravity field.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.CentralGravityFieldSettings` class
	"""

def central_spice() -> GravityFieldSettings:
    """Factory function to create central gravity field settings from Spice settings.
	
	Factory function for settings object, defining a point-mass gravity field model. This function provides the same model as :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.central`), but with gravitational parameter :math:`\\mu` from Spice.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` class of gravity field type ``central_spice``"
	"""

def from_file_spherical_harmonic(file: str, maximum_degree: int, maximum_order: int, associated_reference_frame: str='', gravitational_parameter_index: int=0, reference_radius_index: int=1) -> GravityFieldSettings:
    """Factory function to load a custom spherical harmonics gravity field settings from a file.
	
	Factory function to load a custom spherical harmonics gravity field settings from a file. The file should contain **fully normalized** spherical harmonic coefficients.
	The associated gravitational paramerer and reference radius should be given in m^3/s^2 and m, respectively. The file format should be the same as that used for the files
	in the directories `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_. Specifically, the file should contain
	
	- The first line should be a series of text blocks (typically numerical data). Two of these blocks (by default the first and second one) should be the gravitational parameter and reference radius, respectively. The text block should be separated by spaces, tabs and/or commas
	- Each subsequent line should contain a set of spherical harmonic coefficients (first ordered in ascending order by degree, then in ascending order by order), where the first, second, third and fourth value of the line should be: degree :math:`l`, order :math:`m`, normalized cosine coefficient :math:`\x08ar{C}_{lm}`, normalized sine coefficient :math:`\x08ar{S}_{lm}`. Additional entries (for instance with coefficient uncertainties) are ignored.
	
	
	:param file:
			Full file path and name where th gravity field file is located
	:param maximum_degree:
			Maximum degree of the coefficients that are to be loaded
	:param maximum_order:
			Maximum order of the coefficients that are to be loaded
	:param associated_reference_frame:
			Name of the body-fixed reference frame to which the gravity field is to be fixed. If left empty, this reference frame will automatically be set to the body-fixed frame defined by this body's rotation (see :ref:`\\`\\`rotation_model\\`\\`` for specifying rotation models).
	:param gravitational_parameter_index:
			Index of the values in the file header (first line of file) that contains the gravitational parameter
	:param reference_radius_index:
			Index of the values in the file header (first line of file) that contains the reference radius
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
	"""

def polyhedron_from_density(density: float, vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray, associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> GravityFieldSettings:
    """Factory function for creating a polyhedron gravity field settings object, using the density.
	
	Factory function for settings object, defining a gravity field model through a polyhedron.
	The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
	It represents the frame in which the polyhedron field is defined.
	
	The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
	to Werner and Scheeres [2]_.
	
	This function uses the density to define the gravity field. To instead use the
	gravitational parameter see :func:`~tudatpy.astro.gravitation.polyhedron_from_mu`.
	
	
	:param density:
			Density of the polyhedron.
	
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:param associated_reference_frame:
			Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_constant:
			Newton's gravitational constant G, used to computed the gravitational parameter
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class
	"""

def polyhedron_from_mu(gravitational_parameter: float, vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray, associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> GravityFieldSettings:
    """Factory function for creating a polyhedron gravity field settings object, using the gravitational parameter.
	
	Factory function for settings object, defining a gravity field model through a polyhedron.
	The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
	It represents the frame in which the polyhedron field is defined.
	
	The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
	to Werner and Scheeres [2]_.
	
	This function uses the gravitational parameter to define the gravity field. To instead use the density
	constant see :func:`~tudatpy.astro.gravitation.polyhedron_from_density`.
	
	
	:param gravitational_parameter:
			Gravitational parameter :math:`\\mu` of gravity field.
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:param associated_reference_frame:
			Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_constant:
			Newton's gravitational constant G, used to computed the density
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class
	"""

def predefined_spherical_harmonic(predefined_model: PredefinedSphericalHarmonicsModel, maximum_degree: int=-1) -> GravityFieldSettings:
    """Factory function for spherical harmonics gravity field settings of a predefined model.
	
	Factory function for spherical harmonics gravity field settings of a predefined model
	
	
	:param predefined_model:
			Identified for gravity field model that is to be loaded
	:param maximum_degree:
			Maximum degree and order to which the coefficients are to be loaded. If value is negative, all coefficients for the specified gravity field are loaded
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
	"""

def ring_model(gravitational_parameter: float, ring_radius: float, associated_reference_frame: str, elliptic_integral_s_from_d_and_b: bool) -> GravityFieldSettings:
    ...

def sh_triaxial_ellipsoid_from_density(axis_a: float, axis_b: float, axis_c: float, density: float, maximum_degree: int, maximum_order: int, associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> SphericalHarmonicsGravityFieldSettings:
    """Factory function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the density to define the mass distribution.
	
	Factory function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
	The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
	Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
	The constant mass distribution is defined by the density and gravitational constant (optional).
	The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
	This function implements the models of (see Balmino [1]_).
	
	
	:param axis_a:
			Dimension of largest axis of triaxial ellipsoid.
	:param axis_b:
			Dimension of intermediate axis of triaxial ellipsoid.
	:param axis_c:
			Dimension of smallest axis of triaxial ellipsoid.
	:param density:
			Density of ellipsoid.
	:param maximum_degree:
			Maximum degree of spherical harmonics expansion.
	:param maximum_order:
			Maximum order of spherical harmonics expansion.
	:param associated_reference_frame:
			Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_constant:
			Gravitational constant G of the gravity field.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
	"""

def sh_triaxial_ellipsoid_from_gravitational_parameter(axis_a: float, axis_b: float, axis_c: float, maximum_degree: int, maximum_order: int, associated_reference_frame: str, gravitational_parameter: float) -> SphericalHarmonicsGravityFieldSettings:
    """Factory function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the gravitational parameter to define the mass distribution..
	
	Factory function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
	The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
	Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
	The constant mass distribution is defined by the gravitational parameter.
	The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
	This function implements the models of (see Balmino [1]_).
	
	
	:param axis_a:
			Dimension of largest axis of triaxial ellipsoid.
	:param axis_b:
			Dimension of intermediate axis of triaxial ellipsoid.
	:param axis_c:
			Dimension of smallest axis of triaxial ellipsoid.
	:param maximum_degree:
			Maximum degree of spherical harmonics expansion.
	:param maximum_order:
			Maximum order of spherical harmonics expansion.
	:param associated_reference_frame:
			Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:param gravitational_parameter:
			Gravitational parameter :math:`\\mu` of gravity field.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
	"""

def spherical_harmonic(gravitational_parameter: float, reference_radius: float, normalized_cosine_coefficients: numpy.ndarray, normalized_sine_coefficients: numpy.ndarray, associated_reference_frame: str) -> GravityFieldSettings:
    """Factory function for creating a spherical harmonics gravity field settings object.
	
	Factory function for settings object, defining a gravity field model through spherical harmonic expansion.
	The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
	It represents the frame in which the spherical harmonic field is defined.
	
	The gravitational potential is the defined as:
	
	.. math::
	   U(\\mathbf{r})=\\sum_{l=0}^{l_{max}}\\sum_{m=0}^{l}\\mu\\left(
	
	   heta+\x08ar{S}_{lm}\\sin m  heta
	
	
	with :math:`\\mathbf{r}` the position vector of the evaluation point, measured from the body's center of mass. The angles :math:`\\phi` and :math:`	   heta` are the body-fixed latitude and longitude of the evaluation point, and :math:`\x08ar{P}_{lm}` is the associated Legendre polynomial (at degree/order :math`l/m`).
	
	Note: Spherical harmonic coefficients used for this environment model must *always* be fully normalized.
	To normalize un-normalized spherical harmonic coefficients, see :func:`~tudatpy.astro.gravitation.normalize_spherical_harmonic_coefficients`.
	
	
	:param gravitational_parameter:
			Gravitational parameter :math:`\\mu` of gravity field.
	:param reference_radius:
			Reference radius :math:`R` of spherical harmonic field expansion.
	:param normalized_cosine_coefficients:
			Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\x08ar{C}_{ij}` at degree i and order j.
			As such, note that entry (0,0) of cosine coefficients should be equal to 1.
	
	:param normalized_sine_coefficients:
			Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\x08ar{S}_{ij}`
			at degree i and order j.
	
	:param associated_reference_frame:
			Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
	"""

def spherical_harmonic_triaxial_body(axis_a: float, axis_b: float, axis_c: float, density: float, maximum_degree: int, maximum_order: int, associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> GravityFieldSettings:
    ...
central_gravity: GravityFieldType
central_spice_gravity: GravityFieldType
egm96: PredefinedSphericalHarmonicsModel
gggrx1200: PredefinedSphericalHarmonicsModel
ggm02c: PredefinedSphericalHarmonicsModel
ggm02s: PredefinedSphericalHarmonicsModel
glgm3150: PredefinedSphericalHarmonicsModel
goco05c: PredefinedSphericalHarmonicsModel
jgmess160a: PredefinedSphericalHarmonicsModel
jgmro120d: PredefinedSphericalHarmonicsModel
lpe200: PredefinedSphericalHarmonicsModel
polyhedron_gravity: GravityFieldType
ring_gravity: GravityFieldType
shgj180u: PredefinedSphericalHarmonicsModel
spherical_harmonic_gravity: GravityFieldType