import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['body_fixed_to_inertial_rotation_matrix', 'inertial_to_body_fixed_rotation_matrix', 'inertial_to_rsw_rotation_matrix', 'inertial_to_tnw_rotation_matrix', 'rotate_state_to_frame', 'rsw_to_inertial_rotation_matrix', 'tnw_to_inertial_rotation_matrix', 'transform_cartesian_state_to_frame']

def body_fixed_to_inertial_rotation_matrix(pole_declination: float, pole_right_ascension: float, pole_meridian: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from body-fixed to inertial frame.
    
    
    Function to compute the rotation matrix from body-fixed to inertial
    frame, using typical pole right ascension (:math:`\\alpha`), pole
    declination (:math:`\\delta`), and prime meridian longitude
    (:math:`W`) angles.
    
    
    Parameters
    ----------
    pole_declination : float
        Declination of body pole in inertial frame (:math:`\\delta`).
    
    pole_right_ascension : float
        Right ascension of body pole in inertial frame (:math:`\\alpha`).
    
    prime_meridian_longitude : float
        Longitude of prime meridian w.r.t. intermediate frame
        (:math:`W`).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from body-fixed to inertial frame.
    
    
    
    
    Notes
    -----
    This definition of a body-fixed orientation is used by,
    for instance, the IAU Working Group on Cartographic Coordinates
    and Rotational Elements. Rotation is performed by a successive z-x-z
    Euler angle rotation (see :cite:t:`Archinal2018`)."""

def inertial_to_body_fixed_rotation_matrix(pole_declination: float, pole_right_ascension: float, prime_meridian_longitude: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from inertial to body-fixed frame.
    
    
    Function to compute the rotation matrix from inertial to body-fixed
    frame, using typical pole right ascension (:math:`\\alpha`), pole
    declination (:math:`\\delta`), and prime meridian longitude
    (:math:`W`) angles.
    
    
    Parameters
    ----------
    pole_declination : float
        Declination of body pole in inertial frame (:math:`\\delta`).
    
    pole_right_ascension : float
        Right ascension of body pole in inertial frame (:math:`\\alpha`).
    
    prime_meridian_longitude : float
        Longitude of prime meridian w.r.t. intermediate frame
        (:math:`W`).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from inertial to body-fixed frame
    
    
    
    Notes
    -----
    This definition of a body-fixed orientation is used by, for
    instance, the IAU Working Group on Cartographic Coordinates and
    Rotational Elements. Rotation is performed by a successive z-x-z
    Euler angle rotation (see :cite:t:`Archinal2018`)."""

def inertial_to_rsw_rotation_matrix(inertial_cartesian_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from inertial to RSW frame.
    
    
    Function to compute the rotation matrix from inertial to RSW frame.
    The RSW frame is defined  by the state of a body w.r.t. to some
    central body. The x-axis of the RSW frame points away from the
    origin, and the y-axis lies in the orbital plane, and is positive
    for in the direction of the velocity vector (but is not colinear
    with the velocity vector, except for circular orbits). The z-axis
    is perpendicular to the orbital plane, and completes the
    right-handed coordinate system.
    
    
    Parameters
    ----------
    inertial_cartesian_state : numpy.ndarray
        Cartesian state, in an inertial frame, for which the rotation
        matrix is to be calculated. Note that the RSW frame is defined
        w.r.t. some central body, and this Cartesian state must be
        defined w.r.t. that central body (e.g. central body at the
        origin).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from inertial to RSW frame."""

def inertial_to_tnw_rotation_matrix(inertial_cartesian_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], n_axis_points_away_from_central_body: bool=True) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from inertial to TNW frame.
    
    
    Function to compute the rotation matrix from inertial to TNW frame.
    The TNW frame is defined by the state of a body w.r.t. to some
    central body. The x-axis of the TNW frame points along the velocity
    vector, and the y-axis lies in the orbital plane, and is positive
    in the direction away from the central body (or positive **towards**
    the central body if the ``n_axis_points_away_from_central_body``
    variable is set to false, see below). The z-axis is perpendicular
    to the orbital plane, and completes the right-handed coordinate
    system.
    
    
    Parameters
    ----------
    inertial_cartesian_state : numpy.ndarray
        Cartesian state, in an inertial frame, for which the rotation
        matrix is to be calculated. Note that the RSW frame is defined
        w.r.t. some central body, and this Cartesian state must be
        defined w.r.t. that central body (e.g. central body at the
        origin).
    
    n_axis_points_away_from_central_body : Boolean
        Boolean (default is ``True``) defining whether the N axis of the
        TNW frame points away from the central body (if ``True``) or
        towards the central body (if ``False``).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from inertial to TNW frame."""

def rotate_state_to_frame(original_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], rotation_matrix: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)], rotation_matrix_time_derivative: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]=...) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]:
    """Rotates a Cartesian state (position and velocity) from one frame :math:`B` to another frame :math:`A`, using the rotation matrix :math:`\\mathbf{R}^{(A/B)}` from frame :math:`B` to :math:`A`, and its time derivative
    :math:`\\dot{\\mathbf{R}}^{(A/B)}`.
    This function computes:
    
    .. math::
       \\mathbf{r}^{(A)}=\\mathbf{R}^{(A/B)}\\mathbf{r}^{(B)}+\\dot{\\mathbf{R}}^{(A/B)}\\mathbf{v}^{(B)}\\\\
       \\mathbf{v}^{(A)}=\\mathbf{R}^{(A/B)}\\mathbf{v}^{(B)}
    
    Parameters
    ----------
    original_state : ndarray[numpy.float64[6, 1]]
        Cartesian state vector :math:`\\mathbf{x}^{(B)}=[\\mathbf{r}^{(B)};\\mathbf{v}^{(B)}]` in frame :math:`B`
    rotation_matrix: ndarray[numpy.float64[3, 3]]
        Rotation matrix :math:`\\mathbf{R}^{(A/B)}` from frame :math:`B` to :math:`A`
    rotation_matrix_time_derivative: ndarray[numpy.float64[3, 3]], default = numpy.zeros((3, 3))
        Time derivative of rotation matrix (:math:`\\dot{\\mathbf{R}}^{(A/B)})` from frame :math:`B` to :math:`A`; default zero indicates that frames :math:`A` and :math:`B` have a constant orientation w.r.t. one another.
    Returns
    -------
    numpy.ndarray
        Input state in frame :math:`B`, rotated to frame :math:`A`"""

def rsw_to_inertial_rotation_matrix(inertial_cartesian_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from RSW to inertial frame.
    
    
    Function to compute the rotation matrix from RSW to inertial. The
    RSW frame is defined  by the state of a body w.r.t. to some central
    body. The x-axis of the RSW frame points away from the origin, and
    the y-axis lies in the orbital plane, and is positive for in the
    direction of the velocity vector (but is not colinear with the
    velocity vector, except for circular orbits). The z-axis is
    perpendicular to the orbital plane, and completes the right-handed
    coordinate system.
    
    
    Parameters
    ----------
    inertial_cartesian_state : numpy.ndarray
        Cartesian state, in an inertial frame, for which the rotation
        matrix is to be calculated. Note that the RSW frame is defined
        w.r.t. some central body, and this Cartesian state must be
        defined w.r.t. that central body (e.g. central body at the
        origin).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from RSW to inertial frame."""

def tnw_to_inertial_rotation_matrix(inertial_cartesian_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(6, 1)], n_axis_points_away_from_central_body: bool=True) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Computes the rotation matrix from TNW to inertial frame.
    
    
    Function to compute the rotation matrix from TNW to inertial frame.
    The TNW frame is defined by the state of a body w.r.t. to some
    central body. The x-axis of the TNW frame points along the velocity
    vector, and the y-axis lies in the orbital plane, and is positive
    in the direction away from the central body (or positive **towards**
    the central body if the ``n_axis_points_away_from_central_body``
    variable is set to false, see below). The z-axis is perpendicular
    to the orbital plane, and completes the right-handed coordinate
    system.
    
    
    Parameters
    ----------
    inertial_cartesian_state : numpy.ndarray
        Cartesian state, in an inertial frame, for which the rotation
        matrix is to be calculated. Note that the TNW frame is defined
        w.r.t. some central body, and this Cartesian state must be
        defined w.r.t. that central body (e.g. central body at the
        origin).
    
    n_axis_points_away_from_central_body : bool
        Boolean (default=``True``) defining whether the N axis of the
        TNW frame points away from the central body (if ``True``) or
        towards the central body (if ``False``).
    
    Returns
    -------
    numpy.ndarray
        Rotation matrix from TNW to inertial frame"""
transform_cartesian_state_to_frame = rotate_state_to_frame