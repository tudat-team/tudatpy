#include <string>

namespace tudatpy {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else {
        return "No documentation found.";
    }

}


    
namespace astro {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else {
        return "No documentation found.";
    }

}


    
namespace frame_conversion {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "inertial_to_rsw_rotation_matrix" && variant==0) {
        return R"(
        
Computes the rotation matrix from inertial to RSW frame.


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
    Rotation matrix from inertial to RSW frame.






    )";



    } else if(name == "rsw_to_inertial_rotation_matrix" && variant==0) {
        return R"(
        
Computes the rotation matrix from RSW to inertial frame.


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
    Rotation matrix from RSW to inertial frame.






    )";



    } else if(name == "inertial_to_tnw_rotation_matrix" && variant==0) {
        return R"(
        
Computes the rotation matrix from inertial to TNW frame.


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
    Rotation matrix from inertial to TNW frame.






    )";



    } else if(name == "tnw_to_inertial_rotation_matrix" && variant==0) {
        return R"(
        
Computes the rotation matrix from TNW to inertial frame.


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
    Rotation matrix from TNW to inertial frame






    )";



    } else if(name == "inertial_to_body_fixed_rotation_matrix" && variant==0) {
        return R"(
        
Computes the rotation matrix from inertial to body-fixed frame.


Function to compute the rotation matrix from inertial to body-fixed
frame, using typical pole right ascension (:math:`\alpha`), pole
declination (:math:`\delta`), and prime meridian longitude
(:math:`W`) angles.


Parameters
----------
pole_declination : float
    Declination of body pole in inertial frame (:math:`\delta`).

pole_right_ascension : float
    Right ascension of body pole in inertial frame (:math:`\alpha`).

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
Euler angle rotation (see Archinal et al. [1]_).




    )";



    } else if(name == "body_fixed_to_inertial_rotation_matrix" && variant==0) {
        return R"(
        
Computes the rotation matrix from body-fixed to inertial frame.


Function to compute the rotation matrix from body-fixed to inertial
frame, using typical pole right ascension (:math:`\alpha`), pole
declination (:math:`\delta`), and prime meridian longitude
(:math:`W`) angles.


Parameters
----------
pole_declination : float
    Declination of body pole in inertial frame (:math:`\delta`).

pole_right_ascension : float
    Right ascension of body pole in inertial frame (:math:`\alpha`).

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
Euler angle rotation (see Archinal et al. [1]_).




    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace element_conversion {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "cartesian_to_keplerian" && variant==0) {
        return R"(
        
Convert Cartesian to Keplerian elements.

.. note:: See module level documentation for the standard ordering
          convention of Keplerian elements used.


Parameters
----------
cartesian_elements : numpy.ndarray
    Cartesian state that is to be converted to Keplerian elements
gravitational_parameter : float
    Gravitational parameter of central body used for conversion
Returns
-------
numpy.ndarray
    Keplerian elements, as computed from Cartesian element input.






    )";



    } else if(name == "keplerian_to_cartesian" && variant==0) {
        return R"(
        
Convert Keplerian elements to Cartesian.

.. note:: See module level documentation for the standard ordering
          convention of Keplerian elements used.


Parameters
----------
keplerian_elements : numpy.ndarray
    Keplerian state that is to be converted to Cartesian elements
gravitational_parameter : float
    Gravitational parameter of central body used for conversion
Returns
-------
numpy.ndarray
    Cartesian elements, as computed from Keplerian element input.






    )";



    } else if(name == "keplerian_to_cartesian_elementwise" && variant==0) {
        return R"(
        
Convert Keplerian elements to Cartesian, with elementwise input.

.. note:: The final Keplerian element is always the true anomaly.


Parameters
----------
semi_major_axis : float
    Semi-major axis (except if eccentricity = 1.0, then represents semi-latus rectum)
eccentricity : float
    Eccentricity
inclination : float
    Inclination
argument_of_periapsis : float
    Argument of periapsis
longitude_of_ascending_node : float
    Longitude of ascending node
true_anomaly : float
    True anomaly
gravitational_parameter : float
    Gravitational parameter of central body used for conversion
Returns
-------
numpy.ndarray
    Cartesian elements, as computed from Keplerian element input.






    )";



    } else if(name == "mean_to_true_anomaly" && variant==0) {
        return R"(
        
Convert mean to true anomaly.

Convert the mean anomaly of the orbit to its true anomaly. This conversion first converts mean to eccentric anomaly
(hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1), and subsequently to true anomaly.


Parameters
----------
eccentricity : float
    Value of the orbital eccentricity
mean_anomaly : float
    Value of the mean anomaly
use_default_initial_guess : bool, default = True
    Boolean to determine whether the user-defined initial guess (for mean-to-eccentric anomaly conversion) is used, or an automatically generated one.
non_default_initial_guess : float, default = NaN
    User-defined initial guess for mean-to-eccentric anomaly conversion, to be used only if ``use_default_initial_guess`` is set to ``False``.
root_finder : RootFinder, default = None
    User-defined root finder, overriding default root-finding algorithm for mean-to-eccentric anomaly conversion (default is used if this input is left empty)
Returns
-------
float
    Value of the true anomaly






    )";



    } else if(name == "true_to_mean_anomaly" && variant==0) {
        return R"(
        
Convert true to mean anomaly.

Convert the true anomaly of the orbit to its mean anomaly. This conversion first converts true to eccentric anomaly
(hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1),
and subsequently to mean anomaly.


Parameters
----------
eccentricity : float
    Value of the orbital eccentricity
true_anomaly : float
    Value of the true anomaly
Returns
-------
float
    Value of the mean anomaly






    )";



    } else if(name == "true_to_eccentric_anomaly" && variant==0) {
        return R"(
        
Convert true to eccentric anomaly.


Parameters
----------
eccentricity : float
    Value of the orbital eccentricity
true_anomaly : float
    Value of the true anomaly
Returns
-------
float
    Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1






    )";



    } else if(name == "eccentric_to_true_anomaly" && variant==0) {
        return R"(
        
Convert eccentric to true anomaly.


Parameters
----------
eccentric_anomaly : float
    Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1
eccentricity : float
    Value of the orbital eccentricity
Returns
-------
float
    Value of the true anomaly






    )";



    } else if(name == "eccentric_to_mean_anomaly" && variant==0) {
        return R"(
        
Convert eccentric to mean anomaly.


Parameters
----------
eccentric_anomaly : float
    Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1
eccentricity : float
    Value of the orbital eccentricity
Returns
-------
float
    Value of the mean anomaly






    )";



    } else if(name == "mean_to_eccentric_anomaly" && variant==0) {
        return R"(
        
Convert mean to eccentric anomaly.


Parameters
----------
eccentricity : float
    Value of the orbital eccentricity
mean_anomaly : float
    Value of the mean anomaly
use_default_initial_guess : bool, default = True
    Boolean to determine whether the user-defined initial guess is used for conversion, or an automatically generated one.
non_default_initial_guess : float, default = NaN
    User-defined initial guess for conversion, to be used only if ``use_default_initial_guess`` is set to ``False``.
root_finder : RootFinder, default = None
    User-defined root finder, overriding default root-finding algorithm for conversion (default is used if this input is left empty)
Returns
-------
float
    Value of the eccentric anomaly






    )";



    } else if(name == "elapsed_time_to_delta_mean_anomaly" && variant==0) {
        return R"(
        
Convert elapsed time to the corresponding change in mean anomaly along a Keplerian orbit.


Parameters
----------
elapsed_time : float
    Elapsed time (in seconds)
gravitational_parameter : float
    Gravitational parameter of central body
semi_major_axis : float
    Semi-major axis of orbit
Returns
-------
float
    Total change in mean anomaly along the Kepler orbit, accumulated in the provided time.






    )";



    } else if(name == "delta_mean_anomaly_to_elapsed_time" && variant==0) {
        return R"(
        
Convert change in mean anomaly along a Keplerian orbit to the corresponding elapsed time.


Parameters
----------
mean_anomaly_change : float
    Total change in mean anomaly along the Kepler orbit
gravitational_parameter : float
    Gravitational parameter of central body
semi_major_axis : float
    Semi-major axis of orbit
Returns
-------
float
    Time required for the provided mean anomaly change to be accumulated






    )";



    } else if(name == "mean_motion_to_semi_major_axis" && variant==0) {
        return R"(
        
Convert mean motion to corresponding semi-major axis (in a Keplerian orbit).


Parameters
----------
mean_motion : float
    Orbital mean motion
gravitational_parameter : float
    Gravitational parameter of central body
Returns
-------
float
    Semi-major axis corresponding to mean motion






    )";



    } else if(name == "semi_major_axis_to_mean_motion" && variant==0) {
        return R"(
        
Convert semi-major axis to corresponding mean motion (along a Keplerian orbit).


Parameters
----------
semi_major_axis : float
    Semi-major axis of orbit
gravitational_parameter : float
    Gravitational parameter of central body
Returns
-------
float
    Semi-major axis corresponding to mean motion






    )";



    } else if(name == "keplerian_to_mee_manual_singularity" && variant==0) {
        return R"(
        
Convert Keplerian to Modified equinoctial elements.

Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
element :math:`I` is to be provided manually for this function

.. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


Parameters
----------
keplerian_elements : numpy.ndarray
    Keplerian elements that are to be converted to Modified equinoctial elements
singularity_at_zero_inclination : bool
    Singularity at 0 degrees inclination if ``True``, 180 degrees if ``False``
Returns
-------
numpy.ndarray
    Modified equinoctial elements, as computed from Keplerian element input.






    )";



    } else if(name == "keplerian_to_mee" && variant==0) {
        return R"(
        
Convert Keplerian to Modified equinoctial elements.

Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)

.. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


Parameters
----------
keplerian_elements : numpy.ndarray
    Keplerian elements that are to be converted to Modified equinoctial elements
Returns
-------
numpy.ndarray
    Modified equinoctial elements, as computed from Keplerian element input (with element :math:`I` defined by :func:`flip_mee_singularity`).






    )";



    } else if(name == "flip_mee_singularity" && variant==0) {
        return R"(
        
Function to determine 'optimal' location of the singularity-flipping modified equinoctial element.

Function to determine 'optimal' location of the singularity-flipping modified equinoctial element :math:`I`, if orbit inclination is less than
90 degrees, it puts the singularity at 180 degrees, if it is larger than 90 degrees, it puts it at 0 degrees.


Parameters
----------
keplerian_elements : numpy.ndarray
    Keplerian elements that are to be converted to Modified equinoctial elements
Returns
-------
bool
    Singularity at 0 degrees inclination if false, 180 degrees if true






    )";



    } else if(name == "mee_to_keplerian" && variant==0) {
        return R"(
        
Convert Modified equinoctial to Keplerian elements.

Modified equinoctial elements to Keplerian (without intermediate step to Cartesian elements).

.. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


Parameters
----------
modified_equinoctial_elements : numpy.ndarray
    Modified equinoctial elements that are to be converted to Keplerian elements
singularity_at_zero_inclination : bool
    Singularity at 0 degrees inclination if false, 180 degrees if true
Returns
-------
numpy.ndarray
    Keplerian elements, as computed from Modified equinoctial element input.






    )";



    } else if(name == "cartesian_to_mee" && variant==0) {
        return R"(
        
Convert Cartesian to Modified equinoctial elements.

Convert cartesian to Modified equinoctial elements. The singularity-flipping
element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)

.. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


Parameters
----------
cartesian_elements : numpy.ndarray
    Cartesian elements that are to be converted to Modified equinoctial elements
gravitational_parameter : float
    Gravitational parameter of central body
Returns
-------
numpy.ndarray
    Modified equinoctial elements, as computed from Cartesian element input.






    )";



    } else if(name == "cartesian_to_mee_manual_singularity" && variant==0) {
        return R"(
        
Convert Cartesian to Modified equinoctial elements.

Convert cartesian to Modified equinoctial elements. The singularity-flipping
element :math:`I` is to be provided manually for this function

.. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


Parameters
----------
cartesian_elements : numpy.ndarray
    Cartesian elements that are to be converted to Modified equinoctial elements
gravitational_parameter : float
    Gravitational parameter of central body
singularity_at_zero_inclination : bool
    Singularity at 0 degrees inclination if false, 180 degrees if true
Returns
-------
numpy.ndarray
    Modified equinoctial elements, as computed from Cartesian element input.






    )";



    } else if(name == "mee_to_cartesian" && variant==0) {
        return R"(
        
Convert Modified equinoctial to Cartesian elements.


.. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


Parameters
----------
modified_equinoctial_elements : numpy.ndarray
    Modified equinoctial elements that are to be converted to Cartesian elements
gravitational_parameter : float
    Gravitational parameter of central body
singularity_at_zero_inclination : bool
    Singularity at 0 degrees inclination if false, 180 degrees if true
Returns
-------
numpy.ndarray
    Cartesian elements, as computed from Modified equinoctial element input.






    )";



    } else if(name == "quaternion_entries_to_rotation_matrix" && variant==0) {
        return R"(
        
Converts an array of four quaternion elements to the equivalent rotation matrix.

Function to convert an array of four quaternion elements to the equivalent rotation matrix. These quaternion elements
are for instance used when propagating rotational dynamics in Tudat, and this function can be used to convert the
numerical results to a usable rotation matrix. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html?highlight=rotational%20states#rotational-states>`_ for more details.


Parameters
----------
quaternion_entries : numpy.ndarray
    Quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
Returns
-------
numpy.ndarray
    Rotation matrix defining the equivalent rotation.






    )";



    } else if(name == "rotation_matrix_to_quaternion_entries" && variant==0) {
        return R"(
        
Converts a rotation matrix to the equivalent array of four quaternion elements.

Inverse function of :func:`quaternion_entries_to_rotation_matrix`.


Parameters
----------
rotation_matrix : numpy.ndarray
    Rotation matrix
Returns
-------
numpy.ndarray
    Equivalent quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_






    )";



    } else if(name == "cartesian_to_spherical" && variant==0) {
        return R"(
        
Convert Cartesian to spherical elements.

.. note:: See module level documentation for the standard ordering  convention of spherical state elements used.


Parameters
----------
cartesian_elements : numpy.ndarray
    Cartesian state that is to be converted to spherical elements
Returns
-------
numpy.ndarray
    Spherical elements, as computed from Cartesian element input.






    )";



    } else if(name == "spherical_to_cartesian" && variant==0) {
        return R"(
        
Convert spherical elements to Cartesian.

.. note:: See module level documentation for the standard ordering convention of spherical state elements used.


Parameters
----------
spherical_elements : numpy.ndarray
    Spherical state that is to be converted to Cartesian elements
Returns
-------
numpy.ndarray
    Cartesian elements, as computed from spherical element input.






    )";



    } else if(name == "spherical_to_cartesian_elementwise" && variant==0) {
        return R"(
        
Convert Spherical elements to Cartesian, with elementwise input.


Parameters
----------
radial_distance : float
    Distance from origin of central body
latitude : float
    Central body-fixed latitude
longitude : float
    Central body-fixed longitude
speed : float
    Central body-fixed speed (norm of velocity vector). Note that this is *not* the norm of the inertial velocity
flight_path_angle : float
    Flight-path angle (of central body-fixed velocity vector)
heading_angle : float
    Heading angle (of central body-fixed velocity vector)
Returns
-------
numpy.ndarray
    Cartesian elements, as computed from spherical element input.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace time_conversion {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "DateTime") {
         return R"(

        Class to store a calendar date and time of day, with high resolution.

        Class to store a calendar date and time of day, with high resolution compared to Python datetime.datetime. This class
        stores the seconds as a ``long double`` variable in the C++ implementation, corresponding to about
        16 or 19 digits of precision (depending on the compiler used). In either case, this will be sufficient for sub-femtosecond
        resolution. In addition, this class allows easy conversion to typical time representations in astrodynamics (seconds since J2000,
        Julian day, and modified Julian day).





     )";


    } else if(name == "DateTime.year") {
         return R"(

        Calendar year


        :type: int
     )";


    } else if(name == "DateTime.month") {
         return R"(

        Calendar month (value must be 1-12)     


        :type: int
     )";


    } else if(name == "DateTime.day") {
         return R"(

        Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month       


        :type: int
     )";


    } else if(name == "DateTime.hour") {
         return R"(

        Full hours into the current day (value must be 0-23)


        :type: int
     )";


    } else if(name == "DateTime.minute") {
         return R"(

        Full minutes into the current hour (value must be 0-59)


        :type: int
     )";


    } else if(name == "DateTime.seconds") {
         return R"(

        Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.     


        :type: float
     )";




    } else if(name == "DateTime.ctor" && variant==0) {
            return R"(

        Create a date time from decomposed date (Gregorian calendar) and time (UTC).


        Parameters
        ----------
        year : int
            Calendar year.
        month : int
            Calendar month.
        day : int
            Calendar day.
        hour : int, default = 12
            Hours into the day current day (0-23).
        minute : int, default = 0
            Minutes in the current hour (0-59).
        seconds : float, default = 0.0
            Seconds in the current minute.




    )";



    } else if(name == "DateTime.iso_string" && variant==0) {
            return R"(

        Function to get the ISO-compatible string.


        Function to get the current date and time as an ISO-compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..") where the seconds may be provided with any number of digits. The 'T' entry separating the date from the time may be omitted by setting the ``add_T`` parameter to false


        Parameters
        ----------
        add_T : Bool
            Boolean denoting whether to use a 'T' or a blank space to separate the date from the time   

        number_of_digits_seconds : int, default = 15
            Number of digits to use after the decimal separator (trailing zeros will be truncated)

        Returns
        -------
        str
            ISO-compatible string representing the date and time





    )";



    } else if(name == "DateTime.epoch" && variant==0) {
            return R"(

        Function to get the epoch in seconds since J2000 for the current date and time 


        Returns
        -------
        float
            Current epoch in seconds since J2000





    )";



    } else if(name == "DateTime.julian_day" && variant==0) {
            return R"(

        Function to get the epoch as Julian day for the current date and time 


        Returns
        -------
        float
            Current Julian day





    )";



    } else if(name == "DateTime.modified_julian_day" && variant==0) {
            return R"(

        Function to get the epoch as modified Julian day for the current date and time 


        Returns
        -------
        float
            Current modified Julian day





    )";



    } else if(name == "DateTime.day_of_year" && variant==0) {
            return R"(

        Function to get the day number in the current year


        Returns
        -------
        int
            Day number in the current year





    )";





    } else if(name == "datetime_to_tudat" && variant==0) {
        return R"(
        
Function to convert a Python datetime.datetime object to a Tudat :class:`DateTime` object. The Tudat-native alternative has the advantage of providing sub-femtosecond resolution, as opposed to the microsecond resolution of the Python version


Parameters
----------
datetime : datetime.datetime
    Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
Returns
-------
DateTime
    DateTime object defined in Tudat






    )";



    } else if(name == "datetime_to_python" && variant==0) {
        return R"(
        
Function to convert a Tudat :class:`DateTime` object to a Python datetime.datetime object. This is the inverse of the :func:`datetime_to_tudat` function


Parameters
----------
datetime : DateTime
    Tudat-native Datetime object. Both the date and the time (hour, minutes, and seconds), can be specified, up to sub-femtosecond resolution.
Returns
-------
datetime.datetime
    Datetime object, using the Python datetime library






    )";



    } else if(name == "add_seconds_to_datetime" && variant==0) {
        return R"(
        
Function to create a new Tudat :class:`DateTime` object by adding a number of seconds to an existing Tudat :class:`DateTime` object


Parameters
----------
datetime : DateTime
    Tudat-native Datetime object to which a number of seconds are to be added
seconds_to_add : float
    Number of seconds to add
Returns
-------
DateTime
    Tudat-native Datetime object created by adding the given number of seconds to the original DateTime






    )";



    } else if(name == "add_days_to_datetime" && variant==0) {
        return R"(
        
Function to create a new Tudat :class:`DateTime` object by adding a number of days (86400 seconds) to an existing Tudat :class:`DateTime` object


Parameters
----------
datetime : DateTime
    Tudat-native Datetime object to which a number of days are to be added
days_to_add : float
    Number of days to add
Returns
-------
DateTime
    Tudat-native Datetime object created by adding the given number of days to the original DateTime






    )";



    } else if(name == "calendar_date_to_julian_day" && variant==0) {
        return R"(
        
Convert a calendar date to Julian days.


Parameters
----------
calendar_date : datetime.datetime
    Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified, up to millisecond resolution.
Returns
-------
float
    Julian day number (days since noon January 1st 4713 BC.)





Examples
--------
In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days.

.. code-block:: python

  # Define the calendar date using datetime
  calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
  # Convert the calendar date to Julian days since January 1st 4713 BC 
  julian_date = time_conversion.calendar_date_to_julian_day(calendar_date)
  # Print the converted output
  print(julian_date)  # prints 2459721.0782523146


    )";



    } else if(name == "calendar_date_to_days_since_epoch" && variant==0) {
        return R"(
        
Convert a calendar date to Julian days since a given epoch.


Parameters
----------
calendar_date : datetime.datetime
    Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored.
days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
    Reference epoch (in days) since when the Julian days have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0) corresponding to the 1st of January 2000.
Returns
-------
float
    Date in Julian days since the given epoch.





Examples
--------
In this example, the calendar date of the 21st of May 2022 at 13:52 and 41 seconds is converted to Julian days since J2000 (the 1st of January 2000).

.. code-block:: python

  # Define the calendar date using datetime
  calendar_date = datetime.datetime(2022, 5, 21, 13, 52, 41)
  # Convert the calendar date to Julian days since J2000
  julian_date = time_conversion.calendar_date_to_days_since_epoch(calendar_date)
  # Print the converted output
  print(julian_date)  # prints 8176.07825231459


    )";



    } else if(name == "julian_day_to_calendar_date" && variant==0) {
        return R"(
        
Convert Julian days to a calendar date.

Inverse function of :func:`calendar_date_to_julian_day`.

Parameters
----------
julian_day : float
    Date in Julian days since January 1st 4713 BC.
Returns
-------
datetime.datetime
    Datetime object, using the Python datetime library, containing the date and time corresponding to the Julian date input.





Examples
--------
In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to a calendar date.

.. code-block:: python

  # Define the Julian date in days since January 1st 4713 BC
  julian_date = 2459721.0783
  # Convert the Julian date to a calendar date
  calendar_date = time_conversion.julian_day_to_calendar_date(julian_date)
  # Print the converted output
  print(calendar_date)  # prints datetime.datetime(2022, 5, 21, 13, 52, 45)


    )";



    } else if(name == "julian_day_to_seconds_since_epoch" && variant==0) {
        return R"(
        
Convert Julian days to seconds since a given epoch.


Parameters
----------
julian_day : float
    Date in Julian days since January 1st 4713 BC.
days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
    Reference epoch (in days since January 1st 4713 BC) since when the number of seconds have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
Returns
-------
float
    Seconds since the Julian date and the given epoch.





Examples
--------
In this example, the Julian date `2459721.0783` (in days since January 1st 4713 BC), is converted to seconds since J2000 (January 1st 2000).

.. code-block:: python

  # Define the Julian date in days since January 1st 4713 BC
  julian_date = 2459721.0783
  # Convert the Julian date to the number of seconds since J2000
  seconds_since_J2000 = time_conversion.julian_day_to_seconds_since_epoch(julian_date)
  # Print the converted output
  print(seconds_since_J2000)  # prints 706413165.1200145
  


    )";



    } else if(name == "seconds_since_epoch_to_julian_day" && variant==0) {
        return R"(
        
Convert seconds since a given reference epoch to a Julian day.

Inverse function of :func:`julian_day_to_seconds_since_epoch`.

Parameters
----------
seconds_since_epoch : float
    Seconds since ``days_since_julian_day_zero`` which are to be converted to date in Julian days.
days_since_julian_day_zero : float, default = constants.JULIAN_DAY_ON_J2000
    Reference epoch (in days since January 1st 4713 BC) since when the number of seconds have to be counted. By default, set to `constants.JULIAN_DAY_ON_J2000` (2451545.0), corresponding to the 1st of January 2000.
Returns
-------
float
    Date in Julian days since January 1st 4713 BC, as computed from the input parameters





Examples
--------
In this example, an amount of seconds since J2000 (January 1st 2000) is converted to the Julian date (in days since January 1st 4713 BC).

.. code-block:: python
  # Define the amount of seconds since January 1st 2000
  seconds_since_J2000 = 706413165.1200145
  # Convert the amount of seconds since J2000 to the Julian date
  julian_date = time_conversion.seconds_since_epoch_to_julian_day(seconds_since_J2000)
  # Print the converted output
  print(julian_date)  # prints 2459721.0783


    )";



    } else if(name == "seconds_since_epoch_to_julian_years_since_epoch" && variant==0) {
        return R"(
        
Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch.

Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch. This is equivalent to converting a time interval in seconds to Julian years

Parameters
----------
seconds_since_epoch : float
    Seconds elapsed since a given (unspecified) epoch.
Returns
-------
float
    Julian years since the specified epoch.

    Since this is a float, not a integer, meaning that the fraction of the year is also included.






Examples
--------
In this example, `706413165.12` seconds since a given epoch are converted to Julian years since the same epoch.

.. code-block:: python

  # Define the number of seconds elapsed
  seconds_since_epoch = 706413165.12
  # Convert the number of seconds to Julian years
  julian_years = time_conversion.seconds_since_epoch_to_julian_years_since_epoch(seconds_since_epoch)
  # Print the converted output
  print(julian_years)  # prints 22.38488240930869


    )";



    } else if(name == "seconds_since_epoch_to_julian_centuries_since_epoch" && variant==0) {
        return R"(
        
Convert the number of seconds since a given (unspecified) epoch to Julian centuries since the same epoch.

Convert the number of seconds since a given (unspecified) epoch to Julian years since the same epoch. This is equivalent to converting a time interval in seconds to Julian centuries

Parameters
----------
seconds_since_epoch : float
    Seconds elapsed since a given (unspecified) epoch.
Returns
-------
float
    Julian centuries since the specified epoch.

    Since this is a float, not a integer, meaning that the fraction of the century is also included.






Examples
--------
In this example, `706413165.12` seconds since a given epoch are converted to Julian centuries since the same epoch.

.. code-block:: python

  # Define the number of seconds elapsed
  seconds_since_epoch = 706413165.12
  # Convert the number of seconds to Julian centuries
  julian_centuries = time_conversion.seconds_since_epoch_to_julian_centuries_since_epoch(seconds_since_epoch)
  # Print the converted output
  print(julian_centuries)  # prints 0.2238488240930869


    )";



    } else if(name == "julian_day_to_modified_julian_day" && variant==0) {
        return R"(
        
Convert a Julian day to a Modified Julian day.


Parameters
----------
julian_day : float
    Date in Julian days (number of days since January 1st 4713 BC).
Returns
-------
float
    Date in modified Julian days (number of days since November 17th 1858).





Examples
--------
In this example, the Julian date `2451545.0` (J2000) is converted to a modified Julian date.

.. code-block:: python

  # Convert from Julian Days to Modified Julian Days
  MJD = time_conversion.julian_day_to_modified_julian_day(constants.JULIAN_DAY_ON_J2000)
  # Print the converted output
  print(MJD)  # prints 51544.5


    )";



    } else if(name == "modified_julian_day_to_julian_day" && variant==0) {
        return R"(
        
Convert a Modified Julian day to a Julian day.

Inverse function of :func:`julian_day_to_modified_julian_day`.

Parameters
----------
modified_julian_day : float
    Date in modified Julian days (number of days since November 17th 1858).
Returns
-------
float
    Date in Julian days (number of days since January 1st 4713 BC).





Examples
--------
In this example, the Modified Julian date `51544.5` ( corresponding to J2000) is converted to a modified Julian date.

.. code-block:: python

  # Define J2000 in Modified Julian Days
  J2000_MJD = 51544.5
  # Convert from Modified Julian Days to Julian Days
  J2000 = time_conversion.modified_julian_day_to_julian_day(J2000_MJD)
  # Print the converted output
  print(J2000)  # prints 2451545.0


    )";



    } else if(name == "calendar_date_to_day_of_year" && variant==0) {
        return R"(
        
Determine the number of full days that have passed in the year of a given calendar date.


Parameters
----------
calendar_date : datetime.datetime
    Datetime object, using the Python datetime library. Both the date and the time (hour, minutes, and seconds), can be specified. Milliseconds are ignored."
Returns
-------
int
    Number of full days that have passed in the year at the given calendar date.





Examples
--------
In this example, the number of days that have passed in 2020 when the date is the 2nd of May is computed.

.. code-block:: python

  # Define the 2nd of May 2020
  date = datetime.datetime(2020, 5, 2)
  # Compute the number of full days that have passed in 2020 when at the given date
  days_passed = time_conversion.calendar_date_to_day_of_year(date)
  # Print the converted output
  print(J2000)  # prints 122.0


    )";



    } else if(name == "year_and_days_in_year_to_calendar_date" && variant==0) {
        return R"(
        
Create the calendar date from the year and the number of days in the year.

Can be seen as the inverse function of :func:`calendar_date_to_day_of_year`.

Parameters
----------
year : int
    Calendar year.
days_in_year : int
    Number of days that have passed in the year.
Returns
-------
datetime.datetime
    Corresponding calendar date as a datetime object, using the Python datetime library. .





Examples
--------
In this example, the calendar date corresponding to when 122 days have passed in 2020 is computed.

.. code-block:: python

  # Compute the calendar date when 122 days have passed in 2020
  date = time_conversion.year_and_days_in_year_to_calendar_date(2020, 122)
  # Print the converted output
  print(J2000)  # prints datetime.datetime(2020, 5, 2, 0, 0)


    )";



    } else if(name == "calculate_seconds_in_current_julian_day" && variant==0) {
        return R"(
        
Determine the number of seconds that have elapsed in the given Julian day.


Parameters
----------
julian_day : float
    Date in Julian days (number of days since January 1st 4713 BC).
Returns
-------
float
    Number of seconds that have passed in the given Julian day.





Examples
--------
In this example, the number of seconds that have elapsed at the Julian day `2451545.2` is computed.

.. code-block:: python

  # Compute the number of seconds that have passed in the given Julian day
  seconds_passed = time_conversion.calculate_seconds_in_current_julian_day(constants.JULIAN_DAY_ON_J2000)
  # Print the converted output
  print(seconds_passed)  # prints 43200.0


    )";



    } else if(name == "is_leap_year" && variant==0) {
        return R"(
        
Assess wether a year is a leap year or not.


Parameters
----------
year : int
    Calendar year.
Returns
-------
bool
    A value of True means that the year is a leap year.





Examples
--------
In this example, the first list should contains only `True`, and the second `False`, since the first list uses leap years and the second does not.

.. code-block:: python

  # Check known leap years
  leap_years = [time_conversion.is_leap_year(year) for year in [2020, 2016, 2000, 2400]]
  # Print the converted output
  print(leap_years)  # prints [True, True, True, True]
  # Check known non-leap years
  non_leap_years = [time_conversion.is_leap_year(year) for year in [2021, 2022, 2100, 2001]]
  # Print the converted output
  print(non_leap_years)  # prints [False, False, False, False]


    )";



    } else if(name == "get_days_in_month" && variant==0) {
        return R"(
        
Get the number of days in the month of a given year.


Parameters
----------
month : int
    Calendar month.
year : int
    Calendar year.
Returns
-------
int
    Number of days in the month of the given year.





Examples
--------
In this example, the number of days in February for both 2021 and 2020 are computed.

.. code-block:: python

  # Check the number of days in February 2021
  days_feb_2021 = time_conversion.get_days_in_month(2, 2021)
  # Print the converted output
  print(days_feb_2021)  # prints 28
  # Check the number of days in February 2022
  days_feb_2020 = time_conversion.get_days_in_month(2, 2020)
  # Print the converted output
  print(days_feb_2020)  # prints 29


    )";



    } else if(name == "TCB_to_TDB" && variant==0) {
        return R"(
        
Convert time from the TCB scale to the TDB scale.

The TCB scale is the Barycentric Coordinate Time, and the TDB scale is the Barycentric Dynamical Time.

Parameters
----------
TCB_time : float
    Time in seconds since J2000, in the TCB time scale.
Returns
-------
float
    Time in seconds since J2000, in the TDB time scale.





Examples
--------
In this example, the calendar date of the 17th of February 2022, at 15:41 and 2 seconds is first converted to Julian seconds since J2000.
Then, this date and time is converted from the TCB scale to the TDB scale.

.. code-block:: python

  # Define the date and time
  date = datetime.datetime(2022, 2, 17, 15, 41, 2)
  # Convert it in Julian days since J2000
  date_J2000 = time_conversion.calendar_date_to_julian_day(date)
  # Convert it in Julian seconds since J2000
  date_J2000_sec = time_conversion.julian_day_to_seconds_since_epoch(date_J2000)
  # Check the date from the TCB scale to the TDB scale
  date_TDB_scale = time_conversion.TCB_to_TDB(date_J2000_sec)
  # Print the converted output
  print(date_TDB_scale)  # prints 698384439.9176273


    )";



    } else if(name == "TDB_to_TCB" && variant==0) {
        return R"(
        
Convert time from the TBD scale to the TCB scale.

The TDB scale is the Barycentric Dynamical Time, and the TCB scale is the Barycentric Coordinate Time.

Inverse function of :func:`TCB_to_TDB`.


Parameters
----------
TDB_time : float
    Time in seconds since J2000, in the TDB time scale.
Returns
-------
float
    Time in seconds since J2000, in the TCB time scale.






    )";



    } else if(name == "TCG_to_TT" && variant==0) {
        return R"(
        
Convert time from the TCG scale to the TT scale.

The TCG scale is the Geocentric Coordinate Time, and the TT scale is the Terrestrial Time.

Parameters
----------
TCG_time : float
    Time in seconds since J2000, in the TCG time scale.
Returns
-------
float
    Time in seconds since J2000, in the TT time scale.






    )";



    } else if(name == "TT_to_TCG" && variant==0) {
        return R"(
        
Convert time from the TT scale to the TCG scale.

The TT scale is the Terrestrial Time, and the TCG scale is the Geocentric Coordinate Time.

Inverse function of :func:`TCG_to_TT`.


Parameters
----------
TT_time : float
    Time in seconds since J2000, in the TT time scale.
Returns
-------
float
    Time in seconds since J2000, in the TCG time scale.






    )";



    } else if(name == "TAI_to_TT" && variant==0) {
        return R"(
        
Convert time from the TAI scale to the TT scale.

The TAI scale is the International Atomic Time, and the TT scale is the Terrestrial Time.

Parameters
----------
TAI_time : float
    Time in seconds since J2000, in the TAI time scale.
Returns
-------
float
    Time in seconds since J2000, in the TT time scale.






    )";



    } else if(name == "TT_to_TAI" && variant==0) {
        return R"(
        
Convert time from the TT scale to the TAI scale.

The TT scale is the Terrestrial Time, and the TAI scale is the International Atomic Time.

Inverse function of :func:`TAI_to_TT`.


Parameters
----------
TT_time : float
    Time in seconds since J2000, in the TT time scale.
Returns
-------
float
    Time in seconds since J2000, in the TAI time scale.






    )";



    } else if(name == "TT_to_TDB_approximate" && variant==0) {
        return R"(
        
Approximately convert time from the TT scale to the TDB scale.

The TT scale is the Terrestrial Time, and the TDB scale is the Barycentric Dynamical Time.

Parameters
----------
TT_time : float
    Time in seconds since J2000, in the TT time scale.
Returns
-------
float
    Time in seconds since J2000, in the TDB time scale.






    )";



    } else if(name == "epoch_from_date_time_components" && variant==0) {
        return R"(
        
Computes the epoch as seconds since J2000 from the entries of the current date and time.

Computes the epoch as seconds since J2000. This function is added for convenience, and creates a :class:`DateTime` object, and subsequently calls its ``epoch`` function

Parameters
----------
year : int
    Calendar year

month : int
    Calendar month (value must be 1-12)     

day : int
    Calendar day in current month, value must be larger than 0, and smaller or equal to the number of days in the month       

hour : int
    Full hours into the current day (value must be 0-23)

minute : int
    Full minutes into the current hour (value must be 0-59)

seconds : float
    Number of seconds into the current minute. Note that this value is stored as ``long double`` in Tudat, which may be 64-bit or 80-bit (16 or 19 digits) depending on the compiler used.     

Returns
-------
float
    Time in seconds since J2000.






    )";



    } else if(name == "epoch_from_date_time_iso_string" && variant==0) {
        return R"(
        
Computes the epoch as seconds since J2000 from an ISO datetime string.

Computes the epoch as seconds since J2000. This function is added for convenience, and creates a :class:`DateTime` object, and subsequently calls its ``epoch`` function

Parameters
----------
iso_datetime : str
    Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)

Returns
-------
float
    Time in seconds since J2000.






    )";



    } else if(name == "date_time_from_epoch" && variant==0) {
        return R"(
        
Creates a Tudat-native :class:`DateTime` object from the seconds since J2000.


Parameters
----------
epoch : float
    Seconds since J2000

Returns
-------
DateTime
    Tudat ``DateTime`` object.






    )";



    } else if(name == "date_time_from_iso_string" && variant==0) {
        return R"(
        
Creates a Tudat-native :class:`DateTime` object from an ISO datetime string.


Parameters
----------
iso_datetime : str
    Date and time as ISO compatible string ("YYYY-MM-DDTHH:MM:SS.SSSSS..", where the T may be replaced with a space)

Returns
-------
DateTime
    Tudat ``DateTime`` object.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace two_body_dynamics {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "propagate_kepler_orbit" && variant==0) {
        return R"(
        
Function to propagate Keplerian elements to a later epoch, assuming an unperturbed system.

Function to propagate Keplerian elements to a later epoch, assuming an unperturbed system. This function will
take the initial Keplerian elements, and propagate the true anomaly in time as per the requested input. This
is done by converting true anomaly to mean anomaly, apply the constant rate in mean motion for the requested
time, and converting the result back to true anomaly. Currently both elliptic and hyperbolic orbits are supported. 
Parabolic orbits are not supported and will result in an error message.


Parameters
----------
initial_kepler_elements : numpy.ndarray
    Keplerian elements that are to be propagated (see :ref:`\`\`element_conversion\`\`` for order)
propagation_time : float
    Time for which the elements are to be propagated w.r.t. the initial elements
gravitational_parameter : float
    Gravitational parameter of central body used for propagation
root_finder : RootFinder, default = None
    Root finder used to solve Kepler's equation when converting mean to eccentric anomaly. When no root finder is specified, the default option of the mean to eccentric anomaly function is used (see :func:`~mean_to_eccentric_anomaly').
Returns
-------
numpy.ndarray
    Keplerian elements, propagated in time from initial elements assuming unperturbed dynamics. Note that the true anomaly is returned within the -PI to PI spectrum. If the user desires a different spectrum (possibly including the number of revolutions), these should be added by the user a posteriori.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace polyhedron_utilities {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "surface_area" && variant==0) {
        return R"(
        
Computes the surface area of a polyhedron [1]_.


Parameters
----------
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

Returns
-------
float
    Surface area.






    )";



    } else if(name == "volume" && variant==0) {
        return R"(
        
Computes the volume of a polyhedron [1]_.


Parameters
----------
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

Returns
-------
float
    Volume.






    )";



    } else if(name == "centroid" && variant==0) {
        return R"(
        
Computes the position of the centroid of a polyhedron [1]_.


Parameters
----------
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

Returns
-------
numpy.ndarray
    Position of the centroid.






    )";



    } else if(name == "modify_centroid" && variant==0) {
        return R"(
        
Modifies vertex coordinates of the polyhedron based on the desired position of the centroid.

Modifies the coordinates of the polyhedron vertices, such that the centroid of the modified polyhedron coincides
with the specified position. The centroid is computed according to Dobrovolskis [1]_.


Parameters
----------
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

desired_centroid : numpy.ndarray
    Desired position of the centroid.

Returns
-------
numpy.ndarray
    Vertices coordinates of the modified polyhedron, which has the specified centroid position.






    )";



    } else if(name == "inertia_tensor_from_density" && variant==0) {
        return R"(
        
Compute the inertia tensor of a polyhedron, from the density.

Computes the inertia tensor of a polyhedron, according to Dobrovolskis [1]_.

The mass distribution is defined using the density of the polyhedron. To instead use the gravitational
parameter see :func:`~tudatpy.astro.polyhedron_utilities.inertia_tensor_from_gravitational_parameter`.


Parameters
----------
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

density : float
    Density of the polyhedron

Returns
-------
numpy.ndarray
    Inertia tensor.






    )";



    } else if(name == "inertia_tensor_from_gravitational_parameter" && variant==0) {
        return R"(
        
Compute the inertia tensor of a polyhedron, from the gravitational parameter.

Computes the inertia tensor of a polyhedron, according to Dobrovolskis [1]_.

The mass distribution is defined using the gravitational parameter of the polyhedron. To instead use the density
see :func:`~tudatpy.astro.polyhedron_utilities.inertia_tensor_from_density`.


Parameters
----------
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

gravitational_parameter : float
    Gravitational parameter :math:`\mu` of gravity field.
Returns
-------
numpy.ndarray
    Inertia tensor.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace gravitation {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "legendre_normalization_factor" && variant==0) {
        return R"(
        
Function to calculate the normalization factor for spherical harmonics at a given degree and order

Function to calculate the normalization factor for spherical harmonics at a given degree and order.
Specifically, this function returns the value :math:`\mathcal{N}_{lm}`, computed from:

.. math::
    \mathcal{N}_{lm}=\sqrt{\frac{(2-\delta_{0m})(2l+1)(l-m)!)}{(l+m)!}}

with :math:`\delta_{0m}` is the Kronecker Delta function. The following can be used such that the conversion between unnormalized and fully
normalized spherical harmonic coefficients and Legendre polynomials can be computed from:

.. math::
    [C,S]_{lm}&=\mathcal{N}_{lm}[\bar{C},\bar{S}]_{lm}\\
    \bar{P}_{lm}&=\mathcal{N}_{lm}{P}_{lm}

with :math:`[C,S]_{lm}` denoting the unnormalized cosine or sine spherical harmonic coefficients at degree :math:`l` and order :math:`m`,
and :math:`P_{lm}` and :math:`\bar{P}_{lm}` representing the unnormalized and normalized associated Legendre polynomials at degree :math:`l` and order :math:`m`.


Parameters
----------
degree : int
    Spherical harmonic degree :math:`l`
order : int
    Spherical harmonic order :math:`m`
Returns
-------
float
    Normalization factor :math:`\mathcal{N}_{lm}`






    )";



    } else if(name == "normalize_spherical_harmonic_coefficients" && variant==0) {
        return R"(
        
Function to normalize spherical harmonic coefficients

Function to normalize spherical harmonic coefficients, using the equations provided in the :func:`~tudatpy.gravitation.astro.legendre_normalization_factor` function.

Parameters
----------
unnormalized_cosine_coefficients : numpy.ndarray
    Matrix for which entry :math:`(i,j)` contains the unnormalized cosine coefficient :math:`C_{lm}`
unnormalized_sine_coefficients : numpy.ndarray
    Matrix for which entry :math:`(i,j)` contains the unnormalized sine coefficient :math:`S_{lm}`
Returns
-------
tuple[numpy.ndarray, numpy.ndarray]
    Tuple of two matrices, containing the normalized coefficients :math:`\bar{C}_{lm}` (first) and :math:`\bar{S}_{lm}` (second)






    )";



    } else if(name == "unnormalize_spherical_harmonic_coefficients" && variant==0) {
        return R"(
        
Function to unnormalize spherical harmonic coefficients

Function to unnormalize spherical harmonic coefficients, using the equations provided in the :func:`~tudatpy.gravitation.astro.legendre_normalization_factor` function.

Parameters
----------
normalized_cosine_coefficients : numpy.ndarray
    Matrix for which entry :math:`(i,j)` contains the normalized cosine coefficient :math:`\bar{C}_{lm}`
normalized_sine_coefficients : numpy.ndarray
    Matrix for which entry :math:`(i,j)` contains the normalized sine coefficient :math:`\bar{S}_{lm}`
Returns
-------
tuple[numpy.ndarray, numpy.ndarray]
    Tuple of two matrices, containing the unnormalized coefficients :math:`{C}_{lm}` (first) and :math:`{S}_{lm}` (second)






    )";



    } else if(name == "spherical_harmonic_coefficients_from_inertia" && variant==0) {
        return R"(
        
Function to compute degree-two spherical harmonic coefficients from an inertia tensor

Function to compute degree-two spherical harmonic coefficients :math:`C_{20}`, :math:`C_{21}`, :math:`C_{22}`, :math:`S_{21}`, :math:`S_{22}` and from an inertia tensor :math:`\mathbf{I}`, according to the following relation"

.. math::
    \mathbf{I}=M R^2\left(\left(\begin{array}{ccc} \frac{C_{20}}{3}-2 C_{22} & -2 S_{22} & -C_{21} \\ -2 S_{22} & \frac{C_{20}}{3}+2 C_{22} & -S_{21} \\ -C_{21} & -S_{21} & -\frac{2 C_{20}}{3} \end{array}\right)+\bar{I} \mathbf{1}_{3 \times 3}\right)

with :math:`M` the mass of the body, and :math:`R` the reference radius of the spherical harmonic coefficients. The term :math:`\bar{I}` is the mean moment of inertia. The spherical harmonic
coefficients in the above matrix are unnormalized.


Parameters
----------
inertia tensor : numpy.ndarray[numpy.float64[3, 3]]
    Full inertia tensor :math:`\mathbf{I}` of the body for which spherical harmonic coefficients are to be computed.
gravitational_parameter : float
    Gravitational parameter :math:`\mu` of the body for which spherical harmonic coefficients are to be computed.
reference_radius : float
    Reference radius w.r.t. which spherical harmonic coefficients are to be computed.
output_normalized_coefficients : bool, default=True
    Boolean denoting whether the coefficients computed are normalized (if true) or unnormalized (if false)
Returns
-------
tuple[numpy.ndarray, numpy.ndarray]
    Tuple of two matrices, containing the spherical harmonic coefficients :math:`{C}_{lm}` (first) and :math:`{S}_{lm}` (second) up to degree and order 2.
    The degree-two coefficients are computed from the inertia tensor, the degree-one coefficients are set to zero (and :math:`C_{00}=0`)







    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace interface {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else {
        return "No documentation found.";
    }

}


    
namespace spice {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "SpiceEphemeris") {
         return R"(

        Ephemeris derived class which retrieves the state of a body directly from the SPICE library.

        The body of which the ephemeris is to be retrieved, as well as the origin and orientation
        of the reference frame in which the states are returned, and any corrections that are
        applied, are defined once during object construction.





     )";




    } else if(name == "SpiceEphemeris.ctor" && variant==0) {
            return R"(

        Constructor.

        Sets the input variables for the calls to the spice function to retrieve state.

        Parameters
        ----------
        target_body_name : str
            Name of body of which the ephemeris is to be calculated.
        observer_body_name : str
            Name of body relative to which the ephemeris is to be calculated.
        correct_for_stellar_aberration : bool, default = False
            Boolean whether to correct for stellar Aberration in
            retrieved values of (observed state).

        correct_for_light_time_aberration : bool, default = True
            Boolean whether to correct for light time in
            retrieved values of (observed state).

        converge_light_time_aberration : bool, default = False
            Boolean whether to use single iteration or max. 3
            iterations for calculating light time.

        reference_frame_name : str, default = "ECLIPJ2000"
            Name of the reference frame in which the ephemeris is to be
            calculated.

        reference_julian_day : float, default = constants.JULIAN_DAY_ON_J2000
            Reference julian day w.r.t. which ephemeris is evaluated.





    )";



    } else if(name == "SpiceEphemeris.get_cartesian_state" && variant==0) {
            return R"(

        Get Cartesian state from ephemeris.

         Returns Cartesian state from ephemeris at given Julian day.

        Parameters
        ----------
        seconds_since_epoch : float
            Seconds since epoch at which ephemeris is to be evaluated.




    )";





    } else if(name == "convert_julian_date_to_ephemeris_time" && variant==0) {
        return R"(
        
Convert a Julian date to ephemeris time (equivalent to TDB in Spice).

Function to convert a Julian date to ephemeris time, which is
equivalent to barycentric dynamical time. A leap second kernel
must have been loaded to use this function.


Parameters
----------
julian_date : int
    Julian date that is to be converted to ephemeris time.
Returns
-------
ephemeris_time : float    Julian date calculated from ephemeris time.






    )";



    } else if(name == "convert_ephemeris_time_to_julian_date" && variant==0) {
        return R"(
        
Convert ephemeris time (equivalent to TDB) to a Julian date.

Function to convert ephemeris time, which is nearly equal to
barycentric dynamical time, to the Julian date. A leap second
kernel must have been loaded to use this function.


Parameters
----------
ephemeris_time : float
    Ephemeris time that is to be converted to Julian date.
Returns
-------
julian_date : float    Julian date calculated from ephemeris time.






    )";



    } else if(name == "convert_date_string_to_ephemeris_time" && variant==0) {
        return R"(
        
Converts a date string to ephemeris time.

Function to convert a date string, for instance
1988 June 13, 3:29:48 to ephemeris time, wrapper for `str2et_c`
spice function.


Parameters
----------
date_string : str
    String representing the date. See documentation of spice
    function `str2et_c` for details on supported formats.

Returns
-------
ephemeris_time : str    Ephemeris time corresponding to given date_string.






    )";



    } else if(name == "get_body_cartesian_state_at_epoch" && variant==0) {
        return R"(
        
Get Cartesian state of a body, as observed from another body.

This function returns the state of a body, relative to another
body, in a frame specified by the user. Corrections for light-time
correction and stellar aberration can be applied to obtain the
state of one of the bodies, as observed from the other. Wrapper
for `spkezr_c` spice function.


Parameters
----------
target_body_name : str
    Name of the body of which the state is to be obtained. A kernel
    with the ephemeris of this body must have been loaded. The
    string must be a spice-recognized name or ID.

observer_body_name : str
    Name of the body relative to which the state is to be obtained.
    A kernel with the ephemeris of this body must have been loaded.
    The string must be a spice-recognized name or ID.

reference_frame_name : str
    The spice-recognized name of the reference frame in which the
    state is to be returned. Spice kernel(s) required to perform
    the necessary conversion from the states of the target and
    observer bodies to this frame need to have been loaded.

aberration_corrections : str
    Setting for correction for setting corrections. See Spice
    documentation for extended discussion.
    Short summary:

    - NONE: none
    - LT: light time corrected (one iteration for calculation)
    - CN: light time corrected (multiple iterations, max 3) for calculation
    - S: Stellar aberration corrected.
    - XLT and XCN: can be provided to make the ephemeris time input argument the transmission time, instead of reception time. Arguments can be combined (i.e."LT+S" or "XCN+S").

ephemeris_time : float
    Observation time (or transmission time of observed light, see description
    of aberrationCorrections).

Returns
-------
cartesian_state_vector : np.ndarray[6,]    Cartesian state vector (x,y,z, position+velocity).






    )";



    } else if(name == "get_body_cartesian_position_at_epoch" && variant==0) {
        return R"(
        
Get Cartesian position of a body, as observed from another body.

This function returns the position of a body, relative to another
body, in a frame specified by the user. Corrections for light-time
correction and stellar aberration can be applied to obtain the
state of one of the bodies, as observed from the other. Wrapper
for `spkpos_c` spice function.


Parameters
----------
target_body_name : str
    Name of the body of which the state is to be obtained. A kernel
    with the ephemeris of this body must have been loaded. The
    string must be a spice-recognized name or ID.

observer_body_name : str
    Name of the body relative to which the state is to be obtained.
    A kernel with the ephemeris of this body must have been loaded.
    The string must be a spice-recognized name or ID.

reference_frame_name : str
    The spice-recognized name of the reference frame in which the
    state is to be returned. Spice kernel(s) required to perform
    the necessary conversion from the states of the target and
    observer bodies to this frame need to have been loaded.

aberration_corrections : str
    Setting for correction for setting corrections. See Spice
    documentation for extended discussion.
    Short summary:

    - NONE: none
    - LT: light time corrected (one iteration for calculation)
    - CN: light time corrected (multiple iterations, max 3) for calculation,
    - S: Stellar aberration corrected.
    - XLT and XCN: can be provided to make the ephemeris time input argument the transmission time, instead of reception time. Arguments can be combined (i.e."LT+S" or "XCN+S").

ephemeris_time : float
    Observation time (or transmission time of observed light, see description
    of aberrationCorrections).






    )";



    } else if(name == "get_cartesian_state_from_tle_at_epoch" && variant==0) {
        return R"(
        
Get Cartesian state of a satellite from its two-line element set at a specified epoch.

This function retrieves the state of a satellite at a certain epoch
by propagating the SGP or SDP models (near-Earth resp. deep space)
with the given two-line elements (TLE). This function serves as a
wrapper for the `ev2lin_` function in CSpice.


Parameters
----------
epoch : float
    Time in seconds since J2000 at which the state is to be retrieved.
tle : :class:`~tudatpy.kernel.astro.ephemerides.Tle`
    Shared pointer to a Tle object containing the SGP/SDP model parameters as derived from the element set.
Returns
-------
cartesian_state_vector : np.ndarray[6,]    Cartesian state vector (x,y,z, position+velocity).






    )";



    } else if(name == "compute_rotation_matrix_derivative_between_frames" && variant==0) {
        return R"(
        
Computes time derivative of rotation matrix between two frames.

This function computes the derivative of the rotation matrix
between two frames at a given time instant. kernels defining the
two frames, as well as any required intermediate frames, at the
requested time must have been loaded. Wrapper for (part of) `sxform_c` spice function.


Parameters
----------
original_frame
    Reference frame from which the rotation is made.
new_frame
    Reference frame to which the rotation is made.
ephemeris_time
    Value of ephemeris time at which rotation is to be determined.
Returns
-------
Time derivative of rotation matrix from original to new frame at given time.






    )";



    } else if(name == "get_angular_velocity_vector_of_frame_in_original_frame" && variant==0) {
        return R"(
        
Computes the angular velocity of one frame w.r.t. to another frame.

Computes the angular velocity of one frame w.r.t. to another frame.
at a given time instant. kernels defining the two frames, as well
as any required intermediate frames, at the requested time must
have been loaded. Wrapper for `xf2rav_c`_ spice function (utilizing `sxform_c`_).

.. _`xf2rav_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/xf2rav_c.html
.. _`sxform_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sxform_c.html


Parameters
----------
original_frame
    Reference frame from which the rotation is made.
new_frame
    Reference frame to which the rotation is made.
ephemeris_time
    Value of ephemeris time at which rotation is to be determined.
Returns
-------
Angular velocity of newFrame w.r.t. originalFrame, expressed in originalFrame.






    )";



    } else if(name == "get_body_properties" && variant==0) {
        return R"(
        
Get property of a body from Spice.

Function to retrieve a property of a body from Spice, wraps the bodvrd_c Spice function.


Parameters
----------
body_name
    Name of the body of which the property is to be retrieved.
property
    Name of the property that is to be retrieved. Naming conventions can be found
    in the `bodvrd_c`_ function documentation.

    .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html

maximum_number_of_values : int
    Number of values by which the property is expressed (i.e. 1 for
    gravitational parameter, 3 for tri-axial ellipsoid principal axes).

Returns
-------
Property value(s) expressed in an STL vector of doubles.



Notes
-----
Function returns values with distance unit km, not m!




    )";



    } else if(name == "get_body_gravitational_parameter" && variant==0) {
        return R"(
        
Get gravitational parameter of a body.

This function retrieves the gravitational parameter of a body.
Wraps the `bodvrd_c`_ spice function with "GM" as property type.

.. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html


Parameters
----------
body
    Name of the body of which the parameter is to be retrieved.
Returns
-------
Gravitational parameter of requested body.






    )";



    } else if(name == "get_average_radius" && variant==0) {
        return R"(
        
Get the (arithmetic) mean of the three principal axes of the tri-axial ellipsoid shape.

Returns the (arithmetic) mean of the three principal axes of the
tri-axial ellipsoid shape of the requested body. Uses the `bodvrd_c` spice function with "RADII" as property type.

.. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html


Parameters
----------
body
    Name of the body of which the average radius is to be retrieved.
Returns
-------
Arithmetic mean of principal axes of tri-axial ellipsoid shape model of body.






    )";



    } else if(name == "convert_body_name_to_naif_id" && variant==0) {
        return R"(
        
Convert a body name to its NAIF identification number.

This function converts a body name to its NAIF identification
number. The NAIF id number is required for a number of spice
functions, whereas the name is easily interpretable by the user.
Wrapper for the ``bods2c_c`` function.


Parameters
----------
body_name
    Name of the body for which NAIF id is to be retrieved.
Returns
-------
NAIF id number for the body with bodyName.






    )";



    } else if(name == "check_body_property_in_kernel_pool" && variant==0) {
        return R"(
        
Check if a certain property of a body is in the kernel pool.

This function checks if a certain property of a body is in the
kernel pool. These properties are defined in PCK kernels. Their
names are given in the kernel file, typical names can be found in
the Spice documentation. Wrapper for the `bodfnd_c`_ function.

.. _`bodfnd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodfnd_c.html


Parameters
----------
body_name
    Name of the body of which the property is to be checked.
body_property
    Name of the property of which the presence is to be checked, not case-sensitive.
Returns
-------
bool
    True if property is in pool, false if not.






    )";



    } else if(name == "load_standard_kernels" && variant==0) {
        return R"(
        
Loads the default spice kernels shopped with tudat.

Loads the default spice kernels shopped with tudat. The kernels that are loaded are (in order):

- pck00010.tpc - Orientation and size/shape data for natural bodies, based mainly on  IAU Working Group on Cartographic Coordinates and Rotational Elements, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/>`_
- inpop19a_TDB_m100_p100_spice.tpc - Masses of solar system planets and large asteroids, as determined in INPOP19a ephemerides, obtained from `here <https://www.imcce.fr/recherche/equipes/asd/inpop/download19a>`_
- NOE-4-2020.tpc - Mars and Martian moon masses; Mars rotation model, as determined/used in NOE Martian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/>`_
- NOE-5-2021.tpc - Jupiter and selected Jovian moon (Io, Europa, Ganymede, Callisto, Amalthea) masses; Jupiter rotation model, as determined/used in NOE Jovian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/>`_
- NOE-6-2018-MAIN-v2.tpc - Saturn and selected Saturnian moon (Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus) masses; Saturn rotation model, as determined/used in NOE Saturnian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/>`_
- codes_300ast_20100725.bsp - Ephemerides of 300 of of the largest asteroids, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/>`_
- inpop19a_TDB_m100_p100_spice.bsp - Ephemerides of Solar system planetary system barycenters, Sun, Moon, Earth and Pluto, as determined in INPOP19a ephemerides, obtained from `here <https://www.imcce.fr/recherche/equipes/asd/inpop/download19a>`_
- NOE-4-2020.bsp - Mars, Phobos and Deimos ephemerides (w.r.t. Martian system barycenter), as determined/used in NOE Martian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/>`_
- NOE-5-2021.bsp - Jupiter, Io, Europa, Ganymede, Callisto, Amalthea ephemerides (w.r.t. Jovian system barycenter), as determined/used in NOE Jovian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/>`_
- NOE-6-2018-MAIN-v2.bsp - Saturn, Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus ephemerides (w.r.t. Saturnian system barycenter), as determined/used in NOE Saturnian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/>`_
- naif0012.tls - Leap second kernel, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/>`_


Parameters
----------
kernel_paths : list[str], default = None
    Optional alternative ephemeris kernels to be loaded, instead of the default ephemeris kernels. Note that using this input automatically prevents all of the above .bsp (ephemeris) kernels from loading.





    )";



    } else if(name == "get_total_count_of_kernels_loaded" && variant==0) {
        return R"(
        
Get the number of spice kernels currently loaded.

This function returns the amount of Spice kernels that are loaded
into the kernel pool. The same kernel can be loaded multiple times.
Wrapper for the `ktotal_c`_ function.

.. _`ktotal_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ktotal_c.html

Returns
-------
n_kernels : int    Number of spice kernels currently loaded.






    )";



    } else if(name == "load_kernel" && variant==0) {
        return R"(
        
Loads a Spice kernel into the pool.

This function loads a Spice kernel into the kernel pool, from which
it can be used by the various internal spice routines. Matters
regarding the manner in which Spice handles different kernels
containing the same information can be found in the spice required
reading documentation, kernel section. Wrapper for the `furnsh_c`_
function.

.. _`furnsh_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html


Parameters
----------
file_path : str
    Path to the spice kernel to be loaded.





    )";



    } else if(name == "clear_kernels" && variant==0) {
        return R"(
        
Clear all loaded spice kernels.

This function removes all Spice kernels from the kernel pool.
Wrapper for the `kclear_c`_ function.

.. _`kclear_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/kclear_c.html

Returns
-------
None
    None






    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace math {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else {
        return "No documentation found.";
    }

}


    
namespace interpolators {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "BoundaryInterpolationType") {
         return R"(

        Enumeration of types of behaviour to be used beyond the edges of the interpolation domain.

        Enumeration of types of behaviour to be used beyond the edges of the interpolation domain. For independent variable
        data in the range :math:`[t_{0}...t_{N}]`, this enum is used to define the behaviour of the interpolator at
        :math:`t<t_{0}` and :math:`t>t_{N}`





     )";


    } else if(name == "BoundaryInterpolationType.throw_exception_at_boundary") {
         return R"(
The program will terminate with an error message when the interpolator is interrogated beyond the range :math:`[t_{0}...t_{N}]`
     )";


    } else if(name == "BoundaryInterpolationType.use_boundary_value") {
         return R"(
The value :math:`\mathbf{x}_{0}` is returned for :math:`t<t_{0}` (and :math:`\mathbf{x}_{N}` if :math:`t>t_{N}`)
     )";


    } else if(name == "BoundaryInterpolationType.use_boundary_value_with_warning") {
         return R"(
Same as ``use_boundary_value``, but a warning is printed to the terminal
     )";


    } else if(name == "BoundaryInterpolationType.extrapolate_at_boundary") {
         return R"(
The interpolation scheme is extended beyond the range :math:`t_{0}...t_{N}` without any warning. That is, the mathematical equation used to compute the value of :math:`x` in the range :math:`[t_{0}...t_{1}]` is used without any checks for :math:`t<t_{0}`  (and equivalently for :math:`t>t_{N}`). Warning, using this setting can result in divergent/unrealistic behaviour
     )";


    } else if(name == "BoundaryInterpolationType.extrapolate_at_boundary_with_warning") {
         return R"(
Same as ``extrapolate_at_boundary``, but a warning is printed to the terminal
     )";


    } else if(name == "BoundaryInterpolationType.use_nan_value") {
         return R"(
The program will return an interpolated value filled with NaN entries.
     )";


    } else if(name == "BoundaryInterpolationType.use_nan_value_with_warning") {
         return R"(
Same as ``use_nan_value``, but a warning is printed to the terminal
     )";



    } else if(name == "AvailableLookupScheme") {
         return R"(

        Enumeration of types of behaviour to be used beyond the edges of the interpolation domain.

        When the interpolation is performed, the interpolator scheme will typically start by finding the nearest neighbor of
        the requested value of the independent variable :math:`t` in the data set :math:`[t_{0}...t_{N}]`.
        The choice of lookup scheme can have a significant influence on computational efficiency for large data sets and/or simple
        interpolation algorithms





     )";


    } else if(name == "AvailableLookupScheme.hunting_algorithm") {
         return R"(
With this option, the interpolator 'remembers' which value of :math:`t_{i}` was the nearest neighbor during the previous call to the interpolate function, and starts looking at/near this entry of the data set :math:`[t_{i}]` to find the nearest neighbor.
     )";


    } else if(name == "AvailableLookupScheme.binary_search") {
         return R"(
With this option, the algorithm uses a binary search algorithm to find the nearest neighbor, initially starting with the full data range :math:`[t_{0}...t_{N}]`.
     )";



    } else if(name == "LagrangeInterpolatorBoundaryHandling") {
         return R"(

        Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.

        Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.
        As explained for :func:`lagrange_interpolation`, the algorithm for the Lagrange interpolation breaks down at the edges of
        the interpolation domain. This enum provides the available options a user has to deal with this.





     )";


    } else if(name == "LagrangeInterpolatorBoundaryHandling.lagrange_cubic_spline_boundary_interpolation") {
         return R"(
A cubic-spline interpolator is created from the first and last :math:`\max(m/2-1,4)` data points of the full data set, and these cubic spline interpolators are used when an interpolation at :math:`t<t_{(m/2-1)}` or :math:`t<t_{N-(m/2)}` is called
     )";


    } else if(name == "LagrangeInterpolatorBoundaryHandling.lagrange_cubic_spline_boundary_interpolation_with_warning") {
         return R"(
Same as ``lagrange_cubic_spline_boundary_interpolation``, but a warning is printed to the terminal
     )";


    } else if(name == "LagrangeInterpolatorBoundaryHandling.lagrange_no_boundary_interpolation") {
         return R"(
The program will terminate with an exception when the Lagrange interpolator is interrogated beyond its valid range
     )";


    } else if(name == "LagrangeInterpolatorBoundaryHandling.lagrange_boundary_nan_interpolation") {
         return R"(
The program will return an interpolated value filled with NaN entries.
     )";


    } else if(name == "LagrangeInterpolatorBoundaryHandling.lagrange_boundary_nan_interpolation_with_warning") {
         return R"(
Same as ``lagrange_boundary_nan_interpolation``, but a warning is printed to the terminal
     )";




    } else if(name == "InterpolatorSettings") {
         return R"(

        Base class to define settings for an interpolator.





     )";





    } else if(name == "LagrangeInterpolatorSettings") {
         return R"(

        :class:`InterpolatorSettings`-derived class to define settings for a Lagrange interpolator.





     )";





    } else if(name == "OneDimensionalInterpolatorScalar") {
         return R"(

        Object that performs interpolation for scalar independent, and scalar dependent variables.

        Object that performs interpolation for scalar independent, and scalar dependent variables. This object is
        not created manually, but is set up using the :func:`create_one_dimensional_scalar_interpolator` function.





     )";




    } else if(name == "OneDimensionalInterpolatorScalar.interpolate" && variant==0) {
            return R"(

        This function performs the interpolation at the requested independent variable value.


        Parameters
        ----------
        independent_variable_value : float
            Value of independent variable at which the interpolation is to bse performed.

        Returns
        -------
        float
            Interpolated dependent variable value, using implemented algorithm at requested independent variable value





    )";




    } else if(name == "OneDimensionalInterpolatorVector") {
         return R"(

        Object that performs interpolation for vector independent, and vector dependent variables.

        Object that performs interpolation for vector independent, and vector dependent variables. This object is
        not created manually, but is set up using the :func:`create_one_dimensional_vector_interpolator` function.





     )";




    } else if(name == "OneDimensionalInterpolatorVector.interpolate" && variant==0) {
            return R"(

        This function performs the interpolation at the requested independent variable value.


        Parameters
        ----------
        independent_variable_value : float
            Value of independent variable at which the interpolation is to be performed.

        Returns
        -------
        np.array
            Interpolated dependent variable value, using implemented algorithm at requested independent variable value





    )";




    } else if(name == "OneDimensionalInterpolatorMatrix") {
         return R"(

        Object that performs interpolation for matrix independent, and matrix dependent variables.

        Object that performs interpolation for matrix independent, and matrix dependent variables. This object is
        not created manually, but is set up using the :func:`create_one_dimensional_matrix_interpolator` function.





     )";




    } else if(name == "OneDimensionalInterpolatorMatrix.interpolate" && variant==0) {
            return R"(

        This function performs the interpolation at the requested independent variable value.


        Parameters
        ----------
        independent_variable_value : float
            Value of independent variable at which the interpolation is to be performed.

        Returns
        -------
        np.array
            Interpolated dependent variable value, using implemented algorithm at requested independent variable value





    )";





    } else if(name == "linear_interpolation" && variant==0) {
        return R"(
        
Function to create settings for linear interpolation.

Function to create settings for linear interpolation, where the interpolator
defines a linear curve between each two subsequent intervals of the
independent variable input data.


Parameters
----------
lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
    Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
boundary_interpolation : BoundaryInterpolationType, default=extrapolate_at_boundary_with_warning
    Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
Returns
-------
InterpolatorSettings
    Linear interpolation settings object






    )";



    } else if(name == "piecewise_constant_interpolation" && variant==0) {
        return R"(
        
Function to create settings for piecewise constant interpolation.

Function to create settings for piecewise constant interpolation. If interpolator
is to return the value at :math:`t`, and :math:`t_{i}\le t < t_{i+1}`, the interpolator
returns :math:`\mathbf{x}_{i}`


Parameters
----------
lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
    Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
    Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
Returns
-------
InterpolatorSettings
    Piecewise constant interpolation settings object






    )";



    } else if(name == "cubic_spline_interpolation" && variant==0) {
        return R"(
        
Function to create settings for cubic spline interpolation.

Function to create settings for cubic spline interpolation, where the interpolator
defines a cubic curve polynomial curve between each two subsequent intervals of the
independent variable input data. The curve has continuous value, first derivative and
second derivate between subsequent intervals. As boundary condition, the spline has
a zero second derivative imposed at the upper and lower boundaries of the interpolation
domain.


Parameters
----------
lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
    Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
    Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
Returns
-------
InterpolatorSettings
    Cubic spline settings object






    )";



    } else if(name == "hermite_spline_interpolation" && variant==0) {
        return R"(
        
Function to create settings for cubic Hermite spline interpolation.

Function to create settings for piecewise cubic Hermite spline interpolation. To use this
interpolator, a key-value container of values, and a key-value container of first derivatives,
must be provided to the function creating an interpolator (:func:`create_one_dimensional_scalar_interpolator`,
:func:`create_one_dimensional_vector_interpolator`,  :func:`create_one_dimensional_matrix_interpolator`). The resulting
spline uses the value and first derivatives (four piece of information for each interval) at two subsequent nodes to construct
a cubic polynomial between each two subsequent nodes. The resulting spline has constant values and first derivatives


Parameters
----------
lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
    Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
    Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
Returns
-------
InterpolatorSettings
    Hermite spline interpolation settings object






    )";



    } else if(name == "lagrange_interpolation" && variant==0) {
        return R"(
        
Function to create settings for cubic Lagrange interpolation.

Function to create settings for piecewise cubic Lagrange interpolation.  This is typically the interpolator of highest accuracy that is available.
The Lagrange interpolator uses :math:`m` consecutive points (input to this function) from the input independent variables :math:`[t_{0}...t_{N}]`
to create the polynomial of order :math:`m-1` that interpolates these points. From here on, we assume :math:`m` is even.
The algorithm that is used (see `here <https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html>`_ for mathematical details
on interpolating Lagrange polynomials) works as follows:

* The nearest lower neighbor of the data point :math:`t` at which the state :math:`\mathbf{x}` is to be interpolated is determined, and denoted :math:`t_{i}`.
* An interpolating Lagrange polynomial is constructed from the consecutive data points :math:`[t_{i-(m/2-1)}...t_{i+m}]`
* This resulting interpolating polynomial is *only* used in the interval :math:`[t_{i}...t_{i+1}]`, to prevent `Runge's phenomenon <https://en.wikipedia.org/wiki/Runge%27s_phenomenon>`_.

For instance, if :math:`m=8` we use a :math:`7^{th}` order polynomial that interpolates a contiguous set of
8 data points out of the full data set. Normally, the interpolating polynomial is only used between the
:math:`4^{th}` and :math:`5^{th}` data point, where it will typically be of good accuracy. Consequently,
a separate interpolating polynomial (using data over a span of :math:`m` consecutive points) is used for
each single interval :math:`[t_{i}...t_{i+1}]` (with the exception of the boundaries, see below).

.. warning:: Issues can occur if the data point :math:`t` at which the interpolation is to be
             performed is close to :math:`t_{0}` or :math:`t_{N}`. In those case, there is not sufficient
             data to construct the interpolating polynomial *and* to only use this interpolating polynomial
             between the middle two data points that were used to it. In these cases, the user has a number of
             options (all defined by an entry of the :class:`LagrangeInterpolatorBoundaryHandling` variable,
             used as input to this function). In short, interpolation between the first and last :math:`m/2`
             data points will lead to degraded results, warnings, or termination.


Parameters
----------
number_of_points : int
    Number of consecutive data points that are used to construct a single interpolating polynomial.
lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
    Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
    Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
lagrange_boundary_handling : LagrangeInterpolatorBoundaryHandling, default = lagrange_cubic_spline_boundary_interpolation
    Interpolator behaviour that is to be used at the boundaries of the domain, where the regular algorithm cannot be executed.
Returns
-------
LagrangeInterpolatorSettings
    Lagrange interpolation settings object






    )";



    } else if(name == "create_one_dimensional_scalar_interpolator" && variant==0) {
        return R"(
        
Function to create an interpolator for scalar dependent variables.

Function to create an interpolator for scalar dependent variables, with a single independent
variable. This function takes the interpolator settings, and the data that is to be interpolated,
as input to create the object that can perform the actual interpolation


Parameters
----------
data_to_interpolate : dict[float, float]
    Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
interpolator_settings : InterpolatorSettings
    Settings that define the type of interpolator that is to be used
data_first_derivatives : dict[float, float] = dict()
    Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator
Returns
-------
OneDimensionalInterpolatorScalar
    Interpolator object






    )";



    } else if(name == "create_one_dimensional_vector_interpolator" && variant==0) {
        return R"(
        
Function to create an interpolator for vector dependent variables.

As :func:`create_one_dimensional_scalar_interpolator`, but with vectors as dependent variables


Parameters
----------
data_to_interpolate : dict[float, np.array]
    Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
interpolator_settings : InterpolatorSettings
    Settings that define the type of interpolator that is to be used
data_first_derivatives : dict[float, np.array] = dict()
    Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
Returns
-------
OneDimensionalInterpolatorVector
    Interpolator object






    )";



    } else if(name == "create_one_dimensional_matrix_interpolator" && variant==0) {
        return R"(
        
Function to create an interpolator for matrix dependent variables.

As :func:`create_one_dimensional_scalar_interpolator`, but with matrices (2-dimensional arrays) as dependent variables


Parameters
----------
data_to_interpolate : dict[float, np.array]
    Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
interpolator_settings : InterpolatorSettings
    Settings that define the type of interpolator that is to be used
data_first_derivatives : dict[float, np.array] = dict()
    Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator
Returns
-------
OneDimensionalInterpolatorMatrix
    Interpolator object






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace root_finders {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "MaximumIterationHandling") {
         return R"(

        Enumeration of types of behaviour to be used when the convergence criterion on maximum number of iterations is reached.





     )";


    } else if(name == "MaximumIterationHandling.accept_result") {
         return R"(
The program will accept the root at the final iteration, without any additional output
     )";


    } else if(name == "MaximumIterationHandling.accept_result_with_warning") {
         return R"(
The program will accept the root at the final iteration, but will print a warning to the terminal that the root finder may not have converged
     )";


    } else if(name == "MaximumIterationHandling.throw_exception") {
         return R"(
The program will not accept the root at the final iteration, and will throw an exception
     )";




    } else if(name == "RootFinderSettings") {
         return R"(

        Class to define settings for a root finder.





     )";






    } else if(name == "bisection" && variant==0) {
        return R"(
        
Function to create settings for a bisection root-finder.

Function to create settings for a bisection root finder. This root finder approximates the root by initializing with
two initial guesses :math:`x_{\downarrow,0}` and :math:`x_{\uparrow,0}`, for which it is required that
:math:`f(x_{\downarrow,0}) < 0` and :math:`f(x_{\uparrow,0}) > 0`. At each iteration :math:`i`, the current guess of
the root :math:`x_{i}` is:

.. math::
   x_{i}=\begin{cases}
   x_{\downarrow,i}, & |f(x_{\downarrow,i})|<|f(x_{\uparrow,i})|\\
	  x_{\uparrow,i}, & |f(x_{\downarrow,i})|\ge|f(x_{\uparrow,i})|
	       \end{cases}

The midpoint :math:`x_{m,i}` of :math:`x_{\downarrow,i}` and :math:`x_{\uparrow,i}` is then computed from :math:`x_{m,i}=(x_{\downarrow,i}-x_{\uparrow,i})/2`.
Depending on the sign of :math:`f(x_{m,i})`, it then replaces either :math:`x_{\downarrow,i}` or :math:`x_{\uparrow,i}` (depending on whether
its sign matches :math:`f(x_{\downarrow,i})` for iteration :math:`i+1` and :math:`f(x_{\uparrow,i})`), while the other point from iteration :math:`i` is retained. 

Although slow, the algorithm is ensured to converge to a root, if the two initial guesses indeed have opposite signs (if not, an exception is thrown).


Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Bisection root-finding settings object






    )";



    } else if(name == "newton_raphson" && variant==0) {
        return R"(
        
Function to create settings for a Newton-Raphson root-finder.

Function to create settings for a bisection root finder. This root finder approximates the root by initializing with
a single initial guesses :math:`x_{0}` and requires an analytical formulation for :math:`f(x)` and :math:`f'(x)=\frac{d}{dx}f(x)`.
The algorithm uses the following equation to iterate:
      
.. math::
   x_{i+1}=x_{i}-\frac{f(x_{i})}{f'(x_{i})}


Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Newton-Raphson root-finding settings object






    )";



    } else if(name == "secant" && variant==0) {
        return R"(
        
Function to create settings for a secant method root-finder.

Function to create settings for a root finder using the secant method. This root finder approximates the root by initializing with
two initial guesses :math:`x_{0}` and :math:`x_{1}`. The algorithm uses the following equation to iterate:
      
.. math::
   x_{i+1}=x_{i}-f(x_{i})\frac{x_{i}-x_{i-1}}{f(x_{i})-f(x_{i-1})}
   


Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Secant root-finding settings object






    )";



    } else if(name == "halley" && variant==0) {
        return R"(
        
Function to create settings for a Halley root-finder.

Function to create settings for a Halley root finder. This root finder approximates the root by initializing with
a single initial guesses :math:`x_{0}` and requires an analytical formulation for :math:`f(x)`, :math:`f'(x)=\frac{d}{dx}f(x)` and :math:`f''(x)=\frac{d^{2}}{dx^{2}}f(x)`.
The algorithm uses the following equation to iterate:
      
.. math::
   x_{i+1}=x_{i}-\frac{2f(x_{i})f'(x_{i})}{2(f'(x_{i}))^{2}-f(x_{i})f''(x_{i})}
   
   


Parameters
----------
relative_variable_tolerance : float, default = nan
    Relative tolerance :math:`\epsilon_{r}` (setting not used if nan)
absolute_variable_tolerance : float, default = nan
    Relative absolute :math:`\epsilon_{a}` (setting not used if nan)
root_function_tolerance : float, default = nan
    Root function tolerance :math:`\epsilon_{f}` (setting not used if nan)
maximum_iteration : int, default = 1000
    Maximum number of iterations :math:`N`
maximum_iteration_handling : MaximumIterationHandling, default = throw_exception
    Algorithm behaviour if maximum number of iterations :math:`N` is reache
Returns
-------
RootFinderSettings
    Halley root-finding settings object






    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace numerical_simulation {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "SingleArcSimulator") {
         return R"(

        Class for consolidating single arc dynamics simulation functionality.

        Class for consolidating all functionality required to perform single arc dynamics simulations.





     )";


    } else if(name == "SingleArcSimulator.integrator_settings") {
         return R"(

        Settings to create the numerical integrator that is to be used
        for the integration of the equations of motion


        :type: IntegratorSettings
     )";


    } else if(name == "SingleArcSimulator.state_derivative_function") {
         return R"(

        **read-only**

        Function that performs a single state derivative function evaluation. This function takes the numerically propagated
        state, and current independent variable (time) as input, and returns the derivative of the state that is then used
        by the numerical integration routine. Typically, this function is NOT used directly by users.


        :type: Callable[[float, numpy.ndarray], numpy.ndarray]
     )";


    } else if(name == "SingleArcSimulator.environment_updater") {
         return R"(

        **read-only**

        Object used in the propagation to update the environment, it uses the current time and numerically calculated state
        to update the translational state, rotational state, flight conditions, etc. of all bodies in the simulation to be
        consistent with this time and state.  Typically, this class is NOT used directly by users, but can be useful in specific situations.


        :type: EnvironmentUpdater
     )";


    } else if(name == "SingleArcSimulator.propagation_results") {
         return R"(

        **read-only**

        This function retrieves all the results of the numerical propagation, stored
        in a single wrapper object


        :type: SingleArcSimulationResults
     )";




    } else if(name == "SingleArcSimulator.integrate_equations_of_motion" && variant==0) {
            return R"(

        This function numerically (re-)integrates the equations of
        motion.


        This function numerically (re-)integrates the equations of
        motion, using the settings set through the constructor and a
        new initial state vector provided here.


        Parameters
        ----------
        initial_states : numpy.ndarray
            Initial state vector that is to be used for numerical
            integration. Note that this state should be in the correct
            frame (i.e. relative to central_bodies in
            propagator_settings), and in Cartesian elements (even when
            using a different propagation scheme than Cowell).

            .. note:: When using default settings for the class
                      constructor, this function is called during object
                      creation.





    )";




    } else if(name == "SingleArcVariationalSimulator") {
         return R"(

        Class for consolidating single arc variational dynamics functionality.

        Class for consolidating all functionality required to perform single arc variational dynamics simulations.





     )";


    } else if(name == "SingleArcVariationalSimulator.parameter_vector") {
         return R"(

        Consolidated set of (estimatable) parameters
        w.r.t. the variational dynamics in the Variational Simulator are defined.


        :type: :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
     )";


    } else if(name == "SingleArcVariationalSimulator.variational_equations_history") {
         return R"(

        **read-only**

        List containing the solution of the variational equations, i.e. the
        state transition matrix history (first entry) and sensitivity matrix history (second vector entry).


        :type: list[ dict[float, numpy.ndarray] ]
     )";


    } else if(name == "SingleArcVariationalSimulator.state_transition_matrix_history") {
         return R"(

        **read-only**

        State transition matrix history, given as epoch with propagation epochs as keys.
        This is (alongside the `sensitivity_matrix_history`) the solution of the variational equations.


        :type: dict[float, numpy.ndarray]
     )";


    } else if(name == "SingleArcVariationalSimulator.sensitivity_matrix_history") {
         return R"(

        **read-only**

        Sensitivity matrix history, given as epoch with propagation epochs as keys.
        This is (alongside the `state_transition_matrix_history`) the solution of the variational equations.


        :type: dict[float, numpy.ndarray]
     )";


    } else if(name == "SingleArcVariationalSimulator.state_history") {
         return R"(

        **read-only**

        State history, given as epoch with propagation epochs as keys.
        This is the solution of the (propagated) equations of motion, describing the states along which
        the variational dynamics are solved.


        :type: dict[float, numpy.ndarray]
     )";


    } else if(name == "SingleArcVariationalSimulator.dynamics_simulator") {
         return R"(

        **read-only**

        Simulator object containing all functionality for solving of the (regular) equations of motion.


        :type: :class:`~tudatpy.numerical_simulation.SingleArcSimulator`
     )";




    } else if(name == "SingleArcVariationalSimulator.ctor" && variant==0) {
            return R"(

        Class constructor.

        Constructor through which the user can create instances of this class.
        Defines environment, propagation and integrations models, as well as a number of settings related
        to the (estimatable) parameters, w.r.t. which the variational equations are defined.

        .. note:: When using default settings, creating an object of
                  this type automatically triggers the propagation


        Parameters
        ----------
        bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
            Object defining the physical environment, with all
            properties of artificial and natural bodies.

        integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`
            Settings to create the numerical integrator that is to be used for the integration of the equations of motion.

        propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
            Settings to create the propagator that is to be used for the propagation of the dynamics.

        estimated_parameters : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
            Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
            linked to the environment and acceleration settings of the simulation.

        integrate_equations_concurrently : Bool, default = True
            Boolean defining whether equations of motion and variational equations are to be propagated concurrently
            (if true) or sequentially (of false).

        variational_only_integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`, default = []
            Settings to create the numerical integrator that is to be used for integration the variational equations.
            If none is given (default), the numerical integration settings are taken to be the same as the ones applied
            in the integration of the equations of motions (specified by the `integrator_settings` parameter).

        clear_numerical_solutions : Bool, default = False
            Boolean to determine whether to clear the raw numerical solution member variables
            and to reset the state transition interface after propagation.

        integrate_on_creation : Bool, default = True
            Boolean defining whether the propagation should be performed immediately (default), or at a later time
            (when calling the :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` member function).

        set_integrated_result : Bool, default = True
            Boolean to determine whether to automatically use the integrated results to set ephemerides for the
            propagated bodies.





    )";



    } else if(name == "SingleArcVariationalSimulator.integrate_equations_of_motion_only" && variant==0) {
            return R"(

        Function to trigger the integration of the (regular) equations of motion.


        Function to trigger the integration only of the (regular) equations of motion, resulting in a `state_history`.
        This step does not yet use variational dynamics. In order to also solve the variational equations,
        use the :func:`integrate_full_equations` member function.

        Returns
        -------
        None
            Creates / modifies the `state_history` property of the :class:`~tudatpy.numerical_simulation.SingleArcVariationalSolver` object.





    )";



    } else if(name == "SingleArcVariationalSimulator.integrate_full_equations" && variant==0) {
            return R"(

        Function to trigger the integration of variational and dynamical equations (equations of motion).


        Function to trigger the integration of the (regular) equations of motion as well as the variational equations,
        solving for `state_history` and `variational_equations_history`
        (in its two components `state_transition_matrix_history` & `sensitivity_matrix_history`).


        Parameters
        ----------
        initial_states : numpy.ndarray([m, 1])
            Initial state to be used for the parameters in the equations of motion.

        integrate_equations_concurrently : Bool, default = True
            Boolean defining whether equations of motion and variational equations are to be propagated concurrently
            (if true) or sequentially (of false).

        Returns
        -------
        None
            Creates / modifies the properties of the VariationalSolver object.





    )";




    } else if(name == "Estimator") {
         return R"(

        Class for consolidating all estimation functionality.

        Class for consolidating all functionality required to perform an estimation.





     )";


    } else if(name == "Estimator.observation_simulators") {
         return R"(

        **read-only**

        Observation simulators contained in the Estimator object. A single observation simulator hosts
        the functionality for simulating a given observable over the defined link geometry.


        :type: list[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
     )";


    } else if(name == "Estimator.observation_managers") {
         return R"(

        **read-only**

        Observation managers contained in the Estimator object. A single observation manager can simulate observations and
        calculate observation partials for all link ends involved in the given observable type.


        :type: dict[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`, :class:`~tudatpy.numerical_simulation.estimation.ObservationManager` ]
     )";


    } else if(name == "Estimator.state_transition_interface") {
         return R"(

        **read-only**

        State transition and sensitivity matrix interface, setting the variational equations/dynamics in the
        Estimator object.


        :type: :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
     )";


    } else if(name == "Estimator.variational_solver") {
         return R"(

        **read-only**

        Variational equations solver, which is used to manage and execute the numerical integration of
        equations of motion and variational equations/dynamics in the Estimator object.


        :type: :class:`~tudatpy.numerical_simulation.SingleArcVariationalSolver`
     )";




    } else if(name == "Estimator.ctor" && variant==0) {
            return R"(

        Class constructor.

        Constructor through which the user can create instances of this class.
        Defines environment, propagation and integrations models, as well as a number of settings related
        to the estimatable parameters and observation settings.

        .. note:: When using default settings, creating an object of
                  this type automatically triggers the propagation


        Parameters
        ----------
        bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
            Object defining the physical environment, with all
            properties of artificial and natural bodies.

        estimated_parameters : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
            Object defining a consolidated set of estimatable parameters,
            linked to the environment and acceleration settings of the simulation.

        observation_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
            List of settings objects, each object defining the observation model settings for one
            combination of observable and link geometry that is to be simulated.

        integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`
            Settings to create the numerical integrator that is to be
            used for the integration of the equations of motion

        propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
            Settings to create the propagator that is to be
            used for the propagation of dynamics

        integrate_on_creation : Bool, default = True
            Boolean defining whether the propagation should be
            performed immediately (default), or at a later time
            (when calling the :func:`perform_estimation` member function.





    )";



    } else if(name == "Estimator.compute_covariance" && variant==0) {
            return R"(

        Function to perform a covariance analysis for the given observations and parameters


        Function to perform a covariance analysis for the given observations and parameters. The observations are provided through the 
        ``covariance_analysis_input`` input, as are the weights :math:`\mathbf{W}` and inverse a priori covariance :math:`(\mathbf{P}_{0})^{-1}`.
        Calling this function uses the environment and propagator settings provided to the constructor of this `Estimator` class to simulate
        the dynamics of any relevant bodies for the observations (and associated variational equations). The observations are then
        computed using the observation models created by the settings provided to the constructor of this `Estimator` class, as is the
        associated design matrix :math:`\mathbf{H}`. This function then produces the covariance :math:`\mathbf{P}` (omitting the normalization used
        internally for numerical stability)

        .. math::
           \mathbf{P}=\left(\mathbf{H}^{T}\mathbf{W}\mathbf{H}+(\mathbf{P}_{0})^{-1}\right)^{-1}

        Note that, although the actual observations are formally not required for a covariance analysis, all additional data (e.g. observation time, type, link ends, etc.)
        are. And, as such, the ``covariance_analysis_input`` does require the full set of observations and associated information, for consistency purposes (e.g., same input as
        ``perform_estimation`` function) .


        Parameters
        ----------
        covariance_analysis_input : :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput`
            Object consolidating all relevant settings for the covariance analysis
            This includes foremost the simulated observations, as well as a priori information about the estimatable parameters

        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.vOutput`
            Object containing all outputs from the estimation process.





    )";



    } else if(name == "Estimator.perform_estimation" && variant==0) {
            return R"(

        Function to trigger the parameter estimation.


        Function to trigger the parameter estimation. Much of the process and requirements are similar to those described in the 
        :func:`~tudatpy.numerical_simulation.Estimator.compute_covariance` function. This function uses an iterative least-squares
        estimate process to fit the data (inside ``estimation_input``) to the model defined by the inputs to the ``Estimator`` constructor.s


        Parameters
        ----------
        estimation_input : :class:`~tudatpy.numerical_simulation.estimation.EstimationInput`
            Object consolidating all relevant settings for the estimation
            This includes foremost the simulated observations, as well as a priori information about the estimatable parameters and convergence criteria for the least squares estimation.

        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.EstimationOutput`
            Object containing all outputs from the estimation process.





    )";





    } else if(name == "create_dynamics_simulator" && variant==0) {
        return R"(
        
Function to create object that propagates the dynamics.

Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
Depending on the specific input type (e.g. which factory function from the :ref:`\`\`propagator\`\`` module was used),
a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
function. When using default settings, calling this function will automatically propagate the dynamics.


Parameters
----------
bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object defining the physical environment, with all
    properties of artificial and natural bodies.

propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
    Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

simulate_dynamics_on_creation : Bool, default=True
    Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
    can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function

Returns
-------
:class:`~tudatpy.numerical_simulation.Simulator`
    Object that propagates the dynamics, and processes the results.






    )";



    } else {
        return "No documentation found.";
    }

}


    
namespace estimation_setup {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "print_parameter_names" && variant==0) {
        return R"(
        
Function for printing a list of estimatable parameter names.

Function that allows you to print a verbose list of all parameters that shall be estimated. Consider parameters are listed separately.


Parameters
----------
parameter_set : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet` Instance of :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.
    None





    )";



    } else if(name == "create_parameter_set" && variant==0) {
        return R"(
        
Function for creating a consolidated set of estimatable parameters.

Function for creating a consolidated parameter from the given estimatable parameter settings.
The function checks for consistency between the parameter settings and the models contained in the simulation setup (given by the `bodies` & and `propagator_settings` parameters).


Parameters
----------
parameter_settings : list( :class:`~tudatpy.numerical_simulation.estimation_setup.EstimatableParameterSettings` )
    List of objects that define the settings for the parameters that are to be created. Each entry in this list is typically created by a call to a factory function in the :ref:`\`\`parameter\`\`` module

bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
    Object containing the consolidated propagation settings of the simulation.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
    Instance of :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.







    )";



    } else if(name == "create_observation_simulators" && variant==0) {
        return R"(
        
Function for creating observation simulator objects.

Factory function for creating observation simulator objects from observation settings.
Note that each observation (i.e. combination of observable and link geometry) requires its own observation simulator object.


Parameters
----------
observation_settings : List[ ObservationSettings ]
    List of settings objects, each object defining the observation model settings for one combination of observable and link geometry that is to be simulated.

bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

Returns
-------
List[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
    List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one combination of observable type and link geometry.






    )";



    } else {
        return "No documentation found.";
    }

}


    
namespace observation {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "LinkEndType") {
         return R"(

        Enumeration of available link end types.





     )";


    } else if(name == "LinkEndType.unidentified_link_end") {
         return R"(
     )";


    } else if(name == "LinkEndType.transmitter") {
         return R"(
     )";


    } else if(name == "LinkEndType.reflector1") {
         return R"(
     )";


    } else if(name == "LinkEndType.retransmitter") {
         return R"(
     )";


    } else if(name == "LinkEndType.reflector2") {
         return R"(
     )";


    } else if(name == "LinkEndType.reflector3") {
         return R"(
     )";


    } else if(name == "LinkEndType.reflector4") {
         return R"(
     )";


    } else if(name == "LinkEndType.receiver") {
         return R"(
     )";


    } else if(name == "LinkEndType.observed_body") {
         return R"(
     )";



    } else if(name == "ObservableType") {
         return R"(

        Enumeration of available observable types.





     )";


    } else if(name == "ObservableType.one_way_range_type") {
         return R"(
     )";


    } else if(name == "ObservableType.n_way_range_type") {
         return R"(
     )";


    } else if(name == "ObservableType.angular_position_type") {
         return R"(
     )";


    } else if(name == "ObservableType.relative_angular_position_type") {
         return R"(
     )";


    } else if(name == "ObservableType.position_observable_type") {
         return R"(
     )";


    } else if(name == "ObservableType.velocity_observable_type") {
         return R"(
     )";


    } else if(name == "ObservableType.one_way_instantaneous_doppler_type") {
         return R"(
     )";


    } else if(name == "ObservableType.one_way_averaged_doppler_type") {
         return R"(
     )";


    } else if(name == "ObservableType.two_way_instantaneous_doppler_type") {
         return R"(
     )";


    } else if(name == "ObservableType.n_way_averaged_doppler_type") {
         return R"(
     )";


    } else if(name == "ObservableType.euler_angle_313_observable_type") {
         return R"(
     )";



    } else if(name == "ObservationViabilityType") {
         return R"(

        Enumeration of observation viability criterion types.





     )";


    } else if(name == "ObservationViabilityType.minimum_elevation_angle") {
         return R"(
     )";


    } else if(name == "ObservationViabilityType.body_avoidance_angle") {
         return R"(
     )";


    } else if(name == "ObservationViabilityType.body_occultation") {
         return R"(
     )";



    } else if(name == "LightTimeFailureHandling") {
         return R"(

        Enumeration of behaviour when failing to converge light-time with required settings.





     )";


    } else if(name == "LightTimeFailureHandling.accept_without_warning") {
         return R"(
     )";


    } else if(name == "LightTimeFailureHandling.print_warning_and_accept") {
         return R"(
     )";


    } else if(name == "LightTimeFailureHandling.throw_exception") {
         return R"(
     )";




    } else if(name == "LinkEndId") {
         return R"(

        Object serving as identifier of a specific link end.





     )";





    } else if(name == "LinkDefinition") {
         return R"(

        Object storing the link ends involved in a given observation.





     )";


    } else if(name == "LinkDefinition.link_ends") {
         return R"(

        Dictionary of link ends, with the key denoting the role in the observaton, and the associated value the identifier for the link end.

        :type: dict[LinkEndType,LinkEndId]
     )";





    } else if(name == "DopplerProperTimeRateSettings") {
         return R"(

        Base class to defining proper time rate settings.

        Functional (base) class for settings of proper time rate (at a single link end) for instantaneous Doppler observation model settings.
        Specific proper time rate settings must be defined using an object derived from this class.
        The derived classes are made accessible via dedicated factory functions.





     )";





    } else if(name == "ObservationSettings") {
         return R"(

        Base class for defining observation settings.

        Functional (base) class for settings of observation models.
        Observation model settings define at least the type and geometry of a given observation.
        They can furthermore set observation biases and/or light-time corrections.
        Simple observation models settings that are fully characterised by these elements can be managed by this base class, which can be instantiated through dedicated factory functions, such as
        :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`, :func:`~tudatpy.numerical_simulation.estimation_setup.observation.cartesian_position`, 
        :func:`~tudatpy.numerical_simulation.estimation_setup.observation.angular_position`, etc.
        Model settings for specific observation models that require additional information such as integration time, retransmission time, etc. must be defined using an object derived from this class.
        The derived classes are made accessible through further factory functions.





     )";





    } else if(name == "OneWayDopplerObservationSettings") {
         return R"(

        Class for defining the settings of one-way instantaneous Doppler observation models.

        :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived class for one-way instantaneous Doppler observation model settings.
        Settings object can account for additional observation model aspects such as light time corrections and proper time rate settings.
        Instances of this class can be created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_instantaneous` factory function.





     )";





    } else if(name == "LightTimeCorrectionSettings") {
         return R"(

        Base class to defining light time correction settings.

        Functional (base) class for settings of light time corrections.
        This class is not used for calculations of corrections, but is used for the purpose of defining the light time correction properties.
        Specific light time correction settings must be defined using an object derived from this class.
        The derived classes are made accessible via dedicated factory functions, such as e.g. 
        :func:`~tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction`





     )";





    } else if(name == "LightTimeConvergenceCriteria") {
         return R"(

        Base class to defining light time convergence criteria.

        Functional (base) class for criteria of light time convergence.
        This class is not used for calculations of corrections, but is used for the purpose of defining the light time convergence criteria.
        Specific light time convergence criteria must be defined using an object derived from this class.
        The derived classes are made accessible via :func:`~tudatpy.numerical_simulation.estimation_setup.observation.light_time_convergence_settings`.





     )";





    } else if(name == "ObservationBiasSettings") {
         return R"(

        Base class to defining observation bias settings.

        Functional (base) class for settings of observation bias.
        Specific observation bias settings must be defined using an object derived from this class.
        The derived classes are made accessible via dedicated factory functions.





     )";





    } else if(name == "ObservationSimulationSettings") {
         return R"(

        Base class for defining settings for simulating observations.

        Base class for defining settings for simulating observations.
        This simulation settings object defines observation times, noise and viability criteria, *etc.* at which observations are to be simulated.
        Therefore, one simulation settings object of this type can only refer to one combination of observable type and link geometry (LinkDefinition).
        The user does not interact with this class directly, but defines specific observation simulation settings using an object derived from this class (created through the associated factory function).





     )";





    } else if(name == "TabulatedObservationSimulationSettings") {
         return R"(

        Class for defining settings for simulating observations at a predefined set of times.

        :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived class for defining settings for simulating observations
        at a predefined set of times
        This type defines predefined time epochs at which applicable observations are to be simulated, stored in a rigid, "tabulated" form. 
        Some observation times may be discarded due to the use of viability settings.
        Instances of this class are typicall created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings` 
        and :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings_list` factory functions. 





     )";





    } else if(name == "ObservationViabilitySettings") {
         return R"(

        Class for defining observation viability calculator settings.

        Class for defining the settings for observation viability calculator creation.
        Instances of this class can be created through various dedicated factory functions, such as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.elevation_angle_viability`, :func:`~tudatpy.numerical_simulation.estimation_setup.observation.body_avoidance_viability` and :func:`~tudatpy.numerical_simulation.estimation_setup.observation.body_occultation_viability`





     )";





    } else if(name == "ObservationDependentVariableSettings") {
         return R"(

        Base class for setting observation dependent variables.

        Functional (base) class for setting observation dependent variables as part of the observation output.
        Note: The associated functionality is not yet mature enough for the end user. Class is exposed for development purposes only.





     )";





    } else if(name == "ObservationAncilliarySimulationSettings") {
         return R"(

        Class for holding ancilliary settings for observation simulation.





     )";




    } else if(name == "ObservationAncilliarySimulationSettings.get_float_settings" && variant==0) {
            return R"(


        Parameters
        ----------
        setting_type : ObservationAncilliarySimulationVariable
            Type of the setting for which the value is to be returned

        throw_exception : bool, default = false
            Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as a floating point value.

        Returns
        -------
        float
            Value of the requested ancilliary variable (or NaN if it does not exist)





    )";



    } else if(name == "ObservationAncilliarySimulationSettings.get_float_list_settings" && variant==0) {
            return R"(


        Parameters
        ----------
        setting_type : ObservationAncilliarySimulationVariable
            Type of the setting for which the value is to be returned

        throw_exception : bool, default = false
            Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as list of floating point values.

        Returns
        -------
        list[ float ]
            Value of the requested ancilliary variable (or empty list if it does not exist)





    )";





    } else if(name == "link_definition" && variant==0) {
        return R"(
        
Function to create a link definition object.


Parameters
----------
link_ends : dict[LinkEndType,LinkEndId]
    Dictionary of link ends, with the key denoting the role in the observaton, and the associated value the identifier for the link end.
Returns
-------
LinkDefinition
    The ``LinkDefinition`` object storing the link ends of the observation






    )";



    } else if(name == "body_origin_link_end_id" && variant==0) {
        return R"(
        
Function to create a link end identifier for the origin (typically center of mass) of a body.

Function to create a link end identifier for the origin (typically center of mass) of a body.
Using this option will simulate the origin of a body transmitter, receiving, etc. the observation.
Although this is typically not physically realistic, it can be a useful approximation, in particular for simulation studies.


Parameters
----------
body_name : str
    Name of the body   

Returns
-------
LinkEndId
    A LinkEndId object representing the center of mass of a body







    )";



    } else if(name == "body_reference_point_link_end_id" && variant==0) {
        return R"(
        
Function to create a link end identifier for a reference point on a body.

Function to create a link end identifier for a reference point on a body, where the reference point
is typically the identifier of a ground stations


Parameters
----------
body_name : str
    Name of the body on which the reference point is located  

body_name : str
    Name of the reference point on the body.  

Returns
-------
LinkEndId
    A LinkEndId object representing a reference point on a body







    )";



    } else if(name == "one_way_downlink_link_ends" && variant==0) {
        return R"(
        
Function for defining one-way downlinks via LinkDefinition types.

Function for defining single or multiple one-way downlinks.
Multiple downlinks share the same transmitters, but may each have a different receiver.
For each downlink, the returned list will contain an additional `LinkDefinition` type.


Parameters
----------
transmitter : Tuple[str, str]
    `LinkEndId` type (tuple of strings), where the first entrance identifies the body and the second entry the reference point of the single transmitter link end.

receivers : List[ Tuple[str, str] ]
    List of `LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the receiver link end(s).

Returns
-------
List[ LinkDefinition ]
    List of one or more `LinkDefinition` types, each defining the geometry for one one-way downlink.
    A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a `LinkEndId` type is mapped.







    )";



    } else if(name == "one_way_uplink_link_ends" && variant==0) {
        return R"(
        
Function for defining one-way uplinks via LinkDefinition types.

Function for defining single or multiple one-way uplinks.
Multiple uplinks share the same receiver, but may each have a different transmitter.
For each uplink, the returned list will contain an additional `LinkDefinition` type.


Parameters
----------
transmitters : List[ Tuple[str, str] ]
    List of `LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the transmitter link end(s).

receiver : Tuple[str, str]
    `LinkEndId` type (tuple of strings), where the first entrance identifies the body and the second entry the reference point of the single receiver link end.

Returns
-------
List[ LinkDefinition ]
    List of one or more `LinkDefinition` types, each defining the geometry for one one-way uplink.
    A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a `LinkEndId` type is mapped.







    )";



    } else if(name == "light_time_convergence_settings" && variant==0) {
        return R"(
        
Factory function for creating convergence settings for solving the light-time equation.

Factory function for creating convergence settings for solving the light-time equation. Computing the light time
:math:`s=t_{R}-t_{T}` between two receiver :math:`R` and transmiter :math:`T` requires the iterative
solution of the following equation:

.. math::
   {t_{R}-t_{T}}=c\left(|\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})| + \Delta s(t_{R},t_{T},mathbf{r}_{R}(t_{R}),mathbf{r}_{T}(t_{T}))

where either the reception time :math:`t_{R}` or the transmission time :math:`t_{T}` is kept fixed (reference link end time). The term :math:`\Delta s` contains any
deviations in the light-time from straight-line propagation at speed of light (relativistic corrections, media corrections, etc.). The algorithm starts
at :math:`t_{R}=t_{T}`, and uses this to evaluate the right-hand side of the above equation. This leads to a new value of :math:`t_{R}` or :math:`t_{T}` (depending on which is kept fixed)
and the right-hand side is re-evaluated in a new iteration. The input to this function defines the settings for when the iteration will terminate.


Parameters
----------
iterate_corrections : bool, default = False
    Boolean denoting whether the terms :math:`\Delta s` are recomputed at each iteration or not. If false, the corrections are calculated only on the first iteration. Afterwards, the value
    is kept fixed until convergence. Once preliminarily converged, the algorithm recomputes :math:`\Delta s`, and continues the iteration (until proper convergence) while now recomputing
    :math:`\Delta s` each iteration. Setting this input to false is typically safe, and is computationally more efficient.

maximum_number_of_iterations : int, default = 50
    Maximum number of iterations taken by the algorithm. If this number of iterations is reached without convergence (as defined by ``absolute_tolerance`` input),
    the behaviour of the algorithm is defined by the ``failure_handling`` input.

absolute_tolerance : float, default = nan
    Difference in :math:`t_{R}-t_{T}` between two consecutive iterataions below which the algorithm is considered to be converged. Default value is nan, which means the default value is taken.
    The default value depends on the time representation used (1 ps for float; 1 fs for Time class)

failure_handling : LightTimeFailureHandling, default = accept_without_warning
    Input defines behaviour when failing to converge within the required number of iterations. NOTE: the default value should be overridden for high-accuracy applications

Returns
-------
:class:`LightTimeConvergenceCriteria`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeConvergenceCriteria` with the required settings.






    )";



    } else if(name == "first_order_relativistic_light_time_correction" && variant==0) {
        return R"(
        
Factory function for creating settings for first-order relativistic light-time corrections.

Factory function for creating settings for first-order relativistic light-time corrections: the correction to
the light time of a (set of) stationary point masses, computed up to c−2 according to general relativity as formulated by e.g. Moyer (2000).
One ambiguity in the model is the time at which the states of the perturbing bodies are evaluated. We distinguish two cases:

* In the case where the perturbing body contains a link end of the observation (for instance perturbation due to Earth gravity field,
  with one of the link ends being an Earth-based station), the time at which the Earth’s state is evaluated equals the transmission time if Earth acts as transmitter, and reception time if
  Earth acts as receiver.
* In other cases, where the perturbing body is not involved in the link ends, its state is evaluated at the midpoint time between transmitter and receiver.


Parameters
----------
perturbing_bodies : str
    A list containing the names of the bodies due to which the light-time correction is to be taken into account.

Returns
-------
:class:`FirstOrderRelativisticLightTimeCorrectionSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` derived :class:`FirstOrderRelativisticLightTimeCorrectionSettings` class,
    defining the settings for the light-time corrections







    )";



    } else if(name == "absolute_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for an absolute observation bias.

Factory function for creating settings for an absolute observation bias. When calculating the observable value, applying this setting
will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

.. math::
   \tilde{h}=h+K

where :math:`K` is the `bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the addition is component-wise.


Parameters
----------
bias_value : numpy.ndarray
    A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
    applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

Returns
-------
:class:`ConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ConstantObservationBiasSettings` class, defining the settings for a constant, absolute observation bias.






    )";



    } else if(name == "relative_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for a relative observation bias.

Factory function for creating settings for a relative observation bias. When calculating the observable value, applying this setting
will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

.. math::
   \tilde{h}=h(1+K)

where :math:`K` is the`bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the multiplication is component-wise.


Parameters
----------
bias_value : numpy.ndarray
    A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
    applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

Returns
-------
:class:`ConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ConstantObservationBiasSettings` class,
    defining the settings for a constant, relative observation bias.







    )";



    } else if(name == "arcwise_absolute_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for arc-wise absolute observation biases.

Factory function for creating settings for arc-wise absolute observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
arc_start_times : List[ float ]
    List containing starting times for each arc.

bias_values : List[ numpy.ndarray ]
    List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
    applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ArcWiseConstantObservationBiasSettings` class.






    )";


    } else if(name == "arcwise_absolute_bias" && variant==1) {
        return R"(
        
Factory function for creating settings for arc-wise absolute observation biases.

Factory function for creating settings for arc-wise absolute observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
    Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
    The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";



    } else if(name == "arcwise_absolute_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for arc-wise absolute observation biases.

Factory function for creating settings for arc-wise absolute observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
arc_start_times : List[ float ]
    List containing starting times for each arc.

bias_values : List[ numpy.ndarray ]
    List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
    applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ArcWiseConstantObservationBiasSettings` class.






    )";


    } else if(name == "arcwise_absolute_bias" && variant==1) {
        return R"(
        
Factory function for creating settings for arc-wise absolute observation biases.

Factory function for creating settings for arc-wise absolute observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
    Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
    The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";



    } else if(name == "arcwise_relative_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for arc-wise relative observation biases.

Factory function for creating settings for arc-wise relative observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
arc_start_times : List[ float ]
    List containing starting times for each arc.

bias_values : List[ numpy.ndarray ]
    List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
    applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";


    } else if(name == "arcwise_relative_bias" && variant==1) {
        return R"(
        
Factory function for creating settings for arc-wise relative observation biases.

Factory function for creating settings for arc-wise relative observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
    Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
    The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";



    } else if(name == "arcwise_relative_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for arc-wise relative observation biases.

Factory function for creating settings for arc-wise relative observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
arc_start_times : List[ float ]
    List containing starting times for each arc.

bias_values : List[ numpy.ndarray ]
    List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
    applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";


    } else if(name == "arcwise_relative_bias" && variant==1) {
        return R"(
        
Factory function for creating settings for arc-wise relative observation biases.

Factory function for creating settings for arc-wise relative observation biases.
This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


Parameters
----------
bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
    Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
    The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

reference_link_end_type : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";



    } else if(name == "time_drift_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for a time-drift bias.

TODO


Parameters
----------
bias_value : numpy.ndarray
    Constant time drift bias that is to be considered for the observation time. This vector should be the same size as the observable to which it is
    assigned (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

time_link_end : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used the current time.

ref_epoch : float
    Defines the reference epoch at which the effect of the time drift is initialised.

Returns
-------
:class:`ConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ConstantObservationBiasSettings` class,
    defining the settings for a constant, relative observation bias.







    )";



    } else if(name == "arc_wise_time_drift_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for arc-wise time-drift biases.

TODO


Parameters
----------
bias_value : numpy.ndarray
    Constant time drift bias that is to be considered for the observation time. This vector should be the same size as the observable to which it is
    assigned (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

arc_start_times : List[ float ]
    List containing starting times for each arc.

time_link_end : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used the current time.

ref_epochs : List[ float ]
    List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ArcWiseConstantObservationBiasSettings` class.






    )";


    } else if(name == "arc_wise_time_drift_bias" && variant==1) {
        return R"(
        
Factory function for creating settings for arc-wise time-drift biases.

TODO


Parameters
----------
bias_value_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
    Dictionary, in which the time bias value vectors for each arc are directly mapped to the starting times of the respective arc.
    The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

time_link_end : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used the current time.

ref_epochs : List[ float ]
    List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";



    } else if(name == "arc_wise_time_drift_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for arc-wise time-drift biases.

TODO


Parameters
----------
bias_value : numpy.ndarray
    Constant time drift bias that is to be considered for the observation time. This vector should be the same size as the observable to which it is
    assigned (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

arc_start_times : List[ float ]
    List containing starting times for each arc.

time_link_end : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used the current time.

ref_epochs : List[ float ]
    List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ArcWiseConstantObservationBiasSettings` class.






    )";


    } else if(name == "arc_wise_time_drift_bias" && variant==1) {
        return R"(
        
Factory function for creating settings for arc-wise time-drift biases.

TODO


Parameters
----------
bias_value_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
    Dictionary, in which the time bias value vectors for each arc are directly mapped to the starting times of the respective arc.
    The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

time_link_end : :class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is used the current time.

ref_epochs : List[ float ]
    List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

Returns
-------
:class:`ArcWiseConstantObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.






    )";



    } else if(name == "combined_bias" && variant==0) {
        return R"(
        
Factory function for creating settings for a combined observation bias.

Factory function for creating settings for a combined observation bias, calculating by combining any number of bias types.
Each contribution of the combined bias is computed from the unbiased observable, so when applying both a relative and absolute bias, we get:

.. math::
   \tilde{h}=h+K_{a}+hK_{r}

And, crucially:

.. math::
   \tilde{h}\neq (h+K_{a})(1+K_{r})

where :math:`K_{r}` and :math:`K_{a}` is the relative and absolute bias, respectively.


Parameters
----------
bias_list : List[ class:`ObservationBiasSettings` ]
    A list containing the bias the bias settings that are to be applied to the observable.

Returns
-------
:class:`MultipleObservationBiasSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`MultipleObservationBiasSettings` class, combining the settings for multiple observation biases.






    )";



    } else if(name == "one_way_range" && variant==0) {
        return R"(
        
Factory function for creating settings for a one-way range observable.

Factory function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
a single-valued observable :math:`h_{_{\text{1-range}}}` as follows (in the unbiased case):

.. math::
   h_{_{\text{1-range}}}(t_{R},t_{T})=|\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})| + \Delta s

where :math:`\mathbf{r}_{R}`, :math:`\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time 
of receiver and transmitter. The term :math:`\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the equation:

.. math::
   t_{R}-t_{T}=c\left(|\mathbf{r}_{R}(t_{R})-\mathbf{r}(t_{R})| + \Delta s\right)

 The method for the iterative solution is described in the :func:`light_time_convergence_settings` entry


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the one-way observable.






    )";



    } else if(name == "n_way_range" && variant==0) {
        return R"(
        
Factory function for creating settings for a n-way range observable.

Factory function for creating observation model settings of n-way range type observables, for a single link definition. The associated observation model creates
a single-valued observable :math:`h_{_{\text{N-range}}}` by combining together a series :math:`n` one-way range observations    
(see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`). By default, the reception time of the :math:`i^{th}` one-way range is set as the 
transmission time of the :math:`(i+1)^{th}` one-way range. A retransmission delay may be defined by ancilliary settings (see TODO) when creating observation
simulation setings.

For this factory function, the settings for each constituent one-way range (with the exception of the link end identifiers) are equal.


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
    as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user). 

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
    Note that only one bias setting is applied to the n-way observable.

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`NWayRangeObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.






    )";



    } else if(name == "n_way_range_from_one_way_links" && variant==0) {
        return R"(
        
Factory function for creating settings for a n-way range observable.

Factory function for creating observation model settings of n-way range type observables, for a single link definition. The
implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with the difference
that the constituent one-way ranges may have different settings.s


Parameters
----------
one_way_range_settings : List[ :class:`ObservationModelSettings` ]
    List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way range observable.
    The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
    ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter``(n-1) and ``receiver`` are defined by the
    ``transmitter`` and ``receiver`` of the :math:`n`^{th} entry of this list.

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
    Note that only one bias setting is applied to the n-way observable.

Returns
-------
:class:`NWayRangeObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.






    )";



    } else if(name == "two_way_range" && variant==0) {
        return R"(
        
Factory function for creating settings for a two-way range observable.

Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with :math:`n=2`. This function is provided
for convenience.


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter`, `retransmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
    Note that only one bias setting is applied to the n-way observable.

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`NWayRangeObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.






    )";



    } else if(name == "two_way_range_from_one_way_links" && variant==0) {
        return R"(
        
Factory function for creating settings for a two-way range observable.

Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_from_one_way_links`, with :math:`n=2`. This function is provided
for convenience.


Parameters
----------
one_way_range_settings : List[ :class:`ObservationModelSettings` ]
    List of observation model settings of size two, with the first entry the one-way range settings for the uplink, and the second entry the one-way range settings for the downlink.
    The ``LinkDefinition`` of this two-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
    ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` and ``receiver`` are defined by the
    ``transmitter`` and ``receiver`` of the second entry of this list.

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
    Note that only one bias setting is applied to the n-way observable.

Returns
-------
:class:`NWayRangeObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`NWayRangeObservationSettings` class.






    )";



    } else if(name == "angular_position" && variant==0) {
        return R"(
        
Factory function for creating settings for an angular position observable.

Factory function for creating observation model settings of angular position type observables (as right ascension :math:`\alpha` and declination :math:`\delta`), 
for a single link definition. The associated observation model creates an observable :math:`\mathbf{h}_{_{\text{ang.pos.}}}` of type two as follows (in the unbiased case):

.. math::
   \Delta\mathbf{r}=\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})\\
   \tan\alpha=\frac{\Delta r_{y}}{\Delta r_{x}}\\
   \delta=\frac{\Delta r_{z}}{\sqrt{\Delta r_{x}^{2}+\Delta r_{y}^{2}}}\\
   \mathbf{h}_{_{\text{ang.pos.}}} = [\alpha;\delta]

The relative position vector :math:`\Delta\mathbf{r}` is computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`
The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the 
environment.


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the angular position observable.






    )";



    } else if(name == "relative_angular_position" && variant==0) {
        return R"(
        
Factory function for creating settings for a relative angular position observable.

Factory function for creating observation model settings of relative angular position type observables (as relative right ascension :math:`\Delta\alpha` and relative declination :math:`\Delta\delta`), 
for a single link definition. The associated observation model creates an observable that is the difference of two :func:`~tudatpy.numerical_simulation.estimation_setup.observation.angular_position`
observables :math:`\left(\mathbf{h}_{_{\text{ang.pos.}}}\right)_{2}` and :math:`\left(\mathbf{h}_{_{\text{ang.pos.}}}\right)_{1}`. 
The resulting observable becomes :math:`\left(\mathbf{h}_{_{\text{ang.pos.}}}\right)_{2}-\left(\mathbf{h}_{_{\text{ang.pos.}}}\right)_{2}`


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter`, 'transmitter2` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.
    The `transmitter` and `transmitter2` entries are used to define separate angular position observables (as observed by the `receiver`), from
    which the observable is then computed

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the angular position observable.






    )";



    } else if(name == "one_way_doppler_instantaneous" && variant==0) {
        return R"(
        
Factory function for creating settings for a one-way instantaneous Doppler observable.

Factory function for creating settings for a one-way instantaneous Doppler observable for a single link definition. The associated observation model creates
a single-valued observable :math:`h_{_{\text{1-Dopp.}}}` as follows (in the unbiased case):

.. math::
   h_{_{\text{1-Dopp.}}}=c\left(\frac{d\tau_{T}}{dt_{T}}\frac{t_{T}}{dt_{R}}\frac{dt_{R}}{d\tau_{R}}-1\right)
   
where :math:`t` and :math:`\tau` denote coordinate and proper time of the transmitter T and receiver R, respectively.
The receiver and transmitter position and coordinate time are computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`.
The detailed mathematical implementation are described on TODO.

This observable represents the 'instantaneous (non-integrated)' Doppler observable, as obtained from open-loop observations.
It should *not* be used for the modelling of the typical closed-loop observations used in deep space tracking (for which the 
:func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` should be used)

The coordinate
time derivative :math:`\frac{t_{A}}{dt_{B}}` is always computed when generating this observable. Settings for the proper time
rates :math:`\frac{d\tau}{dt}` can be specified by the user through the ``transmitter_proper_time_rate_settings`` and ``receiver_proper_time_rate_settings``
inputs. Whene these are left empty, the proper time rates are omitted (set to 1.0).

The observable may be non-dimensionalized by the speed of light :math:`c`, which results in the observable being equal to thee received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires that the
    `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

transmitter_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
    Settings for computing the transmitter proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

receiver_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
    Settings for computing the receiver proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

normalized_with_speed_of_light : bool, default = false
    Option to non-dimensionalize the observable with speed of light :math:`c`

Returns
-------
:class:`OneWayDopplerObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`OneWayDopplerObservationSettings` class defining the settings for the one-way open doppler observable observable.






    )";



    } else if(name == "two_way_doppler_instantaneous" && variant==0) {
        return R"(
        
Factory function for creating settings for a two-way instantaneous Doppler observable.


Factory function for creating settings for a two-way instantaneous Doppler observable for a single link definition. The associated observation model creates
a single-valued observable :math:`h_{_{\text{2-Dopp.}}}` as follows (in the unbiased case), by combining the up- and downlink one-way 
instantanenous Doppler observable :math:`h_{_{\text{1-Dopp.}}}` 
(both normalized by speed of light c, see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_instantaneous`) as follows:

.. math::
   h_{_{\text{2-Dopp.}}}=c\left((\left(h_{_{\text{1-Dopp.}}}\right)_{\text{up}}+1)(\left(h_{_{\text{1-Dopp.}}}\right)_{\text{down}}+1)-1\right)

The resulting observable is non-dimensional (but can be converted to an observed range-rate by multiplying with speed of light :math:`c`)

This observable represents the 'instantaneous (non-integrated)' Doppler observable, as obtained from open-loop observations.
It should *not* be used for the modelling of the typical closed-loop observations used in deep space tracking (for which the 
:func:`~tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_averaged` should be used)

The observable may be non-dimensionalized by the speed of light :math:`c`,
which results in the observable being equal to the received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter`, `retransmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
    Note that only one bias setting is applied to the n-way observable.s

transmitter_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
    Settings for computing the transmitter proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

receiver_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
    Settings for computing the receiver proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

normalized_with_speed_of_light : bool, default = false
    Option to non-dimensionalize the observable with speed of light :math:`c`

Returns
-------
:class:`TwoWayDopplerObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`TwoWayDopplerObservationSettings` class defining the settings for the two-way open doppler observable.






    )";



    } else if(name == "two_way_doppler_instantaneous_from_one_way_links" && variant==0) {
        return R"(
        
Factory function for creating settings for a two-way instantaneous Doppler observable.


Factory function for creating settings for a two-way instantaneous Doppler observable for a single link definition. The
implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_instantaneous`, with the difference
that the constituent one-way ranges may have different settings.
   
The observable may be non-dimensionalized by the speed of light :math:`c` (in the constituent one-way Doppler observable settings),
which results in the observable being equal to the received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.


Parameters
----------
uplink_doppler_settings : :class:`OneWayDopplerObservationSettings`
    Settings for uplink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`

downlink_doppler_settings : :class:`OneWayDopplerObservationSettings`
    Settings for downlink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the full observation, default is none (unbiased observation). Note that,
    even if no bias is applied to the two-way observable, the constituent one-way observables may still be biased.

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`TwoWayDopplerObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived :class:`TwoWayDopplerObservationSettings` class defining the settings for the two-way open doppler observable.






    )";



    } else if(name == "one_way_doppler_averaged" && variant==0) {
        return R"(
        
Factory function for creating settings for a one-way averaged Doppler observable.

Factory function for creating observation model settings for one-way averaged Doppler observables, for a single link definition. The associated observation model creates
a single-valued observable :math:`h_{_{\text{1-\bar{Dopp}}}}` as follows (in the unbiased case):

.. math::
   h_{_{\text{1-\bar{Dopp}}}}&=c\int_{t-\Delta t}^{t+\Delta t}\frac{t_{T}}{dt_{R}}d\bar{t}\\
                             &=\frac{h_{_{\text{1-range}}}(t_{R}=t+\Delta t,t_{T})-h_{_{\text{1-range}}}(t_{R}=t,t_{T})}{\Delta t}

where, in the latter formulation (which is the one that is implemented), the observable is referenced to the receiver time. This averaged Doppler observable
is computed as the difference of two one-way range observables (see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`), 
with the reference time shifted by :math:`\Delta t`. As such, it is sensitive to numerical errors for small :math:`\Delta t`

The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires that the
    `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `OneWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






    )";



    } else if(name == "n_way_doppler_averaged" && variant==0) {
        return R"(
        
Factory function for creating settings for an n-way averaged Doppler observable.

Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
    as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user). 

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






    )";



    } else if(name == "n_way_doppler_averaged_from_one_way_links" && variant==0) {
        return R"(
        
Factory function for creating settings for an n-way averaged Doppler observable.

Factory function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. 
The implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged`, with the difference
that the constituent one-way range observables may have different settings.


Parameters
----------
one_way_range_settings : List[ :class:`ObservationModelSettings` ]
    List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way averaged range rate observable.
    The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
    ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter``(n-1) and ``receiver`` are defined by the
    ``transmitter`` and ``receiver`` of the :math:`n`^{th} entry of this list.

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






    )";



    } else if(name == "two_way_doppler_averaged" && variant==0) {
        return R"(
        
Factory function for creating settings for a two-way averaged Doppler observable.

Factory function for creating settings for a two-way averaged Doppler observable. Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged`,       
with :math:`n=2`. This function is provided for convenience.


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires the
    `transmitter`, `retransmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined

light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
    List of corrections for the light-time that are to be used. Default is none, which will result
    in the signal being modelled as moving in a straight line with the speed of light

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
    Settings for convergence of the light-time 

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






    )";



    } else if(name == "two_way_doppler_averaged_from_one_way_links" && variant==0) {
        return R"(
        
Factory function for creating settings for a two-way averaged Doppler observable.

Factory function for creating settings for a two-way averaged Doppler observable. Same as 
:func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged_from_one_way_links`,       
with :math:`n=2`. This function is provided for convenience.


Parameters
----------
one_way_range_settings : List[ :class:`ObservationModelSettings` ]
    List of observation model settings of size two, with the first entry the one-way range settings for the uplink, and the second entry the one-way range settings for the downlink.
    The ``LinkDefinition`` of this two-way range observable is created from this list, with the ``transmitter`` and ``retransmitter1`` defined by the
    ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` and ``receiver`` are defined by the
    ``transmitter`` and ``receiver`` of the second entry of this list.

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` derived `NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






    )";



    } else if(name == "cartesian_position" && variant==0) {
        return R"(
        
Factory function for creating settings for a Cartesian position observable.

Factory function for creating observation model settings of Cartesian position type observables.
Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires that the
    `observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)	

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian position observable.






    )";



    } else if(name == "cartesian_velocity" && variant==0) {
        return R"(
        
Factory function for creating settings for a Cartesian velocity observable.

Factory function for creating observation model settings of Cartesian position type observables.
Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
This observable provides the inertial (w.r.t. global frame origin) Cartesian velocity of the `observed_body` defined by the `link_ends` input.
The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` velocity


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires that the
    `observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)	

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the cartesian velocity observable.






    )";



    } else if(name == "313_euler_angles" && variant==0) {
        return R"(
        
Factory function for creating settings for observable containing the body orientation as Euler angles.

Factory function for creating observation model settings of Euler angle type observables.
This observable can be used for *e.g.* body attitude observations, but can also be very useful as 'synthetic' observable for verification or analysis purposes.
This observable provides the rotation from inertial (defined by the global frame orientation) to body-fixed orientation of the 
body specified by the `observed_body` in the `link_ends` input.  The observable
has size 3, and contains the  3-1-3 (e.g. z-x-z) Euler angles


Parameters
----------
link_ends : LinkDefinition
    Set of link ends that define the geometry of the observation. This observable requires that the
    `observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

bias_settings : :class:`ObservationBiasSettings`, default = None
    Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)	

Returns
-------
:class:`ObservationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationModelSettings` class defining the settings for the Euler angle observable.






    )";



    } else if(name == "elevation_angle_viability" && variant==0) {
        return R"(
        
Factory function for defining single elevation angle viability setting.

Factory function for defining elevation angle viability settings for single link end.
When simulating observations, this setting ensures that any applicable observations, for which the local elevation angle at link end is less than some limit value, will be omitted.


Parameters
----------
link_end_id : Tuple[str,str]
    Link end (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
    To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

elevation_angle : float
    Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
    value must be in radians.

Returns
-------
:class:`ObservationViabilitySettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` class, defining the settings for observation viability






    )";



    } else if(name == "elevation_angle_viability_list" && variant==0) {
        return R"(
        
Factory function for defining list of elevation angle viability settings.

Factory function for defining elevation angle viability settings for multiple link ends.
Each entry in the returned list contains the observation viability settings for one link end.
When simulating observations, these settings ensure that any applicable observations, for which the local elevation angle at a link end is less than some limit value, will be omitted.


Parameters
----------
link_end_ids : List[ Tuple[str,str] ]
    List of individual link ends (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
    To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
    For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
    link end violates the minimum elevation angle constraint.

elevation_angle : float
    Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
    value must be in radians.

Returns
-------
:class:`ObservationViabilitySettings`
    List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






    )";



    } else if(name == "body_avoidance_viability" && variant==0) {
        return R"(
        
Factory function for defining body avoidance observation viability settings.

Factory function for defining body avoidance observation viability settings for single link ends.
When simulating observations, this settings ensures that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
The definition of 'too close' is computed as the angle between:

* The line-of-sight vector from a link end to a given third body
* The line-of-sight between two link ends 

This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


Parameters
----------
link_end_id : Tuple[str,str]
    Link end (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
    To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
    For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
    link end passes too close to the specified body.

body_to_avoid : str
    Name of the body which the signal path should not pass 'too close' to.

avoidance_angle : float
    Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
    value must be in radians.

Returns
-------
:class:`ObservationViabilitySettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.






    )";



    } else if(name == "body_avoidance_viability_list" && variant==0) {
        return R"(
        
Factory function for defining list of body avoidance viability settings.

Factory function for defining body avoidance viability settings for multiple link ends.
Each entry in the returned list contains the observation viability settings for one link end.
When simulating observations, these settings ensure that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
The definition of 'too close' is computed as the angle between:

* The line-of-sight vector from a link end to a given third body
* The line-of-sight between two link ends

This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


Parameters
----------
link_end_ids : List[ Tuple[str,str] ]
    List of individual link ends (as defined by body/reference point pair, see TODO), for which the elevation angle viability setting is to be created.
    To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

body_to_avoid : str
    Name of the body which the signal path should not pass 'too close' to.

avoidance_angle : float
    Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
    value must be in radians.

Returns
-------
:class:`ObservationViabilitySettings`
    List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






    )";



    } else if(name == "body_occultation_viability" && variant==0) {
        return R"(
        
Factory function for defining body occultation viability settings.

Factory function for defining body occultation viability settings for single link ends.
When simulating observations, this setting ensures that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
The occultation is computed using the shape model of the specified body.


Parameters
----------
link_end_id : Tuple[str,str]
    Link end (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
    To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.

body_to_avoid : str
    Name of the body which the signal path should not be occulted by.

Returns
-------
:class:`ObservationViabilitySettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.






    )";



    } else if(name == "body_occultation_viability_list" && variant==0) {
        return R"(
        
Factory function for defining body occultation viability settings.

Factory function for defining body occultation viability settings for multiple link ends.
Each entry in the returned list contains the observation viability settings for one link end.
When simulating observations, these settings ensure that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
The occultation is computed using the shape model of the specified body.


Parameters
----------
link_end_ids : List[ Tuple[str,str] ]
    List of individual link ends (as defined by body/reference point pair, see TODO), for which the viability settings are to be created.
    To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
    For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
    link end is occulted by the specified body.

body_to_avoid : str
    Name of the body which the signal path should not be occulted by.

Returns
-------
:class:`ObservationViabilitySettings`
    List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






    )";



    } else if(name == "doppler_ancilliary_settings" && variant==0) {
        return R"(
        
Factory function for creating ancilliary settings for averaged Doppler observable.

Factory function for creating ancilliary settings for an averaged Doppler observable. Specifically, this
function can be used to create settings for the integration time of the observable. Note: in case no retransmission
delays (or other additional ancilliary settings) are to be defined, this setting may be used for one-, two-, or N-way
averaged Doppler.


Parameters
----------
integration_time : float, default = 60.0
    Integration time that is to be used for the averaged Doppler observable
Returns
-------
:class:`ObservationAncilliarySimulationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






    )";



    } else if(name == "two_way_range_ancilliary_settings" && variant==0) {
        return R"(
        
Factory function for creating ancilliary settings for two-way range observable.

Factory function for creating ancilliary settings for a two-way range observable. Specifically, this
function can be used to create settings for the retransmission delay of the observable. NOTE:
this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_ancilliary_settings`
with a single retransmission delay.


Parameters
----------
retransmission_delay : float, default = 0.0
    Retransmission delay that is to be applied to the simulation of the two-way observable
Returns
-------
:class:`ObservationAncilliarySimulationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






    )";



    } else if(name == "two_way_doppler_ancilliary_settings" && variant==0) {
        return R"(
        
Factory function for creating ancilliary settings for two-way averaged Doppler observable.

Factory function for creating ancilliary settings for a two-way range observable. Specifically, this
function can be used to create settings for the retransmission delay of the observable.  NOTE:
this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_ancilliary_settings`
with a single retransmission delay.


Parameters
----------
integration_time : float, default = 60.0
    Integration time that is to be used for the averaged Doppler observable
retransmission_delay : float, default = 0.0
    Retransmission delay that is to be applied to the simulation of the two-way observable
Returns
-------
:class:`ObservationAncilliarySimulationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






    )";



    } else if(name == "n_way_range_ancilliary_settings" && variant==0) {
        return R"(
        
Factory function for creating ancilliary settings for n-way range observable.

Factory function for creating ancilliary settings for a n-way range observable. Specifically, this
function can be used to create settings for the retransmission delays of the observable, for each of the retransmitters.


Parameters
----------
retransmission_delays : list[ float ], default = None
    Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
Returns
-------
:class:`ObservationAncilliarySimulationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






    )";



    } else if(name == "n_way_doppler_ancilliary_settings" && variant==0) {
        return R"(
        
Factory function for creating ancilliary settings for n-way averaged Doppler observable.

Factory function for creating ancilliary settings for a n-way averaged Doppler observable. Specifically, this
function can be used to create settings for the integration time of the observable, and the  retransmission delays for each of the retransmitters.


Parameters
----------
integration_time : float, default = 60.0
    Integration time that is to be used for the averaged Doppler observable
retransmission_delays : list[ float ], default = None
    Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
Returns
-------
:class:`ObservationAncilliarySimulationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






    )";



    } else if(name == "tabulated_simulation_settings" && variant==0) {
        return R"(
        
Factory function for creating settings object for observation simulation, using a predefined list of observation times.

Factory function for creating single simulation settings object, using a predefined list of observation times.
The list of resulting observations may be reduced compared to the ``simulation_times`` provided here, as
only observations that meet the viability settings are retained during observation simulation (these may be
provide directly here through the ``viability_settings`` input, or added later to the resulting settings object). 


Parameters
----------
observable_type : :class:`ObservableType`
    Observable type of which observations are to be simulated.
link_ends : LinkDefinition
    Link ends for which observations are to be simulated.
simulation_times : List[float]
    List of times at which to perform the observation simulation.
reference_link_end_type : :class:`LinkEndType`, default = :class:`LinkEndType`.receiver
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference time for the observation.
viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
    Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.

noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
    Function providing the observation noise factors as a function of observation time.
Returns
-------
:class:`TabulatedObservationSimulationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.






    )";



    } else if(name == "tabulated_simulation_settings_list" && variant==0) {
        return R"(
        
Factory function for creating a list of settings object for observation simulation, using a predefined list of observation times.

Factory function for creating multiple tabulated observation simulation settings objects in a list. This function is 
equivalent to calling the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings` repeatedly, with the different 
observables and link definition provided here through `link_ends_per_observable`. 
During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.


Parameters
----------
link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
    Link geometry per observable type of which observations are to be simulated.
simulation_times : List[ float ]
    List of times at which to perform the observation simulation.
reference_link_end_type : :class:`LinkEndType`, default = LinkEndType.receiver
    Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
    The single link end specified here will be considered as the reference link end for all simulation settings object created in the function call.

viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
    Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
    The single settings list given here will be considered as potential viability settings for all simulation settings object created in the function call.

Returns
-------
List[ TabulatedObservationSimulationSettings ]
    List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` objects.






    )";



    } else if(name == "get_default_reference_link_end" && variant==0) {
        return R"(
        
Factory function for automatically retrieving the reference link end associated with a given observable type.


Parameters
----------
observable_type : :class:`ObservableType`
    Observable type for which the associated reference link end is to be retrieved.
Returns
-------
:class:`LinkEndType`
    Defines the link end (via the :class:`LinkEndType`) which is typically used as a reference for observation times in *e.g.* :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`.






    )";



    } else if(name == "continuous_arc_simulation_settings" && variant==0) {
        return R"(
        
Factory function for creating settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.

Factory function for creating settings object for observation simulation. Unlike the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`
function, the resulting settings do not define the observation times explicitly. Instead, this settings object determines the observation times adaptively during the
simulation of the observation, with the requirement that observations should be simulated over a set of contiguous arcs (if possible). The exact algorithm meets the following conditions:

* Observations are only simulated within the time span of ``start_time`` and ``end_time``
* A contiguous tracking arc has simulated observations separated by ``interval_between_observations``
* Starting from ``start_time``, an observation is simulated each ``interval_between_observations``. Once an observation is unviable, as defined by
  the ``arc_limiting_constraints`` input, it is checked whether the arc up until that point 
  is longer in duration than ``minimum_arc_duration``. If it is, the arc is added to the simulated observations. If not, the arc is discarded. In either case, a new arc is started once a 
  viable is observation is encountered
* If the current arc reaching a duration greater than ``maximum_arc_duration``, the arc is added to the existing observations, and a new arc is started
* If defined (e.g. if not NaN), the current observation time is incremented by ``minimum_time_between_arcs`` when an arc has been added to the observations.

Nominally, this algorithm ensures that any arc of observations has a minimum and maximum duration. In addition, it ensures that (if desired) there is a minimum time interval 
between two tracking arcs. This behaviour can be modified by adding ``additional_viability_settings``, which are *not* used when computing the tracking arcs, but which are instead only used
to reduce the set of simulated observations afterwards.


Parameters
----------
observable_type : :class:`ObservableType`
    Observable type of which observations are to be simulated.
link_ends : LinkDefinition
    Link ends for which observations are to be simulated.
start_time : float
    First time at which an observation is to be simulated (and checked for viability).
end_time : float
    Maximum time at which an observation is to be simulated (and checked for viability).
interval_between_observations : float
    Cadence (in seconds) of subsequent observations in an arc
arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
    List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
    whether an arc should be terminated.

minimum_arc_duration : float
    Minimum permissible time for a tracking arc
maximum_arc_duration : float
    Maximum permissible time for a tracking arc
minimum_time_between_arc : float, default = NaN
    Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
    Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
    These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.

noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
    Function providing the observation noise factors as a function of observation time.
Returns
-------
:class:`TabulatedObservationSimulationSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.






    )";



    } else if(name == "continuous_arc_simulation_settings_list" && variant==0) {
        return R"(
        
Factory function for creating a list of settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.

Factory function for creating multiple settings objects for observation simulation in a list. This function is 
equivalent to calling the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.continuous_arc_simulation_settings` repeatedly, with the different 
observables and link definition provided here through `link_ends_per_observable`. 
During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.


Parameters
----------
link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
    Link geometry per observable type of which observations are to be simulated.
start_time : float
    First time at which an observation is to be simulated (and checked for viability).
end_time : float
    Maximum time at which an observation is to be simulated (and checked for viability).
interval_between_observations : float
    Cadence (in seconds) of subsequent observations in an arc
arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
    List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
    whether an arc should be terminated.

minimum_arc_duration : float
    Minimum permissible time for a tracking arc
maximum_arc_duration : float
    Maximum permissible time for a tracking arc
minimum_time_between_arc : float, default = NaN
    Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
    Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
    These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.

Returns
-------
List[ :class:`TabulatedObservationSimulationSettings` ]
    List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` objects.






    )";



    } else if(name == "add_gaussian_noise_to_all" && variant==0) {
        return R"(
        
Function for adding gaussian noise function to all existing observation simulation settings.

Function for including simple time-independent and time-uncorrelated Gaussian noise function to the simulation settings of one or more observable(s).
The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings` 
list.

Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function, 
and thus the function does not return anything.


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
noise_amplitude : float
    Standard deviation defining the un-biased Gaussian distribution for the noise.
Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else if(name == "add_gaussian_noise_to_observable" && variant==0) {
        return R"(
        
Function for adding gaussian noise function to existing observation simulation settings of a given observable type.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the 
`observation_simulation_settings` list that matches the specified `observable_type`.  


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
noise_amplitude : float
    Standard deviation defining the un-biased Gaussian distribution for the noise.
observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
    Identifies the observable type in the observation simulation settings to which the noise is to be added.

Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else if(name == "add_gaussian_noise_to_observable_for_link_ends" && variant==0) {
        return R"(
        
Function for adding gaussian noise function to existing observation simulation settings of a given observable type and link definition.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the 
`observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.  


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
noise_amplitude : float
    Standard deviation defining the un-biased Gaussian distribution for the noise.
observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
    Identifies the observable type in the observation simulation settings to which the noise is to be added.

link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
    Identifies the link definition in the observation simulation settings for which the noise is to be added.

Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else if(name == "add_viability_check_to_all" && variant==0) {
        return R"(
        
Function for including viability checks into existing observation simulation settings.

Function for adding viability checks to the observation simulation settings, such that only observations meeting certain conditions are retained.
The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings` 
list.      
Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function, 
and thus the function does not return anything.


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
viability_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` ]
    List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else if(name == "add_viability_check_to_observable" && variant==0) {
        return R"(
        
Function for including viability checks into existing observation simulation settings.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_all`, except that the function only adds viabilitt settings to entries of the 
`observation_simulation_settings` list that matches the specified `observable_type`.  


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
viability_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` ]
    List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

observable_type : :class:`ObservableType`
    Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else if(name == "add_viability_check_to_observable_for_link_ends" && variant==0) {
        return R"(
        
Function for including viability checks into existing observation simulation settings.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_all`, except that the function only adds noise to entries of the 
`observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`. 


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
viability_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` ]
    List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

observable_type : :class:`ObservableType`
    Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
    Identifies the link definition in the observation simulation settings for which the viability checks are to be considered.

Returns
-------
None
    The :class







    )";



    } else if(name == "add_dependent_variables_to_all" && variant==0) {
        return R"(
        
Function for including dependent variables into all existing observation simulation settings.

Function for including the computation and reporting of dependent variables into the observation simulation settings of all observables.
Note: The associated functionality is not yet mature enough for the end user. Function is exposed for development purposes only.

Modifications are applied to all given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s),
matching each :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object with the corresponding :class:`ObservationDependentVariableSettings` entry in the `dependent_variable_settings` parameter.
Note that the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place and thus the function does not return anything.


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

dependent_variable_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` ]
    List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.

bodies : :class:`~tudatpy.numerical_simulation.environment_setup.SystemOfBodies`
    Object consolidating all bodies and environment models that constitute the physical environment.






    )";



    } else if(name == "add_dependent_variables_to_observable" && variant==0) {
        return R"(
        
Function for including dependent variables into selected existing observation simulation settings.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_dependent_variables_to_all`, except that the function only adds includes the 
computation and reporting of dependent variables to entries of the `observation_simulation_settings` list that matches the specified `observable_type`.


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

dependent_variable_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` ]
    List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.

bodies : :class:`~tudatpy.numerical_simulation.environment_setup.SystemOfBodies`
    Object consolidating all bodies and environment models that constitute the physical environment.

observable_type : :class:`ObservableType`
    Identifies the observable type in the observation simulation settings for which the dependent variables are to be included.






    )";



    } else if(name == "add_dependent_variables_to_obs_for_links_end" && variant==0) {
        return R"(
        
Function for including dependent variables into selected existing observation simulation settings for the chosen link ends.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_dependent_variables_to_all`, except that the function only adds includes the 
computation and reporting of dependent variables to entries of the `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

dependent_variable_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` ]
    List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.

bodies : :class:`~tudatpy.numerical_simulation.environment_setup.SystemOfBodies`
    Object consolidating all bodies and environment models that constitute the physical environment.

observable_type : :class:`ObservableType`
    Identifies the observable type in the observation simulation settings for which the dependent variables are to be included.

link_ends : LinkDefinition
    Link ends for which the dependent variables are to be included.






    )";



    } else if(name == "add_noise_function_to_all" && variant==0) {
        return R"(
        
Function for adding a custom noise function to all existing observation simulation settings.

Function for including a custom noise function to the simulation settings of all observables.
The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings` 
list.

Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function, 
and thus the function does not return anything.


Parameters
----------
observation_simulation_settings_list : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
    Function providing the observation noise factors as a function of observation time.

Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else if(name == "add_noise_function_to_observable" && variant==0) {
        return R"(
        
Function for adding a custom noise function to selected existing observation simulation settings of a given observable type.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_all`, except that the function only adds noise to entries of the 
`observation_simulation_settings` list that matches the specified `observable_type`.


Parameters
----------
observation_simulation_settings_list : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
    Function providing the observation noise factors as a function of observation time.

observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
    Identifies the observable type in the observation simulation settings to which the noise is to be added.

Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else if(name == "add_noise_function_to_observable_for_link_ends" && variant==0) {
        return R"(
        
Function for adding a custom noise function to existing observation simulation settings of a given observable type and link definition.

As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_all`, except that the function only adds noise to entries of the 
`observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.  


Parameters
----------
observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
    Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
    Function providing the observation noise factors as a function of observation time.

observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
    Identifies the observable type in the observation simulation settings to which the noise is to be added.

link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
    Identifies the link definition in the observation simulation settings for which the noise is to be added.

Returns
-------
None
    The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace parameter {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "EstimatableParameterTypes") {
         return R"(

        Enumeration of model parameters that are available for estimation.
        In order to establish a parameter estimation settings for a parameter of a certain type, use the factory function dedicated to this parameter type.
        Note that not all of the listed types might be accessible via factory functions in the python interface yet.






     )";


    } else if(name == "EstimatableParameterTypes.arc_wise_initial_body_state_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.initial_body_state_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.initial_rotational_body_state_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.constant_drag_coefficient_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.arc_wise_constant_drag_coefficient_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.radiation_pressure_coefficient_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.arc_wise_radiation_pressure_coefficient_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.empirical_acceleration_coefficients_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.arc_wise_empirical_acceleration_coefficients_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.desaturation_delta_v_values_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.gravitational_parameter_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.spherical_harmonics_cosine_coefficient_block_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.spherical_harmonics_sine_coefficient_block_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.mean_moment_of_inertia_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.constant_rotation_rate_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.rotation_pole_position_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.polar_motion_amplitude_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.core_factor_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.free_core_nutation_rate_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.periodic_spin_variation_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.constant_additive_observation_bias_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.arcwise_constant_additive_observation_bias_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.constant_relative_observation_bias_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.arcwise_constant_relative_observation_bias_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.ground_station_position_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.full_degree_tidal_love_number_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.single_degree_variable_tidal_love_number_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.direct_dissipation_tidal_time_lag_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.ppn_parameter_gamma_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.ppn_parameter_beta_type") {
         return R"(
     )";


    } else if(name == "EstimatableParameterTypes.equivalence_principle_lpi_violation_parameter_type") {
         return R"(
     )";




    } else if(name == "EstimatableParameterSettings") {
         return R"(

        Base class to defining settings of parameter to be estimated.

        Functional (base) class for settings of model parameter to be estimated.
        Settings of simple parameters types are managed via this class, more complex parameter types are handled by specialised derivates of this class.
        Instances of either base or derived class can be created via dedicated factory functions.





     )";






    } else if(name == "initial_states" && variant==0) {
        return R"(
        
Function for defining parameter settings for initial state parameters.

Factory function for creating a (linear sensitivity) parameter settings object for initial state parameters.
The factory function uses the propagator settings to determine which type of initial state parameter (single/multi/hybrid-arc; translational/rotational/... dynamics) is to be estimated,
e.g. if a single-arc translational state propagator is defined, the function will automatically create the parameters for the associated initial state parameter

.. note:: These function return lists of parameter settings objects.
This means that which the return of this function cannot simply be added to the parameter settings objects of single parameters in a list creation statement.
Instead, list concatenation is recommended. Please see the following example:

.. code-block:: python 
   
   # define single estimatable parameters 
   single_parameter_1 = ...  
   single_parameter_2 = ...  
   ...  

   # bad: list creation statement --> will result in nested list, undesired!
   list_of_all_parameters = [estimation_setup.parameter.initial_states(...), single_parameter_1, single_parameter_2, ...] 

   # better: list concatenation --> will result in simple list, desired!
   list_of_all_parameters = estimation_setup.parameter.initial_states(...) + [single_parameter_1, single_parameter_2, ...] 


Parameters
----------
propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
    Object containing the consolidated propagation settings of the simulation in the context of which the given model parameters are to be estimated.

bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies
    Object consolidating all bodies and environment models that constitute the physical environment.

arc_initial_times : List[ float ] = []
    Initial times of arcs, only required if arc-wise propagation settings are passed via the `propagator_settings` object.

Returns
-------
List[ :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` ]
    List of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` objects, one per component of each initial state in the simulation.







    )";



    } else if(name == "constant_drag_coefficient" && variant==0) {
        return R"(
        
Function for defining parameter settings for constant drag coefficients.

Factory function for creating a (linear sensitivity) parameter settings object for constant drag coefficients.
Using the constant drag coefficient as an estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
* The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.aerodynamic` acceleration


Parameters
----------
body : str
    Name of the body, with whose drag acceleration model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's constant drag coefficient.







    )";



    } else if(name == "arcwise_constant_drag_coefficient" && variant==0) {
        return R"(
        
Function for defining parameter settings for arc-wise constant drag coefficients.

Factory function for creating (linear sensitivity) parameter settings object for arc-wise constant drag coefficients (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.constant_drag_coefficient`).
Using the arc-wise constant drag coefficient as an estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
* The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.aerodynamic` acceleration

.. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the drag coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.  


Parameters
----------
body : str
    Name of the body, with whose drag acceleration model the estimatable parameter is associated.

arc_initial_times : List[ float ]
    List of times at which the arcs over which the drag coefficient is to be estimated will start.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseDragCoefficientEstimatableParameterSettings` class
    for arc-wise treatment of the specified body's constant drag coefficient.







    )";



    } else if(name == "radiation_pressure_coefficient" && variant==0) {
        return R"(
        
Function for defining parameter settings for radiation pressure coefficients.

Factory function for creating a (linear sensitivity) parameter settings object for a radiation pressure coefficient.
Using the radiation pressure coefficient as an estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.cannonball` radiation pressure interface to be defined for the body specified by the ``body`` parameter
* The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.cannonball_radiation_pressure` acceleration


Parameters
----------
body : str
    Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's radiation pressure coefficient.







    )";



    } else if(name == "arcwise_radiation_pressure_coefficient" && variant==0) {
        return R"(
        
Function for defining parameter settings for arc-wise radiation pressure coefficients.

Factory function for creating a (linear sensitivity) parameter settings object for arc-wise radiation pressure coefficients (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.radiation_pressure_coefficient`).
Using the radiation pressure coefficient as an estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.cannonball` radiation pressure interface to be defined for the body specified by the ``body`` parameter
* The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.cannonball_radiation_pressure` acceleration
The radiation pressure coefficient is defined according to the universal convention for a cannonball model and thus no further model information is given.

.. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.  


Parameters
----------
body : str
    Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.

arc_initial_times : List[ float ]
    List of times at which the arcs over which the radiation pressure coefficient is to be estimated will start.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseRadiationPressureCoefficientEstimatableParameterSettings` class
    for arc-wise treatment of the specified body's radiation pressure coefficient.







    )";



    } else if(name == "empirical_accelerations" && variant==0) {
        return R"(
        
Function for defining parameter settings for empirical acceleration magnitudes.

Factory function for creating a (linear sensitivity) parameter settings object for empirical acceleration magnitudes.
Using the empirical acceleration terms as estimatable parameters requires:

* The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms


Parameters
----------
body : str
    Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

centralBody : str
    Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

acceleration_components : dict[ EmpiricalAccelerationComponents, list[ EmpiricalAccelerationFunctionalShapes] ]
    Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
    for the specified body's empirical acceleration terms.







    )";



    } else if(name == "arcwise_empirical_accelerations" && variant==0) {
        return R"(
        
Function for defining parameter settings for arc-wise empirical acceleration magnitudes.

Factory function for creating a (linear sensitivity) parameter settings object for arc-wise empirical acceleration magnitudes (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.empirical_accelerations`).
Using the empirical acceleration terms as estimatable parameters requires:

* The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms

.. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.  


Parameters
----------
body : str
    Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

centralBody : str
    Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

acceleration_components : Dict[ EmpiricalAccelerationComponents, List[ EmpiricalAccelerationFunctionalShapes] ]
    Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.

arc_initial_times : List[ float ]
    List of times at which the arcs over which the empirical accelerations are to be estimated will start.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
    for the specified body's arc-wise empirical acceleration terms.







    )";



    } else if(name == "constant_empirical_acceleration_terms" && variant==0) {
        return R"(
        
Function for defining parameter settings for constant empirical acceleration terms.

As :func:`~tudatpy.numerical_simulation.estimation_setup.parameter.empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience


Parameters
----------
body : str
    Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

centralBody : str
    Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
    for the specified body's empirical acceleration terms.







    )";



    } else if(name == "arcwise_constant_empirical_acceleration_terms" && variant==0) {
        return R"(
        
Function for defining parameter settings for arc-wise constant empirical acceleration terms.

As :func:`~tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience


Parameters
----------
body : str
    Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

centralBody : str
    Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

arc_initial_times : List[ float ]
    List of times at which the arcs over which the empirical accelerations are to be estimated will start.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
    for the specified body's arc-wise constant empirical acceleration terms.







    )";



    } else if(name == "quasi_impulsive_shots" && variant==0) {
        return R"(
        
Function for defining parameter settings for quasi-impulsive shots.

Factory function for creating a (linear sensitivity) parameter settings object for so-called 'quasi-impulsive shots', such as desaturation maneuvers. With this parameter, the total :math:`\Delta V` vector of a set of such shots/maneuvers can be estimated.
Using the quasi-impulsive shots as an estimatable parameter requires:

* The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.quasi_impulsive_shots_acceleration` acceleration

.. note:: this parameter considers *all* shots/maneuvers used in the above acceleration model, and estimates the value of the 'delta_v_values' input of this acceleration.


Parameters
----------
body : str
    Name of the body, with which the quasi-impulsive shot estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's quasi-impulsive shots







    )";



    } else if(name == "gravitational_parameter" && variant==0) {
        return R"(
        
Function for defining parameter settings for a massive body's gravitational parameter.

Factory function for creating a (linear sensitivity) parameter settings object for the gravitational parameter of massive bodies.
Using the gravitational parameter as estimatable parameter requires:

* The body specified by the ``body`` parameter to be endowed with a gravity field (see :ref:`\`\`gravity_field\`\`` module for options)
* Any dynamical or observational model to depend on the gravitational parameter of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose gravitational model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's gravitational parameter.







    )";



    } else if(name == "spherical_harmonics_c_coefficients" && variant==0) {
        return R"(
        
Function for defining parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.

Factory function for creating a (linear sensitivity) parameter settings object for the spherical harmonics cosine-coefficients (:math:`\bar{C}_{lm}`) of a hody with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/0, and maximum degree/order 4/4, all spherical harmonic cosine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 0, 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
Using the spherical harmonics cosine coefficients as estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
* Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic` acceleration


Parameters
----------
body : str
    Name of the body, with whose gravitational model the estimatable parameters are associated.

minimum_degree : int
    Minimum degree of c-coefficients to be included.
minimum_order : int
    Minimum order of c-coefficients to be included.
maximum_degree : int
    Maximum degree of c-coefficients to be included.
maximum_order : int
    Maximum order of c-coefficients to be included.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
    for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.







    )";



    } else if(name == "spherical_harmonics_s_coefficients" && variant==0) {
        return R"(
        
Function for defining parameter settings for the sine coefficients of body's spherical harmonics gravitational model.

Factory function for creating a (linear sensitivity) parameter settings object for the spherical harmonics sine-coefficients (:math:`\bar{S}_{lm}`) of a hody with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/1 (there is no order 0 sine coefficient), and maximum degree/order 4/4, all spherical harmonic sine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
Using the spherical harmonics cosine coefficients as estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
* Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic` acceleration


Parameters
----------
body : str
    Name of the body, with whose gravitational model the estimatable parameters are associated.

minimum_degree : int
    Minimum degree of s-coefficients to be included.
minimum_order : int
    Minimum order of s-coefficients to be included.
maximum_degree : int
    Maximum degree of s-coefficients to be included.
maximum_order : int
    Maximum order of s-coefficients to be included.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
    for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.







    )";



    } else if(name == "spherical_harmonics_c_coefficients_block" && variant==0) {
        return R"(
        
Function for defining parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.

As :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_c_coefficients`, but with a manually defined set of coefficients.


Parameters
----------
body : str
    Name of the body, with whose gravitational model the estimatable parameters are associated.

block_indices : List[ Tuple[int, int] ]
    List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
    For each pair, the first value is the degree and the second the order of the coefficient to be included.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
    for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.







    )";



    } else if(name == "spherical_harmonics_s_coefficients_block" && variant==0) {
        return R"(
        
Function for defining parameter settings for the sine coefficients of body's spherical harmonics gravitational model.

As :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_s_coefficients`, but with a manually defined set of coefficients.


Parameters
----------
body : str
    Name of the body, with whose gravitational model the estimatable parameters are associated.

block_indices : List[ Tuple[int, int] ]
    List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
    For each pair, the first value is the degree and the second the order of the coefficient to be included.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
    for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.







    )";



    } else if(name == "constant_rotation_rate" && variant==0) {
        return R"(
        
Function for defining parameter settings for a body's constant rotation rate.

Factory function for creating a (linear sensitivity) parameter settings object for a body's constant rotation rate parameter.
Using the constant rotation rate as estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` or :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
* Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose rotation model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's constant spin rate.







    )";



    } else if(name == "rotation_pole_position" && variant==0) {
        return R"(
        
Function for defining parameter settings for a body's rotation pole position.

Factory function for creating a (linear sensitivity) parameter settings object for a body's rotation pole position, parameterized by the constant pole rotation angles (:math:`alpha` and :math:`\delta`).
Using the rotation pole position as estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` or :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
* Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose rotation model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's rotation pole position.







    )";



    } else if(name == "mean_moment_of_inertia" && variant==0) {
        return R"(
        
Function for defining parameter settings for a body's mean moment of inertia.

Factory function for creating a (linear sensitivity) parameter settings object for a body's mean moment of inertia. 
In most cases, the mean moment of inertia will not influence the dynamics/observation directly and sensitivity to this parameter will not be included. The dynamics/observation will be sensitive to this parameter if the rotational dynamics of a relevant body is estimated.
Using the mean moment of inertia as estimatable parameter requires:

* The estimation of an initial rotational state of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose body model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's mean moment of inertia.







    )";



    } else if(name == "periodic_spin_variations" && variant==0) {
        return R"(
        
Function for defining parameter settings for a body's periodic spin variations.

Factory function for creating a (linear sensitivity) parameter settings object for a body's periodic spin variation parameters.
Using the mean moment of inertia as estimatable parameter requires:

* A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
* Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose rotation model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's periodic spin variations.







    )";



    } else if(name == "polar_motion_amplitudes" && variant==0) {
        return R"(
        
Function for defining parameter settings for a body's polar motion amplitudes.

Factory function for creating a (linear sensitivity) parameter settings object for a body's polar motion amplitudes.
Using the polar motion amplitudes as estimatable parameter requires

* A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
* Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose rotation model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's polar motion amplitudes.







    )";



    } else if(name == "core_factor" && variant==0) {
        return R"(
        
Function for defining parameter settings for a body's core factor.

Factory function for creating a (linear sensitivity) parameter settings object for a body's core factor.
Using the core factor as estimatable parameter requires

* A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
* Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose rotation model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's core factor.







    )";



    } else if(name == "free_core_nutation_rate" && variant==0) {
        return R"(
        
Function for defining parameter settings for a body's free core nutation rate.

Factory function for creating a (linear sensitivity) parameter settings object for a body's free core nutation rate.
Using the free core nutation rate as estimatable parameter requires

* A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
* Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


Parameters
----------
body : str
    Name of the body, with whose rotation model the estimatable parameter is associated.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's free core nutation rate.







    )";



    } else if(name == "absolute_observation_bias" && variant==0) {
        return R"(
        
Function for defining parameter settings for an absolute observation bias.

Factory function for creating a (linear sensitivity) parameter settings object for an observation's absolute bias parameter.
Using the absolute observation bias as estimatable parameter requires:

* The observation model (corresponding to the `link_ends` and `observable_type`) to include an absolute bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias`)


Parameters
----------
link_ends : Dict[:class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType`, Tuple[str, str]
    Set of link ends that define the geometry of the biased observations.

observable_type : ObservableType
    Observable type of the biased observations.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ConstantObservationBiasEstimatableParameterSettings`
    for the specified observation's arc-wise absolute bias.







    )";



    } else if(name == "relative_observation_bias" && variant==0) {
        return R"(
        
Function for defining parameter settings for an relative observation bias.

Factory function for creating a (linear sensitivity) parameter settings object for an observation's relative bias parameter.
Using the relative observation bias as estimatable parameter requires

* The observation model (corresponding to the `link_ends` and `observable_type`) to include a relative bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias`)

.. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.  


Parameters
----------
link_ends : Dict[:class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType`, Tuple[str, str]
    Set of link ends that define the geometry of the biased observations.

observable_type : ObservableType
    Observable type of the biased observations.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ConstantObservationBiasEstimatableParameterSettings`
    for the specified observation's arc-wise relative bias.







    )";



    } else if(name == "arcwise_absolute_observation_bias" && variant==0) {
        return R"(
        
Function for defining parameter settings for arc-wise absolute observation bias.

Factory function for creating a (linear sensitivity) parameter settings object for the arc-wise treatment of an observation's absolute bias parameter.
Using the arc-wise absolute observation bias as estimatable parameter requires

* The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise absolute bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.arcwise_absolute_bias`)


Parameters
----------
link_ends : Dict[:class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType`, Tuple[str, str]
    Set of link ends that define the geometry of the biased observations.

observable_type : ObservableType
    Observable type of the biased observations.
arc_start_times : List[ float ]
    List of times at which the arcs over which the bias is to be estimated will start.
time_link_end : LinkEndType
    The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseConstantObservationBiasEstimatableParameterSettings`
    for the specified observation's arc-wise absolute bias.







    )";



    } else if(name == "arcwise_relative_observation_bias" && variant==0) {
        return R"(
        
Function for defining parameter settings for arc-wise absolute observation bias.

Factory function for creating a (linear sensitivity) parameter settings object for the arc-wise treatment of an observation's relative bias parameter.
Using the arc-wise relative observation bias as estimatable parameter requires

* The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise relative bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.arcwise_relative_bias`)

.. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.  


Parameters
----------
link_ends : Dict[:class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType`, Tuple[str, str]
    Set of link ends that define the geometry of the biased observations.

observable_type : ObservableType
    Observable type of the biased observations.
arc_start_times : List[ float ]
    List of times at which the arcs over which the bias is to be estimated will start.
time_link_end : LinkEndType
    The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseConstantObservationBiasEstimatableParameterSettings`
    for the specified observation's arc-wise relative bias.







    )";



    } else if(name == "ground_station_position" && variant==0) {
        return R"(
        
Function for defining parameter settings for ground station position bias.

Factory function for creating a (linear sensitivity) parameter settings object for a ground station's body-fixed Cartesian position.
Using the ground station position bias as estimatable parameter requires:

* At least one observation model to rely on the specified ground station


Parameters
----------
body : str
    Body name identifying the body, with which the ground station is associated.
ground_station_name : str
    Name which identifies the position-biased ground station.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified ground station's position bias.







    )";



    } else if(name == "ppn_parameter_gamma" && variant==0) {
        return R"(
        
Function for defining parameter settings for post-newtonian gamma parameter.

Factory function for creating a (linear sensitivity) parameter settings object for a global PPN :math:`\gamma` parameter.
Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:

* An acceleration model depending on this parameter, such as :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.relativistic_correction` 
* An observation model with a light-time correction depending on this parameter, such as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction`

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for a global post-newtonian :math:`\gamma` parameter.







    )";



    } else if(name == "ppn_parameter_beta" && variant==0) {
        return R"(
        
Function for defining parameter settings for post-newtonian beta parameter.

Factory function for creating a (linear sensitivity) parameter settings object for a global PPN :math:`\beta` parameter.
Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:
 
* An acceleration model depending on this parameter, such as :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.relativistic_correction` 
* An observation model with a light-time correction depending on this parameter (none yet implemented)

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings`
    :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for a global post-newtonian :math:`\beta` parameter.







    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace environment_setup {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "BodyListSettings") {
         return R"(

        Class for defining settings for the creation of a system of bodies.

        Class for defining settings for the creation of a system of bodies. This object is typically created from default settings, and
        then adapted to the user's specific needs.





     )";


    } else if(name == "BodyListSettings.frame_origin") {
         return R"(

        **read-only**

        Definition of the global frame origin for the bodies

        :type: str
     )";


    } else if(name == "BodyListSettings.frame_orientation") {
         return R"(

        **read-only**

        Definition of the global frame orientation for the bodies

        :type: str
     )";




    } else if(name == "BodyListSettings.get" && variant==0) {
            return R"(

        This function extracts a single BodySettings object .


        Parameters
        ----------
        body_name : str
            Name of the body for which settings are to be retrieved





    )";




    } else if(name == "BodySettings") {
         return R"(

        Class for defining settings for the creation of a single body.

        Class for defining settings for the creation of a single body, this object is typically stored inside a
        :class:`BodyListSettings`, object.





     )";


    } else if(name == "BodySettings.constant_mass") {
         return R"(

        Mass that gets assigned to the vehicle. This mass does *not* automatically define a gravity field
        model, but is instead used for the calculation of non-conservative forces only. When creating a body with a gravity field,
        leave this entry empty. NOTE: this option is a shorthand for assigning a mass-only
        :func:`~tudatpy.numerical_simulation.environment_setup.rigid_body.constant_rigid_body_properties` to ``mass_property_settings``, and will be deprecated.


        :type: float
     )";


    } else if(name == "BodySettings.atmosphere_settings") {
         return R"(

        Object that defines the settings of the atmosphere model that is to be created. Note that wind model settings
        may be defined inside this object. A variable of this type is typically assigned by using a factory function from the
        :ref:`\`\`atmosphere\`\`` module.


        :type: AtmosphereSettings
     )";


    } else if(name == "BodySettings.ephemeris_settings") {
         return R"(

        Object that defines the settings of the ephemeris model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`ephemeris\`\`` module.


        :type: EphemerisSettings
     )";


    } else if(name == "BodySettings.gravity_field_settings") {
         return R"(

        Object that defines the settings of the gravity field model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`gravity_field\`\`` module.


        :type: GravityFieldSettings
     )";


    } else if(name == "BodySettings.rotation_model_settings") {
         return R"(

        Object that defines the settings of the rotation model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`rotation_model\`\`` module.


        :type: RotationModelSettings
     )";


    } else if(name == "BodySettings.shape_settings") {
         return R"(

        Object that defines the settings of the shape model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`shape\`\`` module.


        :type: BodyShapeSettings
     )";


    } else if(name == "BodySettings.aerodynamic_coefficient_settings") {
         return R"(

        Object that defines the settings of the aerodynamic coefficient model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`aerodynamic_coefficients\`\`` module.


        :type: AerodynamicCoefficientSettings
     )";


    } else if(name == "BodySettings.gravity_field_variation_settings") {
         return R"(

        List of objects that define the settings of time variations of the gravity field variation models that are to be created. Variables in this list are typically
        assigned by using a factory function from the :ref:`\`\`gravity_field_variations\`\`` module.


        :type: list[GravityFieldVariationSettings]
     )";


    } else if(name == "BodySettings.shape_deformation_settings") {
         return R"(

        List of objects that define the settings of time variations of the exterior shape of natural bodies are to be created. Variables in this list are typically
        assigned by using a factory function from the :ref:`\`\`shape_deformation\`\`` module.


        :type: list[BodyDeformationSettings]
     )";


    } else if(name == "BodySettings.rigid_body_settings") {
         return R"(

        Object that defines the settings of the body rigid body (mass, center of mass, inertia) properties that are to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`rigid_body\`\`` module. Note that this setting does *not* define
        the gravity field, but rather only the mass, center of mass and inertia tensor.


        :type: RigidBodyPropertiesSettings
     )";


    } else if(name == "BodySettings.radiation_pressure_target_settings") {
         return R"(

        Object that defines the settings of the radiation pressure target model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`radiation_pressure\`\`` module. 


        :type: RadiationPressureTargetModelSettings
     )";


    } else if(name == "BodySettings.radiation_source_settings") {
         return R"(

        Object that defines the settings of the radiation source model that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`radiation_pressure\`\`` module. 


        :type: RadiationSourceModelSettings
     )";


    } else if(name == "BodySettings.vehicle_shape_settings") {
         return R"(

        Object that defines the settings of an exterior panelled vehicle shape that is to be created. A variable of this type is typically
        assigned by using a factory function from the :ref:`\`\`vehicle_systems\`\`` module. 


        :type: FullPanelledBodySettings
     )";






    } else if(name == "get_default_body_settings" && variant==0) {
        return R"(
        
Function that retrieves the default settings for the given set of input bodies.

Function that retrieves the default settings for the given set of input bodies. Default settings are described in
detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models.html>`_ .
Note that if a body is provided as input for which default settings do not exist, an exception is thrown. In addition
to settings for each separate body, this function returns an object that defines the global frame origin and orientation,


Parameters
----------
bodies : list[str]
    List of name of bodies for which default settings are to be retrieved and created.
base_frame_origin : str, default = 'SSB'
    Base frame origin of the set of bodies that is to be created. It defaults to the solar system barycenter (SSB), but it can by any of the bodies in `bodies_to_create` (provided it has an ephemeris defined).
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the set of bodies that is to be created. It can be either ECLIPJ2000 (default) or J2000.
Returns
-------
BodyListSettings
    Object containing the settings for the SystemOfBodies that are to be created






    )";



    } else if(name == "get_default_body_settings_time_limited" && variant==0) {
        return R"(
        
Function that retrieves the default settings for the given set of input bodies, with a limited valid time interval.

Same as :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but with body settings valid over a limited time interval. This makes the
the extraction of states from ephemerides more computationally efficient, at the expense of more RAM usage, and a
constrained time interval over which the ephemerides are valid. See `this page <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models/default_bodies_limited_time_range.html>`_ for more details.


Parameters
----------
bodies : list[str]
    List of name of bodies for which default settings are to be retrieved and created.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_origin : str
    Base frame origin of the set of bodies that is to be created.
base_frame_orientation : str
    Base frame orientation of the set of bodies that is to be created.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodyListSettings
    Object containing the settings for the SystemOfBodies that are to be created






    )";



    } else if(name == "get_default_single_body_settings" && variant==0) {
        return R"(
        
Function that retrieves the default settings for a single body.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but for retrieving default settings of only a single body


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved and created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
Returns
-------
BodySettings
    Object containing the settings for the bodt that is to be created






    )";



    } else if(name == "get_default_single_body_settings_time_limited" && variant==0) {
        return R"(
        
Function that retrieves the default settings for a single body, with a limited valid time interval.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved and created.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodySettings
    Object containing the settings for the bodt that is to be created






    )";


    } else if(name == "get_default_single_body_settings_time_limited" && variant==1) {
        return R"(
        
Function that retrieves the default settings for a single body, with a limited valid time interval.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body,
where the default settings of body with name ``source_body_name`` are retrieved and assigned to a body with name ``body_name``.
For instance, if ``source_body_name`` is set to "Mars", and ````body_name`` is set to "Earth" body name Earth will be created, with all the properties
of Mars


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved.
source_body_name : str
    Name of body for which default settings are to be retrieved, and assigned to a body with name ``body_name``.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodySettings
    Object containing the settings for the bodt that is to be created






    )";



    } else if(name == "get_default_single_alternate_body_settings" && variant==0) {
        return R"(
        
Function that retrieves the default settings for a single body, and assigns them to another body.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but for retrieving default settings of only a single body,
where the default settings of body with name ``source_body_name`` are retrieved and assigned to a body with name ``body_name``.
For instance, if ``source_body_name`` is set to "Mars", and ````body_name`` is set to "Earth" body name Earth will be created, with all the properties
of Mars


Parameters
----------
body_name : str
    Name of body for which default settings are to be created.
source_body_name : str
    Name of body for which default settings are to be retrieved, and assigned to a body with name ``body_name``.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
Returns
-------
BodySettings
    Object containing the settings for the bodt that is to be created






    )";



    } else if(name == "get_default_single_body_settings_time_limited" && variant==0) {
        return R"(
        
Function that retrieves the default settings for a single body, with a limited valid time interval.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved and created.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodySettings
    Object containing the settings for the bodt that is to be created






    )";


    } else if(name == "get_default_single_body_settings_time_limited" && variant==1) {
        return R"(
        
Function that retrieves the default settings for a single body, with a limited valid time interval.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body,
where the default settings of body with name ``source_body_name`` are retrieved and assigned to a body with name ``body_name``.
For instance, if ``source_body_name`` is set to "Mars", and ````body_name`` is set to "Earth" body name Earth will be created, with all the properties
of Mars


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved.
source_body_name : str
    Name of body for which default settings are to be retrieved, and assigned to a body with name ``body_name``.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodySettings
    Object containing the settings for the bodt that is to be created






    )";



    } else if(name == "add_aerodynamic_coefficient_interface" && variant==0) {
        return R"(
        
Function that creates an aerodynamic coefficient interface from settings, and adds it to an existing body.

This function can be used to add an aerodynamic coefficient interface to an existing body. It requires
settings for the aerodynamic coefficients, created using one of the factory functions from the `~tudatpy.numerical_simulation_environment_setup.aerodynamic_coefficient` module.
This function creates the actual coefficient interface from these settings, and assigns it to the
selected body. In addition to the identifier for the body to which it is assigned, this function
requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
inter-body dependencies in the coefficient interface


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the aerodynamic coefficients are to be assigned
coefficient_settings : AerodynamicCoefficientSettings
    Settings defining the coefficient interface that is to be created.





    )";



    } else if(name == "create_system_of_bodies" && variant==0) {
        return R"(
        
Function that creates a System of bodies from associated settings.

Function that creates a System of bodies from associated settings. This function creates the separate :class:`~tudatpy.numerical_simulation.Body`
objects and stores them in a :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` object. This object represents the full
physical environment in the simulation.


Parameters
----------
body_settings : BodyListSettings
    Object defining the physical environment, with all properties of artificial and natural bodies.
Returns
-------
:class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object containing the objects for bodies and environment models constituting the physical environment






    )";



    } else if(name == "create_simplified_system_of_bodies" && variant==0) {
        return R"(
        
Function that creates a simplified System of bodies.

Function that creates a simplified system of bodies. The following bodies are created in this system: the Sun, all planets of the Solar system, and Pluto. 
All bodies in this system use Gtop ephemerides and point mass gravity. The Earth is setup with a spherical shape model and a simple rotation model. 
The reference frame used to setup this simplified system of bodies has its origin at the SSB, and has an ECLIPJ2000 orientation.


Parameters
----------
initial_time : float, optional, default=0
    Initial system time in seconds since J2000.
Returns
-------
:class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object containing the objects for bodies and environment models constituting the physical environment






    )";



    } else if(name == "create_body_ephemeris" && variant==0) {
        return R"(
        
Function that creates an Ephemeris object.

Function that creates an :class:`~tudatpy.numerical_simulation.environment.Ephemeris` object, but does *not*
associate it with any specific body (e.g., it does not go into the environment, but can be used independently of it)


Parameters
----------
ephemeris_settings : EphemerisSettings
    Object defining the ephemeris settings.
body_name : str
    Name of body for which the ephemeris is created. Note that this input is only relevant for some ephemeris settings (for instance, a spice ephemeris setting), and it does *not* imply that the ephemeris object is associated with a Body object of this name.
Returns
-------
:class:`~tudatpy.numerical_simulation.environment.Ephemeris`
    Ephemeris object, created according to the provided settings






    )";



    } else if(name == "add_radiation_pressure_interface" && variant==0) {
        return R"(
        
Function that creates an radiation pressure interface from settings, and adds it to an existing body.

This function can be used to add an radiation pressure interface to an existing body. It requires
settings for the radiation pressure interface, created using one of the factory functions from the :ref:`\`\`radiation_pressure\`\`` module.
This function creates the actual coefficient interface from these settings, and assigns it to the
selected body. In addition to the identifier for the body to which it is assigned, this function
requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
inter-body dependencies in the radiation pressure interface


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the radiation pressure interface is to be assigned
radiation_pressure_settings : RadiationPressureInterfaceSettings
    Settings defining the radiation pressure interface that is to be created.





    )";



    } else if(name == "add_flight_conditions" && variant==0) {
        return R"(
        
Function that creates a flight conditions, and adds it to an existing body.

This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.FlightConditions` object to an existing body.
Typically, the ``FlightConditions`` are created automatically when they are required (for the calulcation of an
aerodynamic acceleration, or the saving of certain dependent variables). However, in some cases it may be useful
to manually trigger their creation, which is done through this function. If the ``central_body_name`` input
denotes a body that is endowed with an :class:`~tudatpy.numerical_simulation.environment.AtmosphereModel`, this function
automically creates an :class:`~tudatpy.numerical_simulation.environment.AtmosphericFlightConditions` object (capable of
calculating density, speed of sound, etc.), instead of the more basic :class:`~tudatpy.numerical_simulation.environment.FlightConditions`
(which is limited to properties such as altitude, latitude, etc.)


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body for which the flight conditions are to be created
central_body_name : str
    Name of the cenral body w.r.t. which the flight conditions are to be created (typically, but not necesarilly, the central body of propagation)/





    )";



    } else if(name == "add_rotation_model" && variant==0) {
        return R"(
        
Function that creates a rotation model, and adds it to an existing body.

This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.RotationalEphemeris` object to an existing body.
Typically, the ``RotationalEphemeris`` is created along with the `~tudatpy.numerical_simulation.environment.Body` itself However, in some cases it may be useful
to create a rotation model after the Body objects have been created. This function requires
settings for the rotation model, created using one of the factory functions from the :ref:`~tudatpy.numerical_simulation_environment_setup.rotation_model` module.
This function creates the actual coefficient interface from these settings, and assigns it to the
selected body. In addition to the identifier for the body to which it is assigned, this function
requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
inter-body dependencies in the radiation model


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the rotation model is to be assigned
rotation_model_settings
    Settings defining the rotation model that is to be created.





    )";



    } else if(name == "add_mass_properties_model" && variant==0) {
        return R"(
        
Function that creates a body mass property model, and adds it to an existing body.

This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.BodyMassProperties` object to an existing body.
Typically, the ``BodyMassProperties`` is created along with the `~tudatpy.numerical_simulation.environment.Body` itself However, in some cases it may be useful
to create body mass properties after the Body objects have been created. This function requires
settings for the body mass properties, created using one of the factory functions from the :ref:`~tudatpy.numerical_simulation_environment_setup.AAAA` module.
This function creates the actual body mass properties from these settings, and assigns it to the
selected body.


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the model is to be assigned
mass_property_settings
    Settings defining the mass properties model that is to be created.





    )";



    } else if(name == "add_engine_model" && variant==0) {
        return R"(
        
Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.

Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body. It creates and
object of class :class:`~tudatpy.numerical_simulation.environment.EngineModel`, and adds it to an existing body. Properties
assigned to this engine model are:
* The (constant) direction in body-fixed frame in which the engine is pointing (e.g. the body-fixed thrust direction when the engine is on)
* Settings for computing the thrust magnitude (as a function of time and/or other parameters), using a suitable function from the :ref:`\`\`thrust\`\`` submodule


Parameters
----------
body_name : str
    Name of the body to which the engine is to be added.
engine_name : str
    Name (e.g. unique identifier) of the engine that is to be added to the body
thrust_magnitude_settings : ThrustMagnitudeSettings
    Settings for computing the thrust magnitude (and specific impulse) as a function of time
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_fixed_thrust_direction : numpy.ndarray[numpy.float64[3, 1]], default = [1,0,0]
    Unit vector along which the thrust from the engine will point in a body-fixed frame





    )";



    } else if(name == "add_variable_direction_engine_model" && variant==0) {
        return R"(
        
Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.

Same as :func:`add_engine_model`, but with a time-variable body-fixed thrust direction


Parameters
----------
body_name : str
    Name of the body to which the engine is to be added.
engine_name : str
    Name (e.g. unique identifier) of the engine that is to be added to the body
thrust_magnitude_settings : ThrustMagnitudeSettings
    Settings for computing the thrust magnitude (and specific impulse) as a function of time
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_fixed_thrust_direction_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
    Function returning a unit vector, as a function of time, along which the thrust from the engine will point in a body-fixed frame





    )";



    } else {
        return "No documentation found.";
    }

}


    
namespace aerodynamic_coefficients {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "AerodynamicCoefficientSettings") {
         return R"(

        Base class for providing settings for aerodynamic interface model.

        Functional (base) class for settings of aerodynamic interface models that require no
        information in addition to their type.
        Aerodynamic interface model settings requiring additional information must be defined using an object derived from this class.





     )";


    } else if(name == "AerodynamicCoefficientSettings.add_force_contribution_to_moments") {
         return R"(

        Variable that toggles whether to add the force contribution to the moment coefficients as:

        .. math::
           \Delta \mathbf{C}_{M} = (\mathbf{r}_{ref}-\mathbf{r}_{com})\times \Delta \mathbf{C}_{F}
           
        where :math:`(\mathbf{r}_{ref}-\mathbf{r}_{com})` is the vector from the center of mass to the moment reference point, and :math:`\mathbf{C}_{F}` and :math:`\mathbf{C}_{M}` is the vector of forc and moment coefficients. Note that, if the force and moment coefficients afre defined in different frames, the relevant frame conversions are automatically performed.
        By default, his boolean is set to false, implicitly assuming that the moment coefficients are provided w.r.t. the (constant) center of mass.


        :type: bool
     )";




    } else if(name == "AerodynamicCoefficientSettings.add_single_control_surface" && variant==0) {
            return R"(

        Function to add settings for a single control surface to the coefficient settings. Note that, in Tudat, the
        control surface aerodynamic database inherits the reference properties (length, area, moment reference point)
        from the ``AerodynamicCoefficientSettings`` to which it is assigned.



        Parameters
        ----------
        control_surface_settings : ControlSurfaceIncrementAerodynamicCoefficientSettings
            Settings for aerodynamic coefficients of a control surface 

        control_surface_name : str
            Name by which the control surface will be identified





    )";




    } else if(name == "ConstantAerodynamicCoefficientSettings") {
         return R"(

        Class for defining model settings from constant aerodynamic coefficients.

        `AerodynamicCoefficientSettings` derived class for aerodynamic interface model settings using only constant aerodynamic coefficients.




     )";






    } else if(name == "constant" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings entirely from constant coefficients.

Factory function for settings object, defining aerodynamic interface model entirely from constant aerodynamic coefficients,
i.e. coefficients are not a function of any independent variables.


Parameters
----------
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
constant_force_coefficient : numpy.ndarray
    Constant force coefficients.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift

Returns
-------
ConstantAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ConstantAerodynamicCoefficientSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for the artificial body "Vehicle", using only constant aerodynamic coefficients:

.. code-block:: python 
  
  # Define the reference area and constant aerodynamic coefficients 
  reference_area = 20.0 
  drag_coefficient = 1.5 
  lift_coefficient = 0.3 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant( 
      reference_area, 
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient], 
      force_coefficients_frame=environment.negative_aerodynamic_frame_coefficients, 
  ) 
  # Assign aerodynamic interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings) 


    )";



    } else if(name == "custom_aerodynamic_force_coefficients" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings from custom coefficients.

Factory function for settings object, defining aerodynamic interface model via a custom force coefficient function
(function of independent variable).


Parameters
----------
force_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[3, 1]]]
    Function that is defining the aerodynamic coefficients as function of an independent variable (see arg independent_variable_names).
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_name : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

Returns
-------
CustomAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.CustomAerodynamicCoefficientSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for the artificial body "Vehicle", using a function based on the mach number: 

.. code-block:: python 
  
  def force_coefficients(variables_list): 
    # Extract the mach number 
    mach_number = variables_list[0] 
    # If the mach number is below 3, use fixed coefficients 
    if mach_number <= 3: 
        return [0.99, 0, 1.08] 
    # Same if the mach number is above 10 
    elif mach_number >= 10: 
        return [0.82, 0, 0.88] 
    # Otherwise, vary linearly between the ones at M=3 and M=10 
    CD = 1.0667-0.02457*mach_number 
    CL = 1.1636-0.02786*mach_number 
    return [CD, 0, CL] 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.custom( 
      force_coefficients, 
      reference_area=1.50, 
      independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent] 
  ) 
  # Assign the aerodynamic coefficient interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings) 


    )";



    } else if(name == "custom_aerodynamic_force_and_moment_coefficients" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings from custom coefficients.

Factory function for settings object, defining aerodynamic interface model via a custom force and moment coefficient function
(function of independent variable).


Parameters
----------
force_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[3, 1]]]
    Function that is defining the aerodynamic force coefficients as function of an independent variable (see arg independent_variable_names).
moment_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[3, 1]]]
    Function that is defining the aerodynamic moment coefficients as function of an independent variable (see arg independent_variable_names).
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
reference_length : float
    Reference length with which aerodynamic moments are non-dimensionalized.
independent_variable_name : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

moment_coefficients_frame : AerodynamicCoefficientFrames, default = positive_body_frame_coefficients
    Variable defining the frame in which the moment coefficients are defined. By default, this is the positive body
    frame, so that the coefficients are roll, pitch and yaw (:math:`C_{l}, C_{m}, C_{n}`)

moment_reference_point : numpy.ndarray[numpy.float64[3, 1]] = np.full([3, 1], np.nan)
    Point w.r.t. aerodynamic moment coefficients are defined. This variable is used to calculate the contribution of the aerodynamic
    force coefficients to the effective moment coefficients. See the ``add_force_contribution_to_moments`` attribute of the 
    :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for more details.
    If the present input is set to NaN (as is the default), the reference point is left undefined, and the aerodynamic moments are computed
    without computing any force coefficient contribution to the moment coefficients.

Returns
-------
CustomAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.CustomAerodynamicCoefficientSettings` class






    )";



    } else if(name == "tabulated" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings from user-defined, 1-d tabulated coefficients.

Factory function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force and moment coefficients
(tabulated w.r.t. independent variable).


Parameters
----------
independent_variables : list[float]
    Values of independent variables at which the coefficients in the input multi vector are defined (size 1).
force_coefficients : list[numpy.ndarray[numpy.float64[3, 1]]]
    Values of force coefficients at independent variables defined by independent_variables.
moment_coefficients : list[numpy.ndarray[numpy.float64[3, 1]]]
    Values of moment coefficients at independent variables defined by independent_variables.
reference_length : float
    Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_name : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

moment_coefficients_frame : AerodynamicCoefficientFrames, default = positive_body_frame_coefficients
    Variable defining the frame in which the moment coefficients are defined. By default, this is the positive body
    frame, so that the coefficients are roll, pitch yaw (:math:`C_{l}, C_{m}, C_{n}`)

moment_reference_point : numpy.ndarray[numpy.float64[3, 1]] = np.full([3, 1], np.nan)
    Point w.r.t. aerodynamic moment coefficients are defined. This variable is used to calculate the contribution of the aerodynamic
    force coefficients to the effective moment coefficients. See the ``add_force_contribution_to_moments`` attribute of the 
    :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for more details.
    If the present input is set to NaN (as is the default), the reference point is left undefined, and the aerodynamic moments are computed
    without computing any force coefficient contribution to the moment coefficients.

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class (via :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettingsBase` class)





Examples
--------
In this example, aerodynamic force and moment coefficients are defined as multi-dimensional arrays. 
The values for the aerodynamic coefficients vary with Mach number, and are defined for Mach numbers of 3, 5, 10, and 15. 
This example also shows how to set the required reference point, lengths, and area. 

.. code-block:: python 
  
  # Define the aerodynamic force coefficients [CD, CS, CL] for different mach numbers 
  aero_coefficients_array_force = [ 
      [0.7647, 0, 0.9722], 
      [0.6729, 0, 0.8461], 
      [0.6240, 0, 0.7838], 
      [0.6246, 0, 0.7841] 
  ] 
  # Define the aerodynamic moment coefficients for different mach numbers 
  aero_coefficients_array_moment = [ 
      [0.45, 0, 0], 
      [0.50, 0, 0], 
      [0.53, 0, 0], 
      [0.55, 0, 0] 
  ] 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated( 
      independent_variables=[3, 5, 10, 15],       # Mach number at which the coefficients are defined 
      force_coefficients=aero_coefficients_array_force, 
      moment_coefficients=aero_coefficients_array_moment, 
      reference_length=0.25, 
      reference_area=1.50, 
      independent_variable_name=environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent 
  ) 
  # Assign the aerodynamic coefficient interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings) 


    )";



    } else if(name == "tabulated_force_only" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings from user-defined, 1-d tabulated force coefficients.

Factory function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force coefficients
(tabulated w.r.t. independent variable).


Parameters
----------
independent_variables : list[float]
    Values of independent variables at which the coefficients in the input multi vector are defined (size 1)
force_coefficients : list[numpy.ndarray[numpy.float64[3, 1]]]
    Values of force coefficients at independent variables defined by independent_variables.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_name : environment.AerodynamicCoefficientsIndependentVariables
    Identifier of the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
    Pointer to an interpolator settings object where the conditions for interpolation of tabulated inputs are saved.

Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class





Examples
--------
In this example, aerodynamic force coefficients are defined as a multi-dimensional array. 
The values for the force coefficients vary with Mach number, and are defined for Mach numbers of 3, 5, 10, and 15. 

.. code-block:: python 
  
  # Define the aerodynamic coefficients [CD, CS, CL] for different mach numbers 
  aero_coefficients_array = [ 
      [0.7647, 0, 0.9722], 
      [0.6729, 0, 0.8461], 
      [0.6240, 0, 0.7838], 
      [0.6246, 0, 0.7841] 
  ] 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_force_only( 
      independent_variables=[3.0, 5.0, 10.0, 15.0],       # Mach number at which the coefficients are defined 
      force_coefficients=aero_coefficients_array, 
      reference_area=1.50, 
      independent_variable_name=environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent 
  ) 
  # Assign the aerodynamic coefficient interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings) 


    )";



    } else if(name == "tabulated_force_only_from_files" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings from tabulated force coefficients from files.

Factory function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force coefficients
(tabulated w.r.t. independent variable), obtained from data files.


Parameters
----------
force_coefficient_files : Dict[int, str]
    Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_names : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.

Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class





Examples
--------
In this example, the drag and lift coefficients of the Space Transport System are defined from two data files. 
Both of these data files contain coefficient values dependent on both the angle of attack and the mach number, 
as shown in the example in the `independent_variable_names` input. 
This example is taken from the `reentry trajectory example <https://github.com/tudat-team/tudatpy-examples/blob/1f8180b0064226175bbe66e3eaf044f229a897f6/propagation/reentry_trajectory.py>`_. 

.. code-block:: python 
  
  # Define the aerodynamic coefficient files (leave C_S empty) 
  aero_coefficients_files = {0: "input/STS_CD.dat", 2:"input/STS_CL.dat"} 
  # Setup the aerodynamic coefficients settings tabulated from the files 
  coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_force_only_from_files( 
      force_coefficient_files=aero_coefficients_files, 
      reference_area=2690.0*0.3048*0.3048, 
      independent_variable_names=[environment.angle_of_attack_dependent, environment.mach_number_dependent] 
  ) 
  # Add predefined aerodynamic coefficients database to the body 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "STS", coefficient_settings) 


    )";



    } else if(name == "tabulated_from_files" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings from tabulated coefficients from files.

Factory function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force and moment coefficients
(tabulated w.r.t. independent variable), obtained from data files.


Parameters
----------
force_coefficient_files : Dict[int, str]
    Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key (0, 1 and 2 a are x-, y- and z-axis of force frame, respectively).
moment_coefficient_files : Dict[int, str]
    Path of the aerodynamic coefficient files corresponding to the moment coefficient of the given dict key (0, 1 and 2 a are x-, y- and z-axis of moment frame, respectively).
reference_length : float
    Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_names : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

moment_coefficients_frame : AerodynamicCoefficientFrames, default = positive_body_frame_coefficients
    Variable defining the frame in which the moment coefficients are defined. By default, this is the positive body
    frame, so that the coefficients are roll, pitch yaw (:math:`C_{l}, C_{m}, C_{n}`)

moment_reference_point : numpy.ndarray[numpy.float64[3, 1]] = np.full([3, 1], np.nan)
    Point w.r.t. aerodynamic moment coefficients are defined. This variable is used to calculate the contribution of the aerodynamic
    force coefficients to the effective moment coefficients. See the ``add_force_contribution_to_moments`` attribute of the 
    :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for more details.
    If the present input is set to NaN (as is the default), the reference point is left undefined, and the aerodynamic moments are computed
    without computing any force coefficient contribution to the moment coefficients.

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved. 

Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class





Examples
--------
This example is very similar to the one for `tabulated_force_only_from_files`, with the distinction that a pitching moment coefficient is added.

.. code-block:: python 
  
  # Define the force coefficient files (leave C_S empty) 
  force_coefficients_files = {0: "input/STS_CD.dat", 2:"input/STS_CL.dat"} 
  # Define the moment coefficient files (leave C_S empty) 
  moment_coefficients_files = {0: "input/STS_CM.dat"} 
  # Setup the aerodynamic coefficients settings tabulated from the files 
  coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_from_files( 
      force_coefficient_files=force_coefficients_files, 
      moment_coefficient_files=moment_coefficients_files, 
      reference_length=11.9, 
      reference_area=2690.0*0.3048*0.3048, 
      independent_variable_names=[environment.angle_of_attack_dependent, environment.mach_number_dependent] 
  ) 
  # Add the predefined aerodynamic coefficients database to the body 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "STS", coefficient_settings) 


    )";



    } else if(name == "scaled_by_constant" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings by applying one constant scaling factor/value to all coefficients of an existing model settings object.

Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by one constant factor or value.
Via the ``is_scaling_absolute`` 
boolean, the user can apply a constant scaling factor or an absolute value to the resulting force and moment coefficients (for instance for an uncertainty analysis).


Parameters
----------
unscaled_coefficient_settings : AerodynamicCoefficientSettings
    Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
force_scaling_constant : float
    Constant scaling factor to be applied to all aerodynamic force coefficients.
moment_scaling_constant : float
    Constant scaling factor to be applied to all aerodynamic moment coefficients.
is_scaling_absolute : bool, default = False
    Boolean indicating whether aerodynamic coefficient scaling is absolute.
    Setting this boolean to true will add the scaling value to the base value,
    instead of the default behaviour of multiplying the base value by the scaling factor.

Returns
-------
ScaledAerodynamicCoefficientInterfaceSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class





Examples
--------
In this example, we first set constant aerodynamic coefficients, like in the earlier example. 
Then, we use the `scaled_by_constant` function to scale the force coefficients by 1.1. 
Since the `is_scaling_absolute` equals `False` by default, the force coefficients are then increased by 10%. 

.. code-block:: python 
  
  # Define the reference area and constant aerodynamic coefficients 
  reference_area = 20.0 
  drag_coefficient = 1.5 
  lift_coefficient = 0.3 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant( 
      reference_area, 
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient] 
  ) 
  # Define scaled aerodynamic coefficient to increase the force coefficients by 10% 
  scaled_aero_coefficient_settings = environment_setup.aerodynamic_coefficients.scaled_by_constant( 
      unscaled_coefficient_settings=aero_coefficient_settings, 
      force_scaling_constant=1.1, 
      moment_scaling_constant=1.0 
  ) 
  # Assign aerodynamic interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", scaled_aero_coefficient_settings) 


    )";



    } else if(name == "scaled_by_vector" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings by applying constant scaling factors/values to the coefficients of an existing model settings object.

Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by constant factors or values.
Via the ``is_scaling_absolute`` boolean, the user can apply one constant scaling factor or an absolute value to each resulting force and moment coefficient (for instance for an uncertainty analysis).


Parameters
----------
unscaled_coefficient_settings : AerodynamicCoefficientSettings
    Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
force_scaling_vector : numpy.ndarray[numpy.float64[3, 1]]
    Constant scaling factors to be applied to each aerodynamic force coefficient.
moment_scaling_vector : numpy.ndarray[numpy.float64[3, 1]]
    Constant scaling factors to be applied to each aerodynamic moment coefficient.
is_scaling_absolute : bool, default = False
    Boolean indicating whether aerodynamic coefficient scaling is absolute.
    Setting this boolean to true will add the scaling value to the base value,
    instead of the default behaviour of multiplying the base value by the scaling factor.

Returns
-------
ScaledAerodynamicCoefficientInterfaceSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class





Examples
--------
In this example, we first set constant aerodynamic coefficients, like in the earlier example. 
Then, we use the `scaled_by_vector` function to scale the drag coefficient by 2. 

.. code-block:: python 
  
  # Define the reference area and constant aerodynamic coefficients 
  reference_area = 20.0 
  drag_coefficient = 1.5 
  lift_coefficient = 0.3 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant( 
      reference_area, 
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient] 
  ) 
  # Define scaled aerodynamic coefficient to increase CD by a factor of 2 
  scaled_aero_coefficient_settings = environment_setup.aerodynamic_coefficients.scaled_by_vector( 
      unscaled_coefficient_settings=aero_coefficient_settings, 
      force_scaling_vector=[2.0, 1.0, 1.0], 
      moment_scaling_vector=[1.0, 1.0, 1.0] 
  ) 
  # Assign aerodynamic interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", scaled_aero_coefficient_settings) 


    )";



    } else if(name == "scaled_by_vector_function" && variant==0) {
        return R"(
        
Factory function for creating aerodynamic interface model settings by applying custom scaling factors/values to the coefficients of an existing model settings object.

Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by custom factors or values.
Via the ``is_scaling_absolute`` boolean, the user can apply the scaling factors or absolute values to each resulting force and moment coefficient (for instance for an uncertainty analysis).


Parameters
----------
unscaled_coefficient_settings : AerodynamicCoefficientSettings
    Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
force_scaling_vector_function : callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
    Custom scaling factors to be applied to each aerodynamic force coefficient.
moment_scaling_vector_function : callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
    Custom scaling factors to be applied to each aerodynamic moment coefficient.
is_scaling_absolute : bool, default = False
    Boolean indicating whether aerodynamic coefficient scaling is absolute.
    Setting this boolean to true will add the scaling value to the base value,
    instead of the default behaviour of multiplying the base value by the scaling factor.

Returns
-------
ScaledAerodynamicCoefficientInterfaceSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class





Examples
--------
In this example, we first set constant aerodynamic coefficients, like in the earlier example. 
Then, we use the `scaled_by_vector_function` function to scale the drag and lift coefficients according to a function that varies with time. 
This scaling function essentially adds noise to the CD and CL following as a sin or cos function. 

.. code-block:: python 
  
  # Define the reference area and constant aerodynamic coefficients 
  reference_area = 20.0 
  drag_coefficient = 1.5 
  lift_coefficient = 0.3 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant( 
      reference_area, 
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient] 
  ) 
  # Define the aerodynamic coefficient scaling as a function of time 
  def aero_coefficient_scaling(time): 
      CD_scale = 1 + 0.25*np.sin(time/10) 
      CL_scale = 1 + 0.25*np.cos(time/15) 
      return [CD_scale, 1.0, CL_scale] 
  # Define scaled aerodynamic coefficient to increase CD by a factor of 2 
  scaled_aero_coefficient_settings = environment_setup.aerodynamic_coefficients.scaled_by_vector_function( 
      unscaled_coefficient_settings=aero_coefficient_settings, 
      force_scaling_vector_function=aero_coefficient_scaling, 
      moment_scaling_vector_function=lambda x: [1.0, 1.0, 1.0] 
  ) 
  # Assign aerodynamic interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", scaled_aero_coefficient_settings) 


    )";



    } else if(name == "custom_control_surface" && variant==0) {
        return R"(
        
Factory function for creating control surface aerodynamic model settings from custom coefficients.

Factory function for settings object, defining control surface aerodynamic interface model via a custom force and moment coefficient function
(function of independent variable). This function is essentically the control-surface equivalent of the
:func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_and_moment_coefficients` function for body coefficient settings.


Parameters
----------
force_and_moment_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[6, 1]]]
    Function that is defining the aerodynamic force (first three entries) and moment (last three entries) coefficients as function of an independent variables (see  ``independent_variable_names``).
independent_variable_names : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the control surface aerodynamic coefficients are defined. Typically, one entry from this list will be ``control_surface_deflection_dependent``
Returns
-------
ControlSurfaceIncrementAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ControlSurfaceIncrementAerodynamicCoefficientSettings` derived class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ControlSurfaceIncrementAerodynamicCoefficientSettings` for the artificial body "Vehicle", using a function based on the mach number: 

.. code-block:: python 
  
  def force_and_moment_coefficients(variables_list): 
    # Extract the mach number 
    mach_number = variables_list[0] 
    # If the mach number is below 3, use fixed coefficients 
    if mach_number <= 3: 
        return [0.99, 0, 1.08, 0.94, 0.35, 0, 0] 
    # Same if the mach number is above 10 
    elif mach_number >= 10: 
        return [0.82, 0, 0.88, 0.55, 0, 0] 
    # Otherwise, vary linearly between the ones at M=3 and M=10 
    CD = 1.0667-0.02457*mach_number 
    CL = 1.1636-0.02786*mach_number 
    Cl = 0.35 + 0.02857*mach_number 
    return [CD, 0, CL, Cl, 0, 0] 
  # Create the aerodynamic interface settings 
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.custom_control_surface( 
      force_and_moment_coefficients, 
      independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent] 
  ) 
  # Assign the aerodynamic coefficient interface to the vehicle 
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings) 
    


    )";



    } else if(name == "tabulated_from_files_control_surface" && variant==0) {
        return R"(
        
Factory function for creating control surface aerodynamic model settings from tabulated coefficients from files.

Factory function for settings object, defining control surface aerodynamic interface model via user-defined, tabulated aerodynamic force and moment coefficients
(tabulated w.r.t. independent variable), obtained from data files.. This function is essentically the control-surface equivalent of the
:func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_from_files` function for body coefficient settings.

Returns
-------
ControlSurfaceIncrementAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ControlSurfaceIncrementAerodynamicCoefficientSettings` derived class






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace atmosphere {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "WindModelSettings") {
         return R"(

        Class for providing settings for wind model.

        Functional (base) class for settings of wind models that require no information in addition to their type.
        Wind model classes requiring additional information must be created using an object derived from this class.





     )";





    } else if(name == "AtmosphereSettings") {
         return R"(

        Base class for providing settings for atmosphere model.

        Functional (base) class for settings of atmosphere models that require no information in addition to their type.
        Atmosphere model classes requiring additional information must be created using an object derived from this class.





     )";


    } else if(name == "AtmosphereSettings.wind_settings") {
         return R"(

        **read-only**

        Wind model settings for the atmosphere model settings object.

        :type: WindModelSettings
     )";





    } else if(name == "ExponentialAtmosphereSettings") {
         return R"(

        Class for providing settings for exponential atmosphere model.

        `AtmosphereSettings` derived class for a defining the settings of an exponential atmosphere model.




     )";






    } else if(name == "constant_wind_model" && variant==0) {
        return R"(
        
Factory function for creating wind model settings with constant wind velocity.

Factory function for settings object, defining wind model entirely from constant wind velocity in a given reference frame.


Parameters
----------
wind_velocity : numpy.ndarray[numpy.float64[3, 1]]
    Constant wind velocity in the specified reference frame.

associated_reference_frame : numerical_simulation.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
    Reference frame in which constant wind velocity is defined.

Returns
-------
ConstantWindModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ConstantWindModelSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings`,
using a constant wind-velocity vector defined in a vertical aerodynamic reference frame:

.. code-block:: python 
  
  # Define the wind in 3 directions in the vertical reference frame 
  wind_Xv = 3     # Meridional wind of +3 m/s (pointing to the North) 
  wind_Yv = 5     # Zonal wind of +5 m/s (pointing to the West) 
  wind_Zv = -11   # Vertical wind of +11 m/s (pointing out of the centre of the Earth) 
  # Create the constant wind settings 
  constant_wind = environment_setup.atmosphere.constant_wind_model( 
    [wind_Xv, wind_Yv, wind_Zv], 
    environment.AerodynamicsReferenceFrames.vertical_frame) 
  # Apply the constant wind settings to the Earth atmosphere settings 
  body_settings.get("Earth").atmosphere_settings.wind_settings = constant_wind 


    )";



    } else if(name == "custom_wind_model" && variant==0) {
        return R"(
        
Factory function for creating wind model settings with custom wind velocity.

Factory function for settings object, defining wind model entirely from custom wind velocity function in a given reference frame.
The custom wind velocity has to be given as a function of altitude, longitude, latitude and time.

.. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
          The altitude is in meters, and the time is a Julian date in seconds since J2000.


Parameters
----------
wind_velocity : callable[[float, float, float, float], numpy.ndarray[numpy.float64[3, 1]]]
    Custom wind velocity function (w.r.t. altitude, longitude, latitude and time) in the specified reference frame.

associated_reference_frame : numerical_simulation.environment.AerodynamicsReferenceFrames, default = AerodynamicsReferenceFrames.vertical_frame
    Reference frame in which wind velocity is defined.

Returns
-------
CustomWindModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomWindModelSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.WindModelSettings`,
using a user-defined wind-velocity function (of altitude, longitude, latitude and time), defined in a vertical aerodynamic reference frame:

.. code-block:: python 
  
  # Define the wind in 3 directions in the vertical reference frame 
  def wind_function(h, lon, lat, time): 
      # Meridional wind (pointing North) depends on latitude [deg] and time [sec since J2000] 
      wind_Xv = lat*10/time 
      # Zonal wind (pointing West) only depends on the longitude [deg] 
      wind_Yv = 5/lon 
      # Vertical wind (pointing out of the centre of the Earth) only depends on the altitude [m] 
      wind_Zv = 1000/h 
      # Return the custom wind 
      return [wind_Xv, wind_Yv, wind_Zv] 
  # Create the custom wind settings 
  custom_wind = environment_setup.atmosphere.custom_wind_model( 
      wind_function, 
      environment.AerodynamicsReferenceFrames.vertical_frame) 
  # Apply the custom wind settings to the Earth atmosphere settings 
  body_settings.get("Earth").atmosphere_settings.wind_settings = custom_wind 


    )";



    } else if(name == "exponential_predefined" && variant==0) {
        return R"(
        
Factory function for creating atmospheric model settings from pre-defined exponential model.

Factory function for settings object, defining atmosphere model from pre-defined exponential model.
The pre-encoded properties are available for Earth and Mars, as can be seen on the table below.
This function creates an instance of an `AtmosphereSettings` derived `ExponentialAtmosphereSettings` object.

.. list-table:: Pre-defined exponential atmosphere model properties
  :widths: 25 25 25 25
  :header-rows: 1

  * - Property
    - Earth
    - Mars
    - Units
  * - Scale Height
    - 7.2
    - 11.1
    - km
  * - Density at Zero Altitude
    - 1.225
    - 0.02
    - kg/m :math:`{}^3`
  * - Constant Temperature
    - 246.0
    - 215.0
    - K
  * - Specific Gas Constant
    - 287.0
    - 197.0
    - J/kg/K
  * - Ratio of Specific Heats
    - 1.4
    - 1.3
    - --


Parameters
----------
body_name : str
    Body for which pre-defined model settings are to be loaded. Available bodies "Earth", "Mars".

Returns
-------
ExponentialAtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Mars, 
using the interface of the predefined exponential model, using pre-encoded values: 

.. code-block:: python 
   
   # Create atmosphere settings and add to body settings of "Mars"  
   body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential_predefined("Mars")  


    )";



    } else if(name == "exponential" && variant==0) {
        return R"(
        
Factory function for creating atmospheric model settings from fully parametrized exponential model.

Factory function for settings object, defining exponential atmosphere model.
The model is solely based on an exponentially decaying density profile with a constant temperature and composition
(i.e. independent of time, latitude and longitude).

The user has access to a fully parametrized model, meaning that in addition to the required input parameters ``scale_height`` and ``surface_density`` (ground-level air density),
the user can specify non-standard values for constant temperature, gas constant and specific heats ratio.


Parameters
----------
scale_height : float
    Scale height for density profile of atmosphere.
surface_density : float
    Atmospheric density at ground level.
constant_temperature : float, default = 288.15
    Constant atmospheric temperature.
specific_gas_constant : float, default = constants.SPECIFIC_GAS_CONSTANT_AIR
    Specific gas constant for (constant) atmospheric chemical composition.
ratio_specific_heats : float, default = 1.4
    Ratio of specific heats for (constant) atmospheric chemical composition.
Returns
-------
ExponentialAtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ExponentialAtmosphereSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
using the minimalist interface to the exponential model and taking parameters with classic values for Earth:

.. code-block:: python 
   
   # define parameters of an invariant exponential atmosphere model  
   density_scale_height = 7.2E3  
   constant_temperature = 290  
   # create atmosphere settings and add to body settings of "Earth"  
   body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.exponential(  
   	density_scale_height, density_at_zero_altitude)  


    )";



    } else if(name == "nrlmsise00" && variant==0) {
        return R"(
        
Factory function for creating NRLMSISE-00 atmospheric model settings.

Factory function for settings object, defining atmosphere model in accordance to the NRLMSISE-00 global reference model for Earth's atmosphere.


Parameters
----------
space_weather_file : str, default = :func:`~tudatpy.io.get_space_weather_path` + 'sw19571001.txt'
    File to be used for space weather characteristics as a function of time (e.g. F10.7, Kp, etc.). The file is typically taken from here `celestrak <https://celestrak.org/SpaceData/sw19571001.txt>`_ (note that the file in your resources path will not be the latest version of this file; download and replace your existing file if required). Documentation on the file is given `here <https://celestrak.org/SpaceData/SpaceWx-format.php>`_
Returns
-------
AtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
using the NRLMSISE-00 global reference model:

.. code-block:: python 
   
   # create atmosphere settings and add to body settings of body "Earth"
   body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.nrlmsise00() 


    )";



    } else if(name == "us76" && variant==0) {
        return R"(
        
Factory function for creating US76 standard atmosphere model settings.

Factory function for creating US76 standard atmosphere model settings. The model is defined using tabulated data for density, pressure and temperature,
from an altitude of -5 km up to 1000 km. Up to 100 km, a data point is provided every 100 m. Above 100 km, a data point is provided every 1 km. The data
are interpolated using a cubic spline interpolator. Note that this model is specific to Earth's atmosphere.

Returns
-------
AtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
using the US76 standard atmosphere model:

.. code-block:: python 
   
   # create atmosphere settings and add to body settings of body "Earth"
   body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.us76() 

   


    )";



    } else if(name == "custom_constant_temperature" && variant==0) {
        return R"(
        
Factory function for creating atmospheric model settings from custom density profile.

Factory function for settings object, defining constant temperature atmosphere model from custom density profile.
The user is specifying the density profile as a function of altitude.
The value of pressure is computed by assuming hydrostatic equilibrium, temperature, gas constant and the ratio of specific heats are modelled as constants.


Parameters
----------
density_function : callable[[float], float]
    Function to retrieve the density at the current altitude.

constant_temperature : float
    Constant atmospheric temperature.
specific_gas_constant : float, default = 287.0
    Specific gas constant for (constant) atmospheric chemical composition.
ratio_specific_heats : float, default = 1.4
    Ratio of specific heats for (constant) atmospheric chemical composition.
Returns
-------
CustomConstantTemperatureAtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
with constant temperature and composition, but a density which varies with altitude according to a user-defined model:

.. code-block:: python 
  
  # Define the density as a function of altitude [in m] 
  def density_function(h): 
      # Return the density according to a modified exponential model 
      return 1.15 * np.exp(-h/7300) 
  # Define parameters for constant temperature and composition 
  constant_temperature = 250.0 
  specific_gas_constant = 300.0 
  ratio_of_specific_heats = 1.4 
  # Create the custom constant temperature atmosphere settings 
  custom_density_settings = environment_setup.atmosphere.custom_constant_temperature( 
      density_function, 
      constant_temperature, 
      specific_gas_constant, 
      ratio_of_specific_heats) 
  # Add the custom density to the body settings of "Earth" 
  body_settings.get("Earth").atmosphere_settings = custom_density_settings 


    )";



    } else if(name == "custom_four_dimensional_constant_temperature" && variant==0) {
        return R"(
        
Factory function for creating atmospheric model settings from custom density profile.

Factory function for settings object, defining constant temperature atmosphere model from custom density profile.
The user is specifying the density profile as a function of altitude, longitude, latitude and time.

.. note:: The longitude and latitude will be passed to the function in **degree** and not in radians.
          The altitude is in meters, and the time is a Julian date in seconds since J2000.


Parameters
----------
density_function : callable[[float, float, float, float], float]
    Function to retrieve the density at the current altitude, longitude, latitude and time.

constant_temperature : float
    Constant atmospheric temperature.
specific_gas_constant : float, default = 287.0
    Specific gas constant for (constant) atmospheric chemical composition.
ratio_specific_heats : float, default = 1.4
    Ratio of specific heats for (constant) atmospheric chemical composition.
Returns
-------
CustomConstantTemperatureAtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.CustomConstantTemperatureAtmosphereSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
with constant temperature and composition (gas constant and ratio of specific heats), but a density which varies with altitude, longitude, latitude and time, according to a user-defined model:

.. code-block:: python 
  
  # Define the density as a function of altitude [m], longitude [deg], latitude [deg], and time [sec since J2000] 
  def density_function(h, lon, lat, time): 
      # Return the density according to an exponential model that varies with time to add noise with a sine (ignore lon/lat) 
      return (1 + 0.15 * np.sin(time/10)) * np.exp(-h/7300) 
  # Define the parameters for constant temperature and composition 
  constant_temperature = 250.0  
  specific_gas_constant = 300.0  
  ratio_of_specific_heats = 1.4  
  # Create the atmosphere settings and add to body settings of "Earth"  
  body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.custom_constant_temperature( 
      density_function, 
      constant_temperature, 
      specific_gas_constant,  
      ratio_of_specific_heats ) 


    )";



    } else if(name == "scaled_by_constant" && variant==0) {
        return R"(
        
Factory function for creating scaled atmospheric model settings.

Factory function for settings object, defining atmospheric model based on an scaling of an existing atmospheric settings object.
The user can apply a scaling factor (or an absolute value) to the air densities of the existing model settings (for instance for an uncertainty analysis).


Parameters
----------
unscaled_atmosphere_settings : AtmosphereSettings
    Sets base settings of atmosphere model to be scaled.
density_scaling : float
    Constant scaling factor to be applied to the entire air density profile.
is_scaling_absolute : bool, default=false
    Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.

Returns
-------
ScaledAtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.



Notes
-----
At present, the scaled atmosphere model only supports scaling of the density value.
For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
this calculation is performed using the `unscaled` density!



Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
by modifying an existing :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` object such that the resulting air density profile is scaled by a constant:

.. code-block:: python 
  
  # define parameter for scaling 
  scaling_constant = 1.5 
  # define variable containing the existing atmosphere model settings  
  unscaled_atmosphere_settings = body_settings.get( "Earth" ).atmosphere_settings 
  # create atmosphere settings and add to body settings of "Earth"  
  body_settings.get( "Earth" ).atmosphere_settings =  environment_setup.atmosphere.scaled_by_constant( 
      unscaled_atmosphere_settings, 
      scaling_constant ) 


    )";



    } else if(name == "scaled_by_function" && variant==0) {
        return R"(
        
Factory function for creating scaled atmospheric model settings.

Factory function for settings object, defining atmospheric model based on scaling an existing atmospheric settings object.
The user can apply custom scaling factors (or absolute values) to the air densities of the existing model settings (for instance for an uncertainty analysis).


Parameters
----------
unscaled_atmosphere_settings : AtmosphereSettings
    Sets base settings of atmosphere model to be scaled.
density_scaling_function : Callable[[float], float]
    Specifies air density scaling factor as a function of time.
is_scaling_absolute : bool, default=false
    Boolean indicating whether density scaling is absolute. Setting this boolean to true will add the scaling value to the baseline density, instead of the default behaviour of multiplying the baseline density by the scaling value.

Returns
-------
ScaledAtmosphereSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.ScaledAtmosphereSettings` class.



Notes
-----
At present, the scaled atmosphere model only supports scaling of the density value.
For cases where the density is used to compute other atmospheric quantities (such as pressure using hydrostatic equilibrium),
this calculation is performed using the `unscaled` density!



Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` for Earth,
by modifying an existing :class:`~tudatpy.numerical_simulation.environment_setup.atmosphere.AtmosphereSettings` object such, that the resulting air density profile is scaled with a user-defined function of time:

.. code-block:: python 
  
  # Define the density scaling as a function of time [sec since J2000] (to add noise with a sine) 
  def scaling_function(time): 
      return 1 + np.sin(time / 50) * 0.25
  # Extract the existing atmosphere model settings 
  unscaled_atmosphere_settings = body_settings.get( "Earth" ).atmosphere_settings 
  # Create the atmosphere settings and add to body settings of "Earth"  
  body_settings.get( "Earth" ).atmosphere_settings =  environment_setup.atmosphere.scaled_by_function( 
      unscaled_atmosphere_settings, 
      scaling_function ) 


    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace ephemeris {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "EphemerisSettings") {
         return R"(

        Base class for providing settings for ephemeris model.

        Functional (base) class for settings of ephemeris models that require no information in addition to their type (and frame origin and orientation).
        Ephemeris model classes requiring additional information must be created using an object derived from this class.





     )";


    } else if(name == "EphemerisSettings.ephemeris_type") {
         return R"(

        **read-only**

        Type of ephemeris that is to be created.

        :type: EphemerisType
     )";


    } else if(name == "EphemerisSettings.frame_origin") {
         return R"(

        Origin of frame in which ephemeris data is to be defined.

        :type: str
     )";


    } else if(name == "EphemerisSettings.frame_orientation") {
         return R"(

        Orientation of frame in which ephemeris data is to be defined.

        :type: str
     )";


    } else if(name == "EphemerisSettings.make_multi_arc_ephemeris") {
         return R"(

        Boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris.

        :type: bool
     )";





    } else if(name == "ScaledEphemerisSettings") {
         return R"(

        Class for defining settings from scaling existing ephemeris settings.

        `EphemerisSettings` derived class for a new ephemeris created from scaling an existing ephemeris settings object. It allows the user to apply a scaling factor to the resulting Cartesian states (for instance for an uncertainty analysis).




     )";





    } else if(name == "DirectSpiceEphemerisSettings") {
         return R"(

        Class for defining settings of an ephemeris linked directly to Spice.

        `EphemerisSettings` derived class for ephemeris which are directly linked to Spice.




     )";


    } else if(name == "DirectSpiceEphemerisSettings.correct_for_stellar_aberration") {
         return R"(

        **read-only**

        Boolean defining whether to correct for stellar aberrations in retrieved values (of observed state).

        :type: bool
     )";


    } else if(name == "DirectSpiceEphemerisSettings.correct_for_light_time_aberration") {
         return R"(

        **read-only**

        Boolean defining whether to correct for light time in retrieved values (of observed state).

        :type: bool
     )";


    } else if(name == "DirectSpiceEphemerisSettings.converge_light_time_aberration") {
         return R"(

        **read-only**

        Boolean defining whether to use single iteration or max. 3 iterations for calculating light time correction.

        :type: bool
     )";





    } else if(name == "InterpolatedSpiceEphemerisSettings") {
         return R"(

        Class for defining settings of an ephemeris interpolated from Spice data.

        `DirectSpiceEphemerisSettings` derived class for setting ephemerides to be created from interpolated Spice ephemeris data.




     )";


    } else if(name == "InterpolatedSpiceEphemerisSettings.initial_time") {
         return R"(

        **read-only**

        Initial time from which interpolated data from Spice should be created.

        :type: float
     )";


    } else if(name == "InterpolatedSpiceEphemerisSettings.final_time") {
         return R"(

        **read-only**

        Final time from which interpolated data from Spice should be created.

        :type: float
     )";


    } else if(name == "InterpolatedSpiceEphemerisSettings.time_step") {
         return R"(

        **read-only**

        Time step setting to be used for the state interpolation.

        :type: float
     )";





    } else if(name == "ApproximateJplEphemerisSettings") {
         return R"(

        Class for creating settings of approximate ephemeris for major planets.

        `EphemerisSettings` derived class for approximate ephemeris for major planets as implemented in ApproximateJplEphemerisSettings class and derived class (described in `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_).




     )";





    } else if(name == "ConstantEphemerisSettings") {
         return R"(

        Class for defining settings of constant ephemerides.

        `EphemerisSettings` derived class for ephemerides producing a constant (time-independent) state.




     )";





    } else if(name == "CustomEphemerisSettings") {
         return R"(

        Class for defining settings of a custom ephemeris.

        `EphemerisSettings` derived class for ephemerides which represent an ideal Kepler orbit.




     )";


    } else if(name == "CustomEphemerisSettings.initial_state_in_keplerian_elements") {
         return R"(

        **read-only**

        Kepler elements at time epochOfInitialState.

        :type: numpy.ndarray[numpy.float64[6, 1]]
     )";


    } else if(name == "CustomEphemerisSettings.epoch_of_initial_state") {
         return R"(

        **read-only**

        Time at which initialStateInKeplerianElements represents the Keplerian state.

        :type: float
     )";


    } else if(name == "CustomEphemerisSettings.central_body_gravitational_parameter") {
         return R"(

        **read-only**

        Gravitational parameter of the central body that is used in the computations.

        :type: float
     )";


    } else if(name == "CustomEphemerisSettings.root_finder_absolute_tolerance") {
         return R"(

        **read-only**

        Convergence tolerance on iterative conversion from mean to eccentric anomaly;
        applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.


        :type: float
     )";


    } else if(name == "CustomEphemerisSettings.root_finder_maximum_number_of_iterations") {
         return R"(

        **read-only**

        Maximum iteration on iterative conversion from mean to eccentric anomaly;
        applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.


        :type: float
     )";





    } else if(name == "TabulatedEphemerisSettings") {
         return R"(

        Class for defining settings of ephemeris to be created from tabulated data.

        `EphemerisSettings` derived class for ephemeris created from tabulated data. The provided data is interpolated into ephemerides.




     )";


    } else if(name == "TabulatedEphemerisSettings.body_state_history") {
         return R"(

        **read-only**

        Dictionary of the discrete state history data from which ephemeris is to be created.

        :type: Dict[[float], numpy.ndarray[numpy.float64[6, 1]]]
     )";


    } else if(name == "TabulatedEphemerisSettings.use_long_double_states") {
         return R"(

        **read-only**

        Boolean defining whether increased numerical precision (long double type) is to be used when creating the ephemeris.

        :type: bool
     )";






    } else if(name == "direct_spice" && variant==0) {
        return R"(
        
Factory function for creating ephemeris model settings entirely from Spice.

Factory function for settings object, defining ephemeris model directly and entirely from Spice.
Requires an appropriate Spice kernel to be loaded.


Parameters
----------
frame_origin : str, default="SSB"
    Origin of frame in which ephemeris data is defined.
frame_orientation : str, default="ECLIPJ2000"
    Orientation of frame in which ephemeris data is defined.
body_name_to_use : str, default = ""
    Body from which Spice ephemeris is to be created.
Returns
-------
DirectSpiceEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.DirectSpiceEphemerisSettings` class





Examples
--------
In this example, we create barycentric (origin: SSB) :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` with axes along J2000, using data directly from spice:

.. code-block:: python 
   
   frame_origin = "SSB" 
   frame_orientation = "J2000" 
   body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.direct_spice( 
   	frame_origin, frame_orientation) 


Alternatively, we can assign the DirectSpiceEphemerisSettings of Jupiter (or any other body for which a direct Spice ephemeris is available) to any custom body:

.. code-block:: python 
   
   frame_origin = "SSB" 
   frame_orientation = "J2000" 
   body_name_to_use =  "Jupiter" 
   # create ephemeris settings from "Jupiter" spice data and add to body settings of body "CustomBody" 
   body_settings.get( "CustomBody" ).ephemeris_settings = environment_setup.ephemeris.direct_spice( 
   	frame_origin, frame_orientation, body_name_to_use ) 


    )";



    } else if(name == "interpolated_spice" && variant==0) {
        return R"(
        
Factory function for creating ephemeris model settings using interpolated Spice data.

Factory function for settings object defining an ephemeris model from interpolated Spice data.
Using this option the state of the body is retrieved from Spice at regular intervals `before` the environment propagation (as opposed to during the propagation).
These data are then used to create an interpolator, which is put into the environment, and called during the propagation.
This option has the downside of being applicable only during a limited time interval and requiring the tabulated data to be stored in RAM,
but may for `some special cases <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models/default_bodies_limited_time_range.html>`_
offer an advantage over a direct Spice ephemeris (:func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.direct_spice`).


Parameters
----------
initial_time : float
    Initial time from which interpolated data from Spice should be created.
final_time : float
    Final time from which interpolated data from Spice should be created.
time_step : float
    Time step with which interpolated data from Spice should be created.
frame_origin : str, default="SSB"
    Origin of frame in which ephemeris data is defined.
frame_orientation : str, default="ECLIPJ2000"
    Orientation of frame in which ephemeris data is defined.
interpolator_settings : std::make_shared< interpolators::InterpolatorSettings >, default=std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 )
    Settings to be used for the state interpolation.
body_name_to_use : str, default = ""
    Body from which Spice ephemeris is to be created.
Returns
-------
InterpolatedSpiceEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.DirectSpiceEphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.InterpolatedSpiceEphemerisSettings` class





Examples
--------
In this example, we define :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for Jupiter by retrieving ephemeris data from Spice at 3600 s intervals between t=0 and t=1.0E8:

.. code-block:: python 
  
  # Define the interpolation settings
  initial_time = 0.0 
  final_time = 1.0E8 
  time_step = 3600.0 
  # Define the ephemeris frame
  frame_origin = "SSB" 
  frame_orientation = "J2000" 
  # create ephemeris settings and add to body settings of body "Jupiter"
  body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.interpolated_spice( 
    initial_time, final_time, time_step, frame_origin, frame_orientation ) 


By default, a 6th order Lagrange interpolator is used (NOTE: the Lagrange interpolator is not reliable at the edges of the interpolation interval, as discussed here: :func:`~tudatpy.math.interpolators.lagrange_interpolation`).
Settings for an alternative interpolator can be use by specifying the optional input argument.
Additionally, as is the case for the :func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.direct_spice` and :func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.approximate_jpl_model` functions, an optional input argument ``body_name_to_use`` allows to use an ephemeris model from Spice for some body and assign it to a custom body.


    )";



    } else if(name == "approximate_jpl_model" && variant==0) {
        return R"(
        
Factory function for creating approximate ephemeris model settings for major planets.

Factory function for settings object, defining approximate ephemeris model for major planets.
In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described in `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_).
Note that this option is only available for solar system planets. For the case of the Earth the approximate ephemeris of the Earth-Moon barycenter is returned.


Parameters
----------
body_name : str
    String that is attempted to be matched to an identifier for the body that the ephemeris is to be created for.
Returns
-------
ApproximateJplEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ApproximateJplEphemerisSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for Jupiter using JPL's approximate planet position model:

.. code-block:: python 
   
   # create ephemeris settings and add to body settings of body "Jupiter" 
   body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.approximate_jpl_model( "Jupiter" ) 


Alternatively, we can assign the ApproximateJplEphemerisSettings of Jupiter (or any other body for which an approximate JPL ephemeris is available) to any custom body:

.. code-block:: python 
   
   # create ephemeris settings (normally used for Jupiter) and add to body settings of body "CustomBody" 
   body_settings.get( "CustomBody" ).ephemeris_settings = environment_setup.ephemeris.approximate_jpl_model( "Jupiter" ) 

   							ephemerides::ApproximatePlanetPositionsBase::jupiter, false );   


    )";



    } else if(name == "constant" && variant==0) {
        return R"(
        
Factory function for creating constant ephemeris model settings.

Factory function for settings object, defining ephemeris model with a constant, time-independent state.


Parameters
----------
constant_state : numpy.ndarray[numpy.float64[6, 1]]
    Constant state that will be provided as output of the ephemeris at all times.
frame_origin : str, default="SSB"
    Origin of frame in which ephemeris data is defined.
frame_orientation : str, default="ECLIPJ2000"
    Orientation of frame in which ephemeris data is defined.
Returns
-------
ConstantEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ConstantEphemerisSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for a time-independent, constant state of Jupiter:

.. code-block:: python 
   
   # Define the constant cartesian state 
   constant_cartesian_state = [100.0e9, 100.0e9, 100.0e9, 10.0e3, 10.0e3, 10.0e3] 
   # Define the ephemeris frame 
   frame_origin = "SSB" 
   frame_orientation = "J2000" 
   # Make the ephemeris settings 
   body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.constant( 
     constant_cartesian_state, 
     frame_origin, frame_orientation) 


    )";



    } else if(name == "custom_ephemeris" && variant==0) {
        return R"(
        
Factory function for creating custom ephemeris model settings.

Factory function for settings object, defining ephemeris model with a custom state.
This allows the user to provide a custom state function as ephemeris model.


Parameters
----------
custom_state_function : Callable[[float], numpy.ndarray[numpy.float64[6, 1]]]
    Function returning the state as a function of time.
frame_origin : str, default="SSB"
    Origin of frame in which ephemeris data is defined.
frame_orientation : str, default="ECLIPJ2000"
    Orientation of frame in which ephemeris data is defined.
Returns
-------
CustomEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.CustomEphemerisSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for Earth from a custom state history function:

.. code-block:: python 
  
  # Define the custom state function for Earth 
  def custom_state_function(time): 
      # Compute what fraction of the year it is 
      frac_year = (time - 2451545) % (365.25*24*3600) 
      # Distance and velocity of the Earth w.r.t. the Sun 
      AU, v_E = 1.496e11, 30e3 
      # Compute the position and velocity of the Earth in a 2D circle 
      x_pos = np.sin(frac_year*np.pi) * AU 
      y_pos = np.cos(frac_year*np.pi) * AU 
      x_vel = np.cos(frac_year*np.pi) * v_E 
      y_vel = np.sin(frac_year*np.pi) * v_E 
      return [x_pos, y_pos, 0, x_vel, y_vel, 0] 
  # Define the ephemeris frame 
  frame_origin = "SSB" 
  frame_orientation = "J2000" 
  # Make the ephemeris settings 
  body_settings.get("Earth").ephemeris_settings = environment_setup.ephemeris.custom( 
      custom_state_function, 
      frame_origin, 
      frame_orientation) 


    )";



    } else if(name == "keplerian" && variant==0) {
        return R"(
        
Factory function for creating Keplerian ephemeris model settings.

Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from the given Kepler elements.
These are taken as the elements at the ``initial_state_epoch`` and propagated to any other time using the provided ``central_body_gravitational_parameter``.
See `Element Types <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#element-types>`_ and the :ref:`\`\`astro\`\`` module for more details on orbital elements in tudat.


Parameters
----------
initial_state_in_keplerian_elements : numpy.ndarray[numpy.float64[6, 1]]
    Kepler elements at epoch given by ``initial_state_epoch``. 

initial_state_epoch : float
    Epoch at which ``initial_state_epoch`` represents the Keplerian state. 

central_body_gravitational_parameter : float
    Effective gravitational parameter of the central body that is used in the computations. Note that when
    the Keplerian orbit is to represent the relative state of two massive bodies, with one of these bodies as the origin
    this values should be the *sum* of the two bodies' gravitational parameters

frame_origin : str, default="SSB"
    Origin of frame in which ephemeris data is defined.
frame_orientation : str, default="ECLIPJ2000"
    Orientation of frame in which ephemeris data is defined.
root_finder_absolute_tolerance : float
    Convergence tolerance on iterative conversion from mean to eccentric anomaly;
    applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

root_finder_maximum_number_of_iterations : float
    Maximum iteration on iterative conversion from mean to eccentric anomaly;
    applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

Returns
-------
KeplerEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.KeplerEphemerisSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for a simple, barycentric (SSB) Kepler orbit of Jupiter:

.. code-block:: python 
  
  # Define the computation of the Kepler orbit ephemeris 
  initial_state_in_keplerian_elements = [100e9, 0.7, 1.0, 2.0, 2.0, 2.0] 
  initial_state_epoch = 12345 
  central_body_gravitational_parameter = 1.3284e20 # (sum of Sun and Jupiter) 
  # Define the ephemeris frame 
  frame_origin = "SSB" 
  frame_orientation = "J2000" 
  # Create ephemeris settings and add to body settings of "Jupiter" 
  body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.keplerian( 
      initial_state_in_keplerian_elements, 
      initial_state_epoch, 
      central_body_gravitational_parameter, 
      frame_origin, frame_orientation ) 


    )";



    } else if(name == "keplerian_from_spice" && variant==0) {
        return R"(
        
Factory function for creating Keplerian ephemeris model settings with initial state from Spice.

Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from an initial state from Spice.
The Kepler elements inferred from the initial state are propagated to any other time using the provided ``central_body_gravitational_parameter``.
See `Element Types <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#element-types>`_ and the :ref:`\`\`astro\`\`` module for more details on orbital elements in tudat.


Parameters
----------
body : str
    Name of body for which to create ephemeris settings and infer initial state from Spice.
initial_state_epoch : float
    Epoch at which ``initial_state_epoch`` represents the Keplerian state. 

central_body_gravitational_parameter : float
    Gravitational parameter of the central body that is used in the computations.
frame_origin : str, default="SSB"
    Origin of frame in which ephemeris data is defined.
frame_orientation : str, default="ECLIPJ2000"
    Orientation of frame in which ephemeris data is defined.
root_finder_absolute_tolerance : float
    Convergence tolerance on iterative conversion from mean to eccentric anomaly;
    applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

root_finder_maximum_number_of_iterations : float
    Maximum iteration on iterative conversion from mean to eccentric anomaly;
    applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

Returns
-------
KeplerEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.KeplerEphemerisSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for a simple, barycentric (SSB) Kepler orbit of Jupiter.
The initial keplerian state is extracted from Spice as the state of ``body_name`` w.r.t. ``frame_origin``

.. code-block:: python 
  
  # Define the parameters for retrieval of the initial Kepler orbit elements from spice 
  body_name = "Jupiter" 
  initial_state_epoch = 12345 
  central_body_gravitational_parameter = 1.3284e20 # (sum of Sun and Jupiter) 
  # Define the ephemeris frame 
  frame_origin = "SSB" 
  frame_orientation = 'J2000' 
  # Make ephemeris the settings and add to body settings of "Jupiter" 
  body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.keplerian_from_spice(  
      body_name, 
      initial_state_epoch, 
      central_body_gravitational_parameter, 
      frame_origin, 
      frame_orientation ) 


Additionally, as is the case for the :func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.direct_spice`, :func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.approximate_jpl_model` and :func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.interpolated_spice` functions, the ephemeris model from Spice can be retrieved for some body and assigned to a custom body.


    )";



    } else if(name == "scaled_by_constant" && variant==0) {
        return R"(
        
Factory function for creating scaled ephemeris model settings.

Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).


Parameters
----------
unscaled_ephemeris_settings : EphemerisSettings
    Sets base settings of ephemeris to be scaled.
scaling_constant : float
    Constant scaling factor to be applied to all elements of the Cartesian state.
is_scaling_absolute : bool, default=false
    Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
Returns
-------
ScaledEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class





Examples
--------
In this example, we create ephemeris settings for Jupiter, by scaling an existing :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettingsObject` with a constant factor:

.. code-block:: python 
   
   # define variable for scaling factor
   scaling_constant = 1.001
   # define variables containing the existing ephemeris settings
   unscaled_ephemeris_settings = body_settings.get( "Jupiter" ).ephemeris_settings
   # make new ephemeris settings
   body_settings.get( "Jupiter" ).ephemeris_settings =  environment_setup.ephemeris.scaled_by_constant(
          unscaled_ephemeris_settings, scaling_constant )

In the above case, the original Jupiter ephemeris setting is taken and each state element (x,y,z position and velocity) from the original ephemeris is multiplied by a factor 1.001.


    )";



    } else if(name == "scaled_by_vector" && variant==0) {
        return R"(
        
Factory function for creating scaled ephemeris model settings.

Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).


Parameters
----------
unscaled_ephemeris_settings : EphemerisSettings
    Sets base settings of ephemeris to be scaled.
scaling_vector : numpy.ndarray[numpy.float64[6, 1]]
    Vector containing scaling factors to be applied to each element of the Cartesian state.
is_scaling_absolute : bool, default=false
    Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
Returns
-------
ScaledEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class





Examples
--------
In this example, we create ephemeris settings for Jupiter, by scaling an existing :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettingsObject` with the constant elements of a vector:

.. code-block:: python 
  
  # Define the scaling vector 
  scaling_vector = [1.01, 0.99, 1, 1, 1, 0] 
  # Extract the unscaled ephemeris settings from Jupiter 
  unscaled_ephemeris_settings = body_settings.get( "Jupiter" ).ephemeris_settings 
  # Create the scaled ephemeris settings and apply to the body "Jupiter" 
  body_settings.get( "Jupiter" ).ephemeris_settings =  environment_setup.ephemeris.scaled_by_vector( 
      unscaled_ephemeris_settings, 
      scaling_vector) 

In the above case, the original Jupiter ephemeris setting is taken and each state element (x,y,z position and velocity) from the original ephemeris is multiplied by the corresponding scaling factor in ``scaling_vector``.


    )";



    } else if(name == "scaled_by_vector_function" && variant==0) {
        return R"(
        
Factory function for creating scaled ephemeris model settings.

Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).


Parameters
----------
unscaled_ephemeris_settings : EphemerisSettings
    Sets base settings of ephemeris to be scaled.
scaling_vector_function : callable[[float], numpy.ndarray[numpy.float64[6, 1]]]
    Function returning a vector with the scaling factors to be applied to each element of the Cartesian state.
is_scaling_absolute : bool, default=false
    Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
Returns
-------
ScaledEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class





Examples
--------
In this example, we create ephemeris settings for Jupiter, by scaling an existing :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` object with factors from a custom function:

.. code-block:: python 
  
  # Define the scaling vector function 
  def scaling_vector_function(time): 
      # Add a wobble in the x and y coordinates 
      wobble = 1 + 0.1 * np.cos(time/50) 
      return [wobble, wobble, 1, 1, 1, 1] 
  # Extract the existing unscaled ephemeris settings 
  unscaled_ephemeris_settings = body_settings.get( "Jupiter" ).ephemeris_settings 
  # Create the scaled ephemeris settings and apply to the body "Jupiter" 
  body_settings.get( "Jupiter" ).ephemeris_settings =  environment_setup.ephemeris.scaled_by_vector_function( 
      unscaled_ephemeris_settings, 
      scaling_vector_function ) 

In the above case, the original Jupiter ephemeris setting is taken and each state element (x,y,z position and velocity) from the original ephemeris is multiplied by the corresponding scaling factor in the vector returned by ``vector_scaling_function``.


    )";



    } else if(name == "tabulated" && variant==0) {
        return R"(
        
Factory function for creating ephemeris model settings from tabulated data.

Factory function for settings object, defining ephemeris model to be created from tabulated data.
Currently the data that is provided gets interpolated by a 6th order Lagrange interpolator (hardcoded).
At the edges of the interpolation interval a cubic spline interpolator is used to suppress the influence of Runge's phenomenon.


Parameters
----------
body_state_history : dict
    Dictionary of the discrete state history data from which ephemeris is to be created. Keys representing the time (float) and values representing Cartesian states (numpy.ndarray).
frame_origin : str, default="SSB"
    Origin of frame in which ephemeris data is defined.
frame_orientation : str, default="ECLIPJ2000"
    Orientation of frame in which ephemeris data is defined.
Returns
-------
TabulatedEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.TabulatedEphemerisSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for Jupiter from tabulated state history data:

.. code-block:: python 
  
  # Define the Dict containing Jupiter's tabulated state history 
  body_state_history = { 
      0: [7.4713e11, 0, 0, 13.5e3, 0, 0], 
      1000: [7.4711e11, 5e9, 0, 13.4998e3, 75, 0], 
      2150: [7.4671e11, 2.5e10, 0, 13.498e3, 200, 0], 
      # ... truncated 
      15650: [7.3899e11, 1.1e11, 0, 13.416e3, 1.5e3, 0] 
  } 
  # Define the ephemeris frame 
  frame_origin = "SSB" 
  frame_orientation = "J2000" 
  # Create the tabulated ephemeris settings and add them to the body "Jupiter" 
  body_settings.get( "Jupiter" ).ephemeris_settings = environment_setup.ephemeris.tabulated( body_state_history, 
      frame_origin, 
      frame_orientation ) 


    )";



    } else if(name == "tabulated_from_existing" && variant==0) {
        return R"(
        
Factory function for creating tabulated ephemeris model settings from existing ephemeris.

Factory function for creating tabulated ephemeris model settings from existing ephemeris.
The ephemeris that is provided gets tabulated in a given time frame, for a given time step.
When called, this tabulated ephemeris will use interpolation, when needed, from the specified interpolator.

.. note:: Creating tabulated ephemeris from existing ephemeris can for instance be used when combined with estimation.
          This is because estimation needs the ephemeris to be tabulated to work.


Parameters
----------
ephemeris_settings : tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings
    Existing ephemeris settings that have to be tabulated.
start_time : float
    Initial time for which to create the tabulated ephemeris.
end_time : float
    Final time for which to create the tabulated ephemeris.
time_step : float
    Time step to use to tabulate the existing ephemeris.
interpolator_settings : tudatpy.math.interpolators.InterpolatorSettings, default=tudatpy.math.interpolators.lagrange_interpolation(8)
    Interpolator settings to use when interpolating between two tabulated ephemeris.
Returns
-------
TabulatedEphemerisSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.TabulatedEphemerisSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` for Io.
First, we extract the existing ephemeris. Then, we define new tabulated ephemeris settings, from the original settings.

.. code-block:: python 
  
  # Get the original ephemeris settings 
  original_io_ephemeris_settings = body_settings.get( "Io" ).ephemeris_settings 
  # Apply new tabulated ephemeris settings 
  body_settings.get( "Io" ).ephemeris_settings =  environment_setup.ephemeris.tabulated_from_existing( 
    original_io_ephemeris_settings, 
    initial_time, 
    final_time, 
    time_step ) 


    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace gravity_field {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "GravityFieldType") {
         return R"(

        Enumeration of gravity field types.

        Enumeration of gravity field types supported by tudat.





     )";


    } else if(name == "GravityFieldType.polyhedron") {
         return R"(
     )";


    } else if(name == "GravityFieldType.central_gravity") {
         return R"(
     )";


    } else if(name == "GravityFieldType.central_spice_gravity") {
         return R"(
     )";


    } else if(name == "GravityFieldType.spherical_harmonic_gravity") {
         return R"(
     )";


    } else if(name == "GravityFieldType.polyhedron_gravity") {
         return R"(
     )";



    } else if(name == "PredefinedSphericalHarmonicsModel") {
         return R"(

        Enumeration of predefined spherical harmonics models.

        Enumeration of predefined spherical harmonics models supported by tudat, for which thee coefficient files are automatically available (downloaded from
        `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_). The directory where these files are stored can be 
        extracted using the :func:`~tudatpy.io.get_gravity_models_path` function.





     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.egm96") {
         return R"(
Coefficients for EGM96 Earth gravity field up to degree and order 200, (see `link <https://cddis.gsfc.nasa.gov/926/egm96/egm96.html>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.ggm02c") {
         return R"(
Coefficients for the combined GGM02 Earth gravity field up to degree and order 200, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.ggm02s") {
         return R"(
Coefficients for the GRACE-only GGM02 Earth gravity field up to degree and order 160, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.goco05c") {
         return R"(
Coefficients for the GOCO05c combined Earth gravity field up to degree and order 719, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.glgm3150") {
         return R"(
Coefficients for the GLGM3150 Moon gravity field up to degree and order 150, (see `link <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.lpe200") {
         return R"(
Coefficients for the LPE200 Moon gravity field up to degree and order 200, (see `link <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.gggrx1200") {
         return R"(
Coefficients for the GRGM1200A Moon gravity field up to degree and order 1199, (see `link <https://pgda.gsfc.nasa.gov/products/50>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.jgmro120d") {
         return R"(
Coefficients for the MRO120D Moon gravity field up to degree and order 120, (see `link <https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.jgmess160a") {
         return R"(
Coefficients for the MESS160A Moon gravity field up to degree and order 160, (see `link <https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.lbl>`_ )
     )";


    } else if(name == "PredefinedSphericalHarmonicsModel.shgj180u") {
         return R"(
Coefficients for the SHGJ180U Moon gravity field up to degree and order 180, (see `link <https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/shgj120u.lbl>`_ )
     )";




    } else if(name == "GravityFieldSettings") {
         return R"(

        Base class for providing settings for automatic gravity field model creation.

        This class is a functional base class for settings of gravity field models that require no information in addition to their type.
        Gravity field model classes requiring additional information must be created using an object derived from this class.





     )";


    } else if(name == "GravityFieldSettings.gravity_field_type") {
         return R"(

        **read-only**

        Type of gravity field model that is to be created.

        :type: GravityFieldType
     )";




    } else if(name == "GravityFieldSettings.__init__" && variant==0) {
            return R"(





    )";




    } else if(name == "CentralGravityFieldSettings") {
         return R"(

        `GravityFieldSettings` derived class defining settings of point mass gravity field.

        Derived class of `GravityFieldSettings` for central gravity fields, which are defined by a single gravitational parameter.





     )";


    } else if(name == "CentralGravityFieldSettings.gravitational_parameter") {
         return R"(

        Gravitational parameter of central gravity field.

        :type: float
     )";





    } else if(name == "SphericalHarmonicsGravityFieldSettings") {
         return R"(

        `GravityFieldSettings` derived class defining settings of spherical harmonic gravity field representation.

        Derived class of `GravityFieldSettings` for gravity fields, which are defined by a spherical harmonic gravity field representation.





     )";


    } else if(name == "SphericalHarmonicsGravityFieldSettings.gravitational_parameter") {
         return R"(

        Gravitational parameter of gravity field.

        :type: float
     )";


    } else if(name == "SphericalHarmonicsGravityFieldSettings.reference_radius") {
         return R"(

        **read-only**

        Reference radius of spherical harmonic field expansion.

        :type: float
     )";


    } else if(name == "SphericalHarmonicsGravityFieldSettings.normalized_cosine_coefficients") {
         return R"(

        Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.

        :type: numpy.ndarray
     )";


    } else if(name == "SphericalHarmonicsGravityFieldSettings.normalized_sine_coefficients") {
         return R"(

        Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.

        :type: numpy.ndarray
     )";


    } else if(name == "SphericalHarmonicsGravityFieldSettings.associated_reference_frame") {
         return R"(

        Identifier for body-fixed reference frame with which the coefficients are associated.

        :type: str
     )";


    } else if(name == "SphericalHarmonicsGravityFieldSettings.create_time_dependent_field") {
         return R"(

        Boolean that denotes whether the field should be created as time-dependent (even if no variations are imposed initially).

        :type: bool
     )";


    } else if(name == "SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia") {
         return R"(

        Value of the scaled mean moment of inertia :math:`I_{xx}+I_{yy}+I_{zz}/(MR^{2})`. This value does not influence the gravity field itself,
        but together with the degree 2 gravity field coefficients defines the body's inertia tensor.


        :type: float
     )";





    } else if(name == "PolyhedronGravityFieldSettings") {
         return R"(

        `GravityFieldSettings` derived class defining settings of a polyhedron gravity field representation.

        Derived class of `GravityFieldSettings` for gravity fields, which are defined by a polyhedron gravity field representation.





     )";


    } else if(name == "PolyhedronGravityFieldSettings.gravitational_parameter") {
         return R"(

        Gravitational parameter of gravity field.

        :type: float
     )";


    } else if(name == "PolyhedronGravityFieldSettings.density") {
         return R"(

        Density of the polyhedron.

        :type: float
     )";


    } else if(name == "PolyhedronGravityFieldSettings.associated_reference_frame") {
         return R"(

        Identifier for body-fixed reference frame with which the vertices coordinates are associated.

        :type: str
     )";


    } else if(name == "PolyhedronGravityFieldSettings.vertices_coordinates") {
         return R"(

        Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
        row per vertex, 3 columns). 


        :type: numpy.ndarray
     )";


    } else if(name == "PolyhedronGravityFieldSettings.vertices_defining_each_facet") {
         return R"(

        Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
        the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
        when seen from the outside of the polyhedron. 


        :type: numpy.ndarray
     )";






    } else if(name == "central" && variant==0) {
        return R"(
        
Factory function for central gravity field settings object.

Factory function for settings object, defining a point-mass gravity field model with user-defined gravitational parameter :math:`\mu`. The gravitational potential is the defined as:

.. math::
   U(\mathbf{r})=\frac{\mu}{||\mathbf{r}||}

with :math:`\mathbf{r}` the position vector measured from the body's center of mass.


Parameters
----------
gravitational_parameter : float
    Gravitational parameter defining the point-mass gravity field.
Returns
-------
CentralGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.CentralGravityFieldSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using a simple central gravity field model:

.. code-block:: python 
   
   # define parameters describing central gravity model 
   gravitational_parameter = 3.986e14 
   # create gravity field settings 
   body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.central( gravitational_parameter ) 


    )";



    } else if(name == "central_spice" && variant==0) {
        return R"(
        
Factory function to create central gravity field settings from Spice settings.

Factory function for settings object, defining a point-mass gravity field model. This function provides the same model as :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.central`), but with gravitational parameter :math:`\mu` from Spice.

Returns
-------
GravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` class of gravity field type ``central_spice``





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using a simple central gravity field model and data from Spice:

.. code-block:: python 
   
   # create gravity field settings 
   body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.central_spice( )


    )";



    } else if(name == "spherical_harmonic" && variant==0) {
        return R"(
        
Factory function for creating a spherical harmonics gravity field settings object.

Factory function for settings object, defining a gravity field model through spherical harmonic expansion.
The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
It represents the frame in which the spherical harmonic field is defined.

The gravitational potential is the defined as:

.. math::
   U(\mathbf{r})=\sum_{l=0}^{l_{max}}\sum_{m=0}^{l}\mu\left(\frac{{R}^{l}}{r^{l+1}}\right)\bar{P}_{lm}(\sin\phi)\left(\bar{C}_{lm}\cos m\theta+\bar{S}_{lm}\sin m\theta\right)

with :math:`\mathbf{r}` the position vector of the evaluation point, measured from the body's center of mass. The angles :math:`\phi` and :math:`\theta` are the body-fixed latitude and longitude of the evaluation point, and :math:`\bar{P}_{lm}` is the associated Legendre polynomial (at degree/order :math:`l/m`).

For the spherical harmonic gravity field (including other spherical harmonic functions), the normalized mean moment of inertia must be set by the user, to allow an inertia tensor to be computed. This is done using the :attr:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia` attribute of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class, as in the example below

.. code-block:: python 

  # Add gravity field model settings to body of spherical harmonic type 
  body_settings.get( "Mars" ).gravity_field = ... 

  # Add setting for moment of inertia 
  body_settings.get( "Mars" ).gravity_field.scaled_mean_moment_of_inertia = 0.365 
       
This code snippet will automatically create a rigid body properties for Mars, with the inertia tensor computed from this value of 0.365 and the degree 2 gravity field coefficients. Note that, if gravity field variations are used for the body, time-variability of the degree 1- and 2- coefficients will be reflected in time-variability of the body's center of mass and inertia tensor.

Note: Spherical harmonic coefficients used for this environment model must *always* be fully normalized.
To normalize un-normalized spherical harmonic coefficients, see :func:`~tudatpy.astro.gravitation.normalize_spherical_harmonic_coefficients`.


Parameters
----------
gravitational_parameter : float
    Gravitational parameter :math:`\mu` of gravity field.
reference_radius : float
    Reference radius :math:`R` of spherical harmonic field expansion.
normalized_cosine_coefficients : numpy.ndarray
    Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\bar{C}_{ij}` at degree i and order j.
    As such, note that entry (0,0) of cosine coefficients should be equal to 1.

normalized_sine_coefficients : numpy.ndarray
    Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\bar{S}_{ij}`
    at degree i and order j.

associated_reference_frame : str
    Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
Returns
-------
SphericalHarmonicsGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using a spherical harmonics gravity model:

.. code-block:: python 
  
  # Define the spherical harmonics gravity model 
  gravitational_parameter = 3986004.415E+8 
  reference_radius = 6378136.3 
  # Normalized coefficients taken from https://cddis.nasa.gov/archive/egm96/general_info/egm96_to360.ascii 
  # The above file is described in https://cddis.nasa.gov/archive/egm96/general_info/readme.egm96 
  normalized_cosine_coefficients = [ 
      [1,                   0,                   0,                   0], 
      [0,                   0,                   0,                   0], 
      [-0.484165371736E-03, -0.186987635955E-09, 0.243914352398E-05,  0], 
      [0.957254173792E-06,  0.202998882184E-05,  0.904627768605E-06,  0.721072657057E-06] 
  ] 
  normalized_sine_coefficients = [ 
      [0,                   0,                   0,                   0], 
      [0,                   0,                   0,                   0], 
      [0,                   0.119528012031E-08,  -0.140016683654E-05, 0], 
      [0,                   0.248513158716E-06,  -0.619025944205E-06, 0.141435626958E-05] 
  ] 
  associated_reference_frame = "IAU_Earth" 
  # Create the gravity field settings and add them to the body "Earth" 
  body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.spherical_harmonic( 
      gravitational_parameter, 
      reference_radius, 
      normalized_cosine_coefficients, 
      normalized_sine_coefficients, 
      associated_reference_frame ) 


    )";



    } else if(name == "sh_triaxial_ellipsoid_from_density" && variant==0) {
        return R"(
        
Factory function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the density to define the mass distribution.

Factory function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
The constant mass distribution is defined by the density and gravitational constant (optional).
The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
This function implements the models of (see Balmino [1]_).


Parameters
----------
axis_a : float
    Dimension of largest axis of triaxial ellipsoid.
axis_b : float
    Dimension of intermediate axis of triaxial ellipsoid.
axis_c : float
    Dimension of smallest axis of triaxial ellipsoid.
density : float
    Density of ellipsoid.
maximum_degree : int
    Maximum degree of spherical harmonics expansion.
maximum_order : int
    Maximum order of spherical harmonics expansion.
associated_reference_frame : str
    Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
gravitational_constant : float, default=physical_constants::GRAVITATIONAL_CONSTANT
    Gravitational constant G of the gravity field.
Returns
-------
SphericalHarmonicsGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using the expansion of a homogeneous triaxial ellipsoid into a spherical harmonics gravity model:

.. code-block:: python 
  
  # Create the gravity field settings for Earth with Spherical Harmonics from a triaxial ellipsoid 
  body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.spherical_harmonic_triaxial_ellipsoid_from_density( 
      axis_a=6378171.88, 
      axis_b=6378102.03, 
      axis_c=6356752.24, 
      density=5520, 
      maximum_degree=50, 
      maximum_order=50, 
      associated_reference_frame="IAU_Earth" ) 


    )";



    } else if(name == "sh_triaxial_ellipsoid_from_gravitational_parameter" && variant==0) {
        return R"(
        
Factory function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the gravitational parameter to define the mass distribution..

Factory function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`
The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
The constant mass distribution is defined by the gravitational parameter.
The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
This function implements the models of (see Balmino [1]_).


Parameters
----------
axis_a : float
    Dimension of largest axis of triaxial ellipsoid.
axis_b : float
    Dimension of intermediate axis of triaxial ellipsoid.
axis_c : float
    Dimension of smallest axis of triaxial ellipsoid.
maximum_degree : int
    Maximum degree of spherical harmonics expansion.
maximum_order : int
    Maximum order of spherical harmonics expansion.
associated_reference_frame : str
    Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
gravitational_parameter : float
    Gravitational parameter :math:`\mu` of gravity field.
Returns
-------
SphericalHarmonicsGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class






    )";



    } else if(name == "from_file_spherical_harmonic" && variant==0) {
        return R"(
        
Factory function to load a custom spherical harmonics gravity field settings from a file.

Factory function to load a custom spherical harmonics gravity field settings from a file. The file should contain **fully normalized** spherical harmonic coefficients.
The associated gravitational paramerer and reference radius should be given in m^3/s^2 and m, respectively. The file format should be the same as that used for the files
in the directories `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_. Specifically, the file should contain

- The first line should be a series of text blocks (typically numerical data). Two of these blocks (by default the first and second one) should be the gravitational parameter and reference radius, respectively. The text block should be separated by spaces, tabs and/or commas
- Each subsequent line should contain a set of spherical harmonic coefficients (first ordered in ascending order by degree, then in ascending order by order), where the first, second, third and fourth value of the line should be: degree :math:`l`, order :math:`m`, normalized cosine coefficient :math:`\bar{C}_{lm}`, normalized sine coefficient :math:`\bar{S}_{lm}`. Additional entries (for instance with coefficient uncertainties) are ignored. 


Parameters
----------
file : str
    Full file path and name where th gravity field file is located
maximum_degree : int
    Maximum degree of the coefficients that are to be loaded
maximum_order : int
    Maximum order of the coefficients that are to be loaded
associated_reference_frame : str, default = ""
    Name of the body-fixed reference frame to which the gravity field is to be fixed. If left empty, this reference frame will automatically be set to the body-fixed frame defined by this body's rotation (see :ref:`\`\`rotation_model\`\`` for specifying rotation models).
gravitational_parameter_index : int, default = 0
    Index of the values in the file header (first line of file) that contains the gravitational parameter
reference_radius_index : int, default = 1
    Index of the values in the file header (first line of file) that contains the reference radius
Returns
-------
SphericalHarmonicsGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for the Moon using GGGRX spherical harmonics gravity model, up to degree and order 300:

.. code-block:: python 
  
  # Create the gravity field settings for Moon with Spherical Harmonics loaded in from a Spherical Harmonics file 
  body_settings.get("Moon").gravity_field_settings = environment_setup.gravity_field.from_file_spherical_harmonic( 
      r"...\.tudat\resource\gravity_models\Moon\gggrx_1200l_sha.tab", 300, 300) 


    )";



    } else if(name == "predefined_spherical_harmonic" && variant==0) {
        return R"(
        
Factory function for spherical harmonics gravity field settings of a predefined model.

Factory function for spherical harmonics gravity field settings of a predefined model


Parameters
----------
predefined_model : PredefinedSphericalHarmonicsModel
    Identified for gravity field model that is to be loaded
maximum_degree : int, default = -1
    Maximum degree and order to which the coefficients are to be loaded. If value is negative, all coefficients for the specified gravity field are loaded
Returns
-------
SphericalHarmonicsGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` for Earth using EGM96 spherical harmonics gravity model, up to degree and order 32:

.. code-block:: python 
  
  # Create the gravity field settings for Earth with Spherical Harmonics from a triaxial ellipsoid 
  body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.predefined_spherical_harmonic( 
      environment_setup.gravity_field.egm96, 32 ) 


    )";



    } else if(name == "polyhedron_from_mu" && variant==0) {
        return R"(
        
Factory function for creating a polyhedron gravity field settings object, using the gravitational parameter.

Factory function for settings object, defining a gravity field model through a polyhedron.
The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
It represents the frame in which the polyhedron field is defined. 

The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
to Werner and Scheeres [2]_.

This function uses the gravitational parameter to define the gravity field. To instead use the density 
constant see :func:`~tudatpy.astro.gravitation.polyhedron_from_density`. Since both models tend to be computationally intensive, 
it is recommended to use polyhedra with the lowest number of facets that allows meeting the desired accuracy. The number of facets of a polyhedron
model can be reduced using any mesh processing software, for example `PyMeshLab <https://pymeshlab.readthedocs.io/en/latest/>`_.
Additionally, different functions to process a polyhedron are available in `Polyhedron utilities <https://py.api.tudat.space/en/latest/polyhedron_utilities.html>`_.


Parameters
----------
gravitational_parameter : float
    Gravitational parameter :math:`\mu` of gravity field.
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

associated_reference_frame : str
    Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
gravitational_constant : float, default=GRAVITATIONAL_CONSTANT
    Newton's gravitational constant G, used to computed the density

Returns
-------
PolyhedronGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class






    )";



    } else if(name == "polyhedron_from_density" && variant==0) {
        return R"(
        
Factory function for creating a polyhedron gravity field settings object, using the density.

Factory function for settings object, defining a gravity field model through a polyhedron.
The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
It represents the frame in which the polyhedron field is defined.

The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
to Werner and Scheeres [2]_.

This function uses the density to define the gravity field. To instead use the
gravitational parameter see :func:`~tudatpy.astro.gravitation.polyhedron_from_mu`.


Parameters
----------
density : float, default=TUDAT_NAN
    Density of the polyhedron.

vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron. 

associated_reference_frame : str
    Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
gravitational_constant : float, default=GRAVITATIONAL_CONSTANT
    Newton's gravitational constant G, used to computed the gravitational parameter

Returns
-------
PolyhedronGravityFieldSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace gravity_field_variation {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "GravityFieldVariationSettings") {
         return R"(

        Base class for providing settings for gravity field variations.





     )";





    } else if(name == "BasicSolidBodyGravityFieldVariationSettings") {
         return R"(

        Class for providing settings for solid body tidal gravity field variations, derived from GravityFieldVariationSettings.





     )";






    } else if(name == "solid_body_tide" && variant==0) {
        return R"(
        
Factory function for creating solid body tides.

Factory function for creating solid body tides, using a single real Love number at a single degree (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.). This function evaluates Eq. (6.6) from the IERS Conventions 2010, with real :math:`k_{l}=k_{lm}`, a single value of :math:`l` and a single tide-raising body :math:`j`.


Parameters
----------
tide_raising_body : str
    Name of body raising the tide.
love_number : float
    Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
degree : int
    Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
Returns
-------
BasicSolidBodyGravityFieldVariationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class





Examples
--------
In this example, we create gravity field variations of Earth for a tide raised by the Moon, with a single Love number :math:`k_{2}` of 0.301, and add it to the list of gravity field variations

.. code-block:: python 
   
   tide_raising_body = "Moon" 
   degree = 2 
   love_number = 0.301 
   gravity_field_variation_list = list() 
   gravity_field_variation_list.append( environment_setup.gravity_field_variation.solid_body_tide( 
   	tide_raising_body, love_number, degree ) 
   body_settings.get( "Earth" ).gravity_field_variation_settings = gravity_field_variation_list 


    )";



    } else if(name == "solid_body_tide_complex_k" && variant==0) {
        return R"(
        
Factory function for creating solid body tides.

As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide`, but with complex value for the Love number.


Parameters
----------
tide_raising_body : str
    Name of body raising the tide.
love_number : complex
    Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
degree : int
    Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
Returns
-------
BasicSolidBodyGravityFieldVariationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class






    )";



    } else if(name == "solid_body_tide_degree_variable_k" && variant==0) {
        return R"(
        
Factory function for creating solid body tides.

Factory function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.). This output of this function is effectively identical to a list of outputs to :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide`, with differing degrees and associated Love numbers.  This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{l}=k_{lm}`, at a set of values of :math:`l` and a single tide-raising body :math:`j`.


Parameters
----------
tide_raising_body : str
    Name of body raising the tide.
love_number_per_degree : dict( int, float )
    Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself
Returns
-------
BasicSolidBodyGravityFieldVariationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class





Examples
--------
In this example, we create gravity field variations of Earth for a tide raised by the Moon, with a Love numbers :math:`k_{2}=0.301`, and :math:`k_{3}=0.09`, and add it to the list of gravity field variations

.. code-block:: python 
   
   tide_raising_body = "Moon"
   love_numbers = dict( ) 
   love_numbers[ 2 ] = 0.301 
   love_numbers[ 3 ] = 0.09 
   gravity_field_variation_list = list() 
   gravity_field_variation_list.append( environment_setup.gravity_field_variation.solid_body_tide_degree_variable_k( 
   	tide_raising_body, love_numbers ) 
   body_settings.get( "Earth" ).gravity_field_variation_settings = gravity_field_variation_list 


    )";



    } else if(name == "solid_body_tide_degree_variable_complex_k" && variant==0) {
        return R"(
        
Factory function for creating solid body tides.

As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_variable_k`, but with complex values for the Love numbers.


Parameters
----------
tide_raising_body : str
    Name of body raising the tide.
love_number_per_degree : dict( int, complex )
    Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself.
Returns
-------
BasicSolidBodyGravityFieldVariationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class






    )";



    } else if(name == "solid_body_tide_degree_order_variable_k" && variant==0) {
        return R"(
        
Factory function for creating solid body tides.

Factory function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees and orders (e.g. :math:`k_{20}`, :math:`k_{21}`, :math:`k_{22}`, :math:`k_{30}`, etc.).  This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{lm}`, at a set of values of :math:`l` and a single tide-raising body :math:`j`.


Parameters
----------
tide_raising_body : str
    Name of body raising the tide.
love_number_per_degree_and_order : dict( int, list( float ) )
    Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`... :math:`k_{ll}`.
Returns
-------
BasicSolidBodyGravityFieldVariationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class





Examples
--------
In this example, we create gravity field variations of the Moon, for a tide raised by Earth, with a Love numbers :math:`k_{20}=0.024615`, :math:`k_{21}=0.023915` and :math:`k_{21}=0.024852`, and add it to the list of gravity field variations

.. code-block:: python 
   
   tide_raising_body = "Earth"
   love_numbers = dict( ) 
   love_numbers[ 2 ] = list( ) 
   love_numbers[ 2 ].append( 0.024615 ) 
   love_numbers[ 2 ].append( 0	.023915 ) 
   love_numbers[ 2 ].append( 0.024852 ) 
   gravity_field_variation_list = list() 
   gravity_field_variation_list.append( environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k( 
   	tide_raising_body, love_numbers ) 
   body_settings.get( "Earth" ).gravity_field_variation_settings = gravity_field_variation_list 


    )";



    } else if(name == "solid_body_tide_degree_order_variable_complex_k" && variant==0) {
        return R"(
        
Factory function for creating solid body tides.

As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`, but with complex values for the Love number.


Parameters
----------
tide_raising_body : str
    Name of body raising the tide.
love_number_per_degree : dict( int, list( complex ) )
    Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
Returns
-------
BasicSolidBodyGravityFieldVariationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace ground_station {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "GroundStationSettings") {
         return R"(

        Base class for providing settings for the creation of a ground station.





     )";





    } else if(name == "GroundStationMotionSettings") {
         return R"(

        Base class for providing settings for the motion of a single ground station.

        Non-functional base class for settings for the motion of a single ground station
        Station motion settings requiring additional information must be defined using an object derived from this class.





     )";





    } else if(name == "LinearGroundStationMotionSettings") {
         return R"(

        Class for defining linear motion (in an Earth-fixed frame) in time of a ground station.

        `GroundStationMotionSettings` derived class for time-linear station motion




     )";





    } else if(name == "PiecewiseConstantGroundStationMotionSettings") {
         return R"(

        Class for defining piecewise-constant position (e.g. instantaneous change in position at given epochs) of a ground station.

        `GroundStationMotionSettings` derived class for piecewise-constant position of a ground station




     )";





    } else if(name == "CustomGroundStationMotionSettings") {
         return R"(

        Class for defining custom time-dependent motion of a ground station.

        `CustomGroundStationMotionSettings` derived class for custom time-dependent motion of a ground station




     )";






    } else if(name == "basic_station" && variant==0) {
        return R"(
        
Factory function for creating settings for a ground station

Factory function for creating settings for a ground station, defining only its name, body-fixed position, and (optionally) time-variations of its position


Parameters
----------
station_name : string
    Name (unique identifier) by which the station is to be known.
station_position_element_type : PositionElementTypes, default = cartesian_position
    Type of elements for ``station_nominal_position``. Choose between cartesian_position, spherical_position and geodetic_position
station_nominal_position : np.ndarray([3,1])
    Nominal position of the station in a body-fixed frame. Depending on the choice of ``station_position_element_type`` input, this vector must contain
    * Cartesian - :math:`[x,y,z]`, denoting :math:`x-`, :math:`y-` and :math:`z-` components of body-fixed position (w.r.t body-fixed frame origin, typically center of mass) * Spherical - :math:`[r,\phi',\theta]`, denoting distance from body-fixed frame origin (typically center of mass), latitude and longitude * Geodetic - :math:`[h,\phi,\theta]`, denoting the altitude w.r.t. the body shape model, geodetic latitude and longitude
station_motion_settings : list[ GroundStationMotionSettings ], default = None
    List of settings defining time-variations of the individual ground station
Returns
-------
GroundStationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationSettings` defining settings of the to be created ground station





Examples
--------
In this example, we create a station using geodetic coordinates at the approximate location of the city of Delft, and no motion settings:

.. code-block:: python  

  # Define the position of the ground station on Earth 
  station_altitude = 0.0 
  delft_latitude = np.deg2rad(52.00667) 
  delft_longitude = np.deg2rad(4.35556) 
  
  # Create ground station settings
  ground_station_settings = environment_setup.ground_station.basic_station( 
      "TrackingStation", 
       [station_altitude, delft_latitude, delft_longitude], 
       element_conversion.geodetic_position_type) 
  
  # Append station settings to existing (default is empty) list 
  		body_settings.get( "Earth" ).ground_station_settings.append( ground_station_settings ) 


    )";



    } else if(name == "dsn_stations" && variant==0) {
        return R"(
        
Factory function for creating settings for all DSN stations

Factory function for creating settings for all DSN stations, defined by nominal positions and linear velocities, as defined
by Cartesian elements in *DSN No. 810-005, 301, Rev. K*,  see `this link <https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301K.pdf>`_.
Note that calling these settings will use the Cartesian elements provided in this document (in ITRF93) and apply them to the Earth-fixed
station positions, regardless of the selected Earth rotation model.

Returns
-------
list[ GroundStationSettings ]
    List of settings to create DSN stations






    )";



    } else if(name == "linear_station_motion" && variant==0) {
        return R"(
        
Factory function for creating settings for a linear station motion

Factory function for creating settings for a linear station motion, implementing :math:`\Delta \mathbf{r}=\dot{\mathbf{r}}(t-t_{0})`.


Parameters
----------
linear_velocity : np.ndarray([3,1])
    Linear velocity :math:`\dot{\mathbf{r}}` of the station (in m/s)
reference_epoch : float, default = 0.0
    Reference epoch :math:`t_{0}`, in seconds since J2000 epoch
Returns
-------
GroundStationMotionSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.LinearGroundStationMotionSettings` class






    )";



    } else if(name == "piecewise_constant_station_motion" && variant==0) {
        return R"(
        
Factory function for creating settings for a piecewise constant ground station position variation

Factory function for creating settings for a piecewise constant ground station position. Using this model, the added station velocity in a body-fixed frame :math:`\dot{\mathbf{r}}` is
always zero, but its displacement :math:`\Delta\mathbf{r}` is set according to the input list, which contains a list of times and displacments :math:`[t_{i},\Delta\mathbf{r}_{i}]`. 
When the resulting model is queried at a given time :math:`t`, the nearest lower neighbour :math:`t_{i}` from this list is found, and the associated :math:`\Delta\mathbf{r}_{i}` is applied.


Parameters
----------
displacement_list : dict[float,np.ndarray([3,1])]
    Dictionary with the epochs :math:`t_{i}` as values, and the associated displacement :math:`\Delta\mathbf{r}_{i}` as value
Returns
-------
GroundStationMotionSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.PiecewiseConstantGroundStationMotionSettings` class






    )";



    } else if(name == "custom_station_motion" && variant==0) {
        return R"(
        
Factory function for creating settings for a custom ground station position variation

Factory function for creating settings for a custom ground station position. An arbitrary user-defined function of the signature :math:`\Delta\mathbf{r}=\Delta\mathbf{r}(t)` is provided and
applied to the station position


Parameters
----------
custom_displacement_function : dict[float,np.ndarray([3,1])]
    Function returning :math:`\Delta\mathbf{r}`, with the time :math:`t` as input.
Returns
-------
GroundStationMotionSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.GroundStationMotionSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ground_station.CustomGroundStationMotionSettings` class






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace radiation_pressure {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "KnockeTypeSurfacePropertyDistributionModel") {
         return R"(

        Enumeration of available 'Knocke-type' surface distribution coefficient sets (see :func:`~knocke_type_surface_property_distribution`)





     )";


    } else if(name == "KnockeTypeSurfacePropertyDistributionModel.albedo_knocke") {
         return R"(
Coefficients for Earth surface albedo model from Knocke (1989)
     )";


    } else if(name == "KnockeTypeSurfacePropertyDistributionModel.emissivity_knocke") {
         return R"(
Coefficients for Earth surface emissivity model from Knocke (1989)
     )";



    } else if(name == "SphericalHarmonicSurfacePropertyDistribution") {
         return R"(

        Enumeration of available spherical harmonic surface distribution models coefficient sets (see :func:`~spherical_harmonic_surface_property_distribution`)





     )";


    } else if(name == "SphericalHarmonicSurfacePropertyDistribution.albedo_dlam1") {
         return R"(
Coefficients for DLAM lunar albedo model by FLoberhgen (1999)
     )";




    } else if(name == "LuminosityModelSettings") {
         return R"(

        Base class for providing settings for body source luminosity settings, to be used (typically but not necesarilly) for defining the Sun's luminosity.





     )";





    } else if(name == "SurfacePropertyDistributionSettings") {
         return R"(

        Base class for providing settings for body surface property distribution settings, to be used (typically but not necesarilly) for defining surface distribution of albedo and emissivity of solar system bodies for calculations of albedo and planetary radiation pressure.Note that not all albedo/emissivity models require this type of distribution model





     )";





    } else if(name == "PanelRadiosityModelSettings") {
         return R"(

        Base class for providing settings for body panel radiosity models, to be used (typically but not necesarilly) for defining surface radiosoty of a panelled solar system body as a result of albedo and/or planetary radiation pressure





     )";





    } else if(name == "BodyPanelReflectionLawSettings") {
         return R"(

        Base class for providing settings for body panel relfection law models, to be used for defining spacecraft surface properties relevant for the compuation of radiation pressure acting on a macromodel.





     )";





    } else if(name == "RadiationSourceModelSettings") {
         return R"(

        Base class for providing settings for properties of a radiation source (e.g. Sun), to be used in the context of (for instance) calculation of radiation pressure on spacecraft





     )";





    } else if(name == "RadiationPressureTargetModelSettings") {
         return R"(

        Base class for providing settings for properties of a radiation target (e.g. spacecraft), to be used in the context of (for instance) calculation of radiation pressure on spacecraft





     )";






    } else if(name == "constant_luminosity" && variant==0) {
        return R"(
        
Factory function for creating constant radiation source luminosity settings.

Factory function for creating constant radiation source luminosity settings, defining the total
radiated power (in Watts) of a given source. With this function, the source luminosity is constant, 
and is assumed to emit radiation isotropically.


Parameters
----------
luminosity : float
    Constant source luminosity (in Watt)
Returns
-------
LuminosityModelSettings
    Object defining settings for source luminosity






    )";



    } else if(name == "cannonball_radiation_target" && variant==0) {
        return R"(
        
Factory function for cannonball radtiation target


Parameters
----------
reference_area : float
    Cross-sectional area of cannonball [:math:`m^{2}`]
radiation_pressure_coefficient : float
    Radiation pressure coefficient [-]
per_source_occulting_bodies : Dict[str, List[str]]
    Names of bodies to occult the source as seen from this target
Returns
-------
CannonballRadiationPressureTargetModelSettings
    Object defining settings for a cannonball radiation pressure target model






    )";



    } else if(name == "irradiance_based_constant_luminosity" && variant==0) {
        return R"(
        
Factory function for creating source luminosity settings based on the irradiance at a reference distance.

Factory function for creating source luminosity based on the irradiance at a reference distance. For instance,
one can provide the solar irradiance at 1 AU, and this will be translated to the Sun's luminosity. With this function,
the source luminosity is constant, and is assumed to emit radiation isotropically.


Parameters
----------
constant_irradiance : float
    Irradiance at reference distance from center of source (in :math:`W/m^{2}`)
reference_distance : float
    Distance from center of source at which the irradiance is defined
Returns
-------
LuminosityModelSettings
    Object defining settings for source luminosity






    )";



    } else if(name == "time_variable_luminosity" && variant==0) {
        return R"(
        
Factory function for creating time-variable radiation source luminosity settings.

Factory function for creating time-variable radiation source luminosity settings, defining the total
radiated power (in Watts) of a given source as a function of time. With this function, the source 
is assumed to emit radiation isotropically.


Parameters
----------
luminosity_function : Callable[[float], float]
    Function returning source luminosity (in Watt) as a function of time
Returns
-------
LuminosityModelSettings
    Object defining settings for source luminosity






    )";



    } else if(name == "irradiance_based_time_variable_luminosity" && variant==0) {
        return R"(
        
Factory function for creating time-variable source luminosity settings based on the irradiance at a reference distance.

Factory function for creating source time-variable luminosity based on the irradiance at a reference distance. For instance,
one can provide the solar irradiance at 1 AU as a function of time, and this will be translated to the Sun's luminosity.
With this function, the source is assumed to emit radiation isotropically.


Parameters
----------
irradiance_function : Callable[[float], float]
    Function returning irradiance at reference distance from center of source (in :math:`W/m^{2}`) as a function fo time
reference_distance : float
    Distance from center of source at which the irradiance is defined
Returns
-------
LuminosityModelSettings
    Object defining settings for source luminosity






    )";



    } else if(name == "constant_surface_property_distribution" && variant==0) {
        return R"(
        
Factory function for creating constant radiative surface property distribution settings.

Factory function for creating constant radiative surface property (e.g. albedo, emmisivitiy, etc.) distribution settings.


Parameters
----------
constant_value : float
    Constant surface property value
Returns
-------
SurfacePropertyDistributionSettings
    Object defining settings for surface property distribution






    )";



    } else if(name == "spherical_harmonic_surface_property_distribution" && variant==0) {
        return R"(
        
Factory function for creating radiative surface property distribution settings according to a spherical harmonic model.

Factory function for creating radiative surface property (e.g. albedo, emmisivitiy, etc.) distribution settings 
according to a spherical harmonic model. The user provides unnormalized cosine and sine coefficients :math:`C_{lm}` and :math:`S_{lm}`,
from which the surface property :math:`k` is computed from:

.. math::
   k(\phi,\theta)=\sum_{l=0}^{l_{max}}\sum_{m=0}^{l}{P}_{lm}(\sin\phi)\left({C}_{lm}\cos m\theta+{S}_{lm}\sin m\theta\right)

with the angles :math:`\phi` and :math:`\theta` the body-fixed latitude and longitude of the evaluation point.


Parameters
----------
cosine_coefficients : numpy.ndarray
    Cosine coefficients of surface distribution. Entry (i,j) denotes coefficient :math:`{C}_{ij}` at degree i and order j.
sine_coefficients : numpy.ndarray
    Sine coefficients of surface distribution. Entry (i,j) denotes coefficient :math:`{C}_{ij}` at degree i and order j.
Returns
-------
SurfacePropertyDistributionSettings
    Object defining settings for surface property distribution






    )";



    } else if(name == "predefined_spherical_harmonic_surface_property_distribution" && variant==0) {
        return R"(
        
Factory function for creating radiative surface property distribution settings according to a predefined spherical harmonic model.

As :func:`spherical_harmonic_surface_property_distribution`, but with a predefined spherical harmonic distribution.


Parameters
----------
predefined_model : SphericalHarmonicsSurfacePropertyDistributionModel
    Identifier for predefined spherical harmonic surface property model.
Returns
-------
SurfacePropertyDistributionSettings
    Object defining settings for surface property distribution






    )";



    } else if(name == "knocke_type_surface_property_distribution" && variant==0) {
        return R"(
        
Factory function for creating radiative surface property distribution settings according to 'Knocke-type' model

Factory function for creating radiative surface property (e.g. albedo, emmisivitiy, etc.) distribution settings 
according to a model such as the one used by Knocke (1988). This model uses a degree two zonal spherical harmonic model, with 
a sinusoidal variation in the degree one coefficient. The surface property :math:`k` is computed from:

.. math::
   k(\phi,t)=a_{0}+a_{1}P_{1}(\sin\phi)+a_{2}P_{2}(\sin\phi)
.. math::
   a_{1}=c_{0}+c_{1}\cos\left(\frac{2\pi(t-t_{0})}{T}\right)+c_{2}\sin\left(\frac{2\pi(t-t_{0})}{T}\right)

with the angle :math:`\phi` denotes the body-fixed latitude of the evaluation point, and :math:`t`, :math:`t_{0}` and :math:`T` define the current time, 
reference time and period of the variation, respectively. The coefficients :math:`a_{0}, a_{2}, c_{0}, c_{1}, c_{2}` are provided by the user.


Parameters
----------
constant_contribution : float
    Value of :math:`a_{0}` in above formulation.
constant_degree_one_contribution : float
    Value of :math:`c_{0}` in above formulation.
cosine_periodic_degree_one_contribution : float
    Value of :math:`c_{1}` in above formulation.
sine_periodic_degree_one_contribution : float
    Value of :math:`c_{2}` in above formulation.
constant_degree_two_contribution : float
    Value of :math:`a_{2}` in above formulation.
reference_epoch : float
    Reference epoch :math:`t_{0}` of the periodic variation.
period : float
    Period :math:`T` of the periodic variation.
Returns
-------
SurfacePropertyDistributionSettings
    Object defining settings for surface property distribution






    )";



    } else if(name == "predefined_knocke_type_surface_property_distribution" && variant==0) {
        return R"(
        
Factory function for creating radiative surface property distribution settings according to a predefined 'Knocke-type` model.

As :func:`spherical_harmonic_surface_property_distribution`, but with a predefined spherical harmonic distribution.


Parameters
----------
predefined_model : KnockeTypeSurfacePropertyDistributionModel
    Identifier for predefined Knocke-type surface property model.
Returns
-------
SurfacePropertyDistributionSettings
    Object defining settings for surface property distribution






    )";



    } else if(name == "custom_surface_property_distribution" && variant==0) {
        return R"(
        
Factory function for creating radiative surface property distribution settings according to a custom user-defined model.

Factory function for creating radiative surface property (e.g. albedo, emmisivitiy, etc.) distribution settings
according to a custom user-defined model, as a function of latitude, longitude and time.


Parameters
----------
custom_function : Callable[[float, float, float], float]
    Function providing surface property as a function of latitude, longitude and time (in that order).
Returns
-------
SurfacePropertyDistributionSettings
    Object defining settings for surface property distribution






    )";



    } else if(name == "constant_radiosity" && variant==0) {
        return R"(
        
Factory function for creating settings for surface constant surface radiosity of an extended source

Factory function for creating settings for surface radiosity of an extended source, using constant Lambertian radiosity :math:`J` (in :math:`W/m^{2}`).
For a surface panel normal of :math:`\hat{\mathbf{n}}` and a vector :math:`\mathbf{r}` from the surface element to the target, the resulting
irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if :math:`\theta>0`, or in other words if the panel is visible from the target):

.. math::
   \Phi=J\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

with :math:`A` the panel area, :math:`\theta` is the angle between :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.      


Parameters
----------
radiosity : float
    Constant Lambertian radiosity from surface in :math:`W/m^{2}`.
Returns
-------
PanelRadiosityModelSettings
    Object defining settings for source panel radiosity






    )";



    } else if(name == "constant_albedo_surface_radiosity" && variant==0) {
        return R"(
        
Factory function for creating settings for surface constant albedo surface radiosity of an extended source

Factory function for creating settings for surface radiosity of an extended source, with surface radiation the result
of albedo using a Lambertian scattering law, and a constant albedo value over the surface.
For a surface panel normal of :math:`\hat{\mathbf{n}}`, a vector :math:`\mathbf{r}` from the surface element to the target, and a
vector :math:`\mathbf{r}_{s}` from the surface element to the original source (typically the Sun),
the resulting irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):

.. math::
   \Phi=\cos\theta_{s}\Phi_{s}\frac{a}{\pi}\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

with :math:`\theta_{s}` the angle between :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r_{s}}`, :math:`\Phi_{s}` the irradiance from the original source
at the panel of the reflecting body, :math:`a` is the albedo coefficient, :math:`A` the panel area, :math:`\theta` is the angle between :math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.


Parameters
----------
constant_albedo : float
    Constant value of the albedo coefficient :math:`a`.
original_source_name : str
    Name of the original source from which the radiation is reflection to the target.
Returns
-------
PanelRadiosityModelSettings
    Object defining settings for source panel radiosity






    )";



    } else if(name == "variable_albedo_surface_radiosity" && variant==0) {
        return R"(
        
Factory function for creating settings for surface variable albedo surface radiosity of an extended source

As :func:`constant_albedo_surface_radiosity`, but with the surface albedo :math:`a` defined by a surface distribution model.


Parameters
----------
albedo_distribution_settings : SurfacePropertyDistributionSettings
    Model for the surface distribution of the albedo :math:`a`.
original_source_name : str
    Name of the original source from which the radiation is reflection to the target.
Returns
-------
PanelRadiosityModelSettings
    Object defining settings for source panel radiosity






    )";



    } else if(name == "thermal_emission_blackbody_constant_emissivity" && variant==0) {
        return R"(
        
Factory function for creating settings for surface radiosity of an extended source from an isotropically heated body with constant emmisivity

Factory function for creating settings for surface radiosity of an extended source from an isotropically heated body (e.g. IR radiation) with constant surface
emissivity,
where the emitted power of the body is computed from the assumption that all heat absorbed from an original source is
emitted isotropically by the body. For instance, for Earth with Sun as original source, this model is equivalent to
assuming that a given fraction of all heat incident of the Sun on the Earth is absorbed and causes the full Earth surface to
heat to a constant temperature, which then results in the body emitting infrared radiation from its surface.

For a surface panel normal of :math:`\hat{\mathbf{n}}`, a vector :math:`\mathbf{r}` from the surface element to the target,
the resulting irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):

.. math::
   \Phi=\frac{\epsilon\Phi_{s}}{4}\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

with :math:`\epsilon` the emissivity, :math:`\Phi_{s}` the irradiance from the original source,  :math:`A` the panel area, :math:`\theta` is the angle between
:math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.


Parameters
----------
constant_emissivity : float
    Constant emissivity of the surface :math:`\epsilon`.
original_source_name : str
    Name of the original source from which the radiation is reflection to the target.
Returns
-------
PanelRadiosityModelSettings
    Object defining settings for source panel radiosity






    )";



    } else if(name == "thermal_emission_blackbody_variable_emissivity" && variant==0) {
        return R"(
        
Factory function for creating settings for surface radiosity of an extended source from an isotropically heated body with variable emmisivity

As :func:`thermal_emission_blackbody_constant_emissivity`, but with the surface emmisivity :math:`\epsilon` defined by a surface distribution model.


Parameters
----------
emissivity_distribution_model : SurfacePropertyDistributionSettings
    Model for the surface distribution of the emissivity :math:`\epsilon`.
original_source_name : str
    Name of the original source from which the radiation is reflection to the target.
Returns
-------
PanelRadiosityModelSettings
    Object defining settings for source panel radiosity






    )";



    } else if(name == "thermal_emission_angle_based_radiosity" && variant==0) {
        return R"(
        
Factory function for creating settings for surface radiosity of an extended source with surface temperature from Lemoine (2013)

Factory function for creating settings for surface radiosity of an extended source from an isotropically heated body (e.g. IR radiation)
with surface temperature :math:`T` computed from the angle of the surface normal and the original source as follows:

.. math::
   T=\max\left(T_{max}(\cos\phi_{s})^{1/4},T_{min} \right)

with :math:`phi_{s}` the angle along a great cirlce arc from the panel to the subsolar (for the Sun as original source) point; for
a circular body equivalent to the angle of the vector to the original source and the surface normal. The minimum and
maximum temperatures are user parameters.

For a surface panel normal of :math:`\hat{\mathbf{n}}`, a vector :math:`\mathbf{r}` from the surface element to the target,
the resulting irradiance :math:`\Phi` (in :math:`W/m^{2}`) at the target is (if the panel is visible from the target and the original source):

.. math::
   \Phi=\epsilon kT^{4}\frac{A\cos\theta}{\pi ||\mathbf{r}||^{2}}

with :math:`\epsilon` the emissivity, :math:`k` the Stefan-Boltzmann constant, :math:`A` the panel area, :math:`\theta` is the angle between
:math:`\hat{\mathbf{n}}` and :math:`\mathbf{r}`.


Parameters
----------
minimum_temperature : float
    Minimum surface temperature :math:`T_{min}`.
maximum_temperature : float
    Maximum surface temperature :math:`T_{min}`.
constant_emissivity : float
    Constant emissivity of the surface :math:`\epsilon`.
original_source_name : str
    Name of the original source from which the radiation is reflection to the target.
Returns
-------
PanelRadiosityModelSettings
    Object defining settings for source panel radiosity






    )";



    } else if(name == "specular_diffuse_body_panel_reflection" && variant==0) {
        return R"(
        
Factory function for creating settings for target panel reflection law using a specular-diffuse model

Factory function for creating settings for target panel reflection law used for a radiation pressure target, with a
specular diffuse model. The details of the implementation are given by Montenbruck et al. (2015). The reflection
law is defined by the absorption coefficient :math:`\alpha`, diffuse reflectivity :math:`\delta` and specular reflectivity
:math:`\rho`, which must meet the condition :math:`\alpha+\delta+\rho=1`. For the model definition, the user provides
:math:`\alpha` and :math:`\delta` (and :math:`\rho` is calculated). The reaction vector :math:`\hat{\mathbf{f}}` for a panel
with surface normal :math:`\hat{\mathbf{n}}`, and unit vector from panel surface to source :math:`\hat{\mathbf{r}}` then becomes:

.. math::
   \hat{\mathbf{f}}=\cos\theta\left((\alpha+\delta)\hat{\mathbf{r}}+(\frac{2}{3}\delta+2\rho\cos\theta)\hat{\mathbf{n}} \right)

In addition, it can be specified whether the absorbed radiation is also instantaneously retransmitted (according to Lambert's law), in
which case the above is modified to become:

.. math::
   \hat{\mathbf{f}}=\cos\theta\left((\alpha+\delta)\left(\hat{\mathbf{r}}+\frac{2}{3}\hat{\mathbf{n}}\right)+2\rho\cos\theta\hat{\mathbf{n}} \right)


Parameters
----------
specular_reflectivity : float
    Specular reflectivity :math:`\rho`.
diffuse_reflectivity : float
    Diffuse reflectivity :math:`\delta`.
with_instantaneous_reradiation : bool
    Boolean denoting whether absorbed radiation is instantaneously retransmitted (yes, if true).
Returns
-------
BodyPanelReflectionLawSettings
    Object defining settings for target panel reflection law






    )";



    } else if(name == "lambertian_body_panel_reflection" && variant==0) {
        return R"(
        
Factory function for creating settings for target panel reflection law using a Lambertian model

Factory function for creating settings for target panel reflection law used for a radiation pressure target, with a
purely Lambertian model. The implementation is as :func:`specular_diffuse_body_panel_reflection`, with
:math:`\rho=0` and no instantaneous reradiation. The only free parameter is the reflectivity :math:`\delta`, such that 
:math:`\alpha=1-\delta`.


Parameters
----------
reflectivity : float
    Reflectivity :math:`\delta`
Returns
-------
BodyPanelReflectionLawSettings
    Object defining settings for target panel reflection law






    )";



    } else if(name == "isotropic_radiation_source" && variant==0) {
        return R"(
        
Factory function for creating settings for an isotropic radiation source

Factory function for creating settings for a radiation source that emits isotropically. The source is provided
with a luminosity model :math:`L(t)` as a (possible) function of time :math:`t`. The irradiance :math:`\Phi` at a relative position
:math:`\mathbf{r}` from the source's center is then computed from:

.. math::
   \Phi=\frac{L}{4\pi||\mathbf{r}||^{2}}


Parameters
----------
luminosity_model : LuminosityModelSettings
    Settings for the luminosity model.
Returns
-------
RadiationSourceModelSettings
    Object defining settings for source model irradiance






    )";



    } else if(name == "panelled_radiation_target" && variant==0) {
        return R"(
        
Factory function for creating settings for a paneled radiation pressure target model

Factory function for creating settings for a paneled radiation pressure target model. Each source can have
its own set of occulting bodies.


Parameters
----------
source_to_target_occulting_bodies : Dict[str, List[str]]
    Map (source name -> list of occulting body names) of bodies to occult sources as seen from this target.
Returns
-------
RadiationPressureTargetModelSettings
    Object defining settings for a radiation pressure target model






    )";



    } else if(name == "panelled_extended_radiation_source" && variant==0) {
        return R"(
        
Factory function for creating settings for a dynamically panelled extended radiation source

Factory function for creating settings for a radiation source that is modelled as an anisotropic extended source,
such as a source due to albedo or planetary IR. The model can combined any number of superimposed surface panel radiosity models
(e.g. albedo, direct radiation), each of which may or may not involve an 'original source' (e.g. the Sun for albedo).
Each time the radiation at a given target is computed, the surface of the body is re-panelled, using the algorithm described by
Knocke (1989). In short, a single panel is placed at the zenith of the evaluation point, with any number of rings around it, each of
which has any number of (equispaced) panels on it. The width of each ring is such that all panels have the same projected, attenuated area. 
The panelling settings are defined by the user to this function. The
The irradiance :math:`\Phi` at a relative position :math:`\mathbf{r}` from the source's center is then computed from all 
:math:`N` panels, each of which evaluated :math:`M` panel radiosity models

.. math::
   \Phi=\sum_{i=1}^{N}\sum_{j=1}\Phi_{i,j}

where :math:`\Phi_{i,j}` denotes the contribution to the total irradiance of panel radiosity model :math:`j` on panel :math:`i`.


Parameters
----------
luminosity_model : LuminosityModelSettings
    Settings for the luminosity model.
Returns
-------
RadiationSourceModelSettings
    Object defining settings for source model irradiance






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace rigid_body {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "RigidBodyPropertiesSettings") {
         return R"(

        Base class for providing settings for rigid body model creation.

        This class is a functional base class for settings of gravity field models that require no information in addition to their type.
        Gravity field model classes requiring additional information must be created using an object derived from this class.





     )";


    } else if(name == "RigidBodyPropertiesSettings.body_mass_property_type") {
         return R"(

        **read-only**

        Type of rigid body model that is to be created.

        :type: RigidBodyPropertiesType
     )";






    } else if(name == "constant_rigid_body_properties" && variant==0) {
        return R"(
        
Factory function for creating constant rigid body properties.

Factory function for creating constant rigid body properties (mass, center of mass, inertia tensor). The center of mass and/or inertia tensor can be left empty by setting them
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
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






    )";



    } else if(name == "custom_time_dependent_rigid_body_properties" && variant==0) {
        return R"(
        
Factory function for creating custom (time-dependent) rigid body properties.

Factory function for creating custom rigid body properties, where the mass, center of mass and inertia tensor are defined by user-defined functions (as a function of time).
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
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






    )";



    } else if(name == "custom_mass_dependent_rigid_body_properties" && variant==0) {
        return R"(
        
Factory function for creating custom (time-dependent) rigid body properties.

Factory function for creating custom rigid body properties, center of mass and inertia tensor are defined by user-defined functions as a function of mass.
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
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace rotation_model {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "RotationModelType") {
         return R"(

        Enumeration of rotation model types.

        Enumeration of rotation model types supported by tudat.





     )";


    } else if(name == "RotationModelType.simple_rotation_model") {
         return R"(
     )";


    } else if(name == "RotationModelType.spice_rotation_model") {
         return R"(
     )";


    } else if(name == "RotationModelType.gcrs_to_itrs_rotation_model") {
         return R"(
     )";


    } else if(name == "RotationModelType.synchronous_rotation_model") {
         return R"(
     )";


    } else if(name == "RotationModelType.planetary_rotation_model") {
         return R"(
     )";



    } else if(name == "IAUConventions") {
         return R"(

        Enumeration of IAU conventions for Earth rotation.

        Enumeration of IAU conventions for Earth rotation supported by tudat.





     )";


    } else if(name == "IAUConventions.iau_2000_a") {
         return R"(
     )";


    } else if(name == "IAUConventions.iau_2000_b") {
         return R"(
     )";


    } else if(name == "IAUConventions.iau_2006") {
         return R"(
     )";




    } else if(name == "RotationModelSettings") {
         return R"(

        Base class for providing settings for automatic rotation model creation.

        This class is a functional base class for settings of rotation models that require no information in addition to their type.
        Basic rotation model has constant orientation of the rotation axis (body-fixed z-axis) and constant rotation rate about this axis.
        Rotation models requiring additional information must be created using the factory functions which create the specific object derived from this base class.





     )";


    } else if(name == "RotationModelSettings.rotation_type") {
         return R"(

        **read-only**

        Type of rotation model that is to be created.

        :type: RotationModelType
     )";


    } else if(name == "RotationModelSettings.base_frame") {
         return R"(

        Name of the base frame of rotation model.

        :type: str
     )";


    } else if(name == "RotationModelSettings.target_frame") {
         return R"(

        **read-only**

        Name of the target frame of rotation model.

        :type: str
     )";






    } else if(name == "simple" && variant==0) {
        return R"(
        
Factory function for creating simple rotation model settings.

Factory function for settings object, defining a basic rotation model with constant orientation of the rotation axis and constant rotation rate about this axis.
Rotation from original (inertial) to target (body-fixed) frame at some reference time ``initial_time`` (:math:`t_{0}`) is defined by the ``initial_orientation`` (:math:`\mathbf{R}^{(B/I)}(t_{0})`) rotation matrix.
Rotation about the body-fixed z-axis is defined by the ``rotation_rate`` (:math:`\omega`) float variable (in rad/s). The rotation matrix is computed from:

.. math::
   \mathbf{R}^{(B/I)}(t)=\mathbf{R}_{z}(\omega(t-t_{0}))(t_{0})\mathbf{R}^{(B/I)}(t_{0})

where :math:`\mathbf{R}^{(B/I)}` denotes the rotation matrix from inertial to body-fixed frame, and :math:`\mathbf{R}_{z}` denotes a rotaion matrix about the z-axis.

The matrix :math:`\mathbf{R}^{(B/I)}(t_{0})` is sometimes parameterized by pole right ascension and declination (:math:`\alpha` and :math:`\delta`), as well as the meridian of date :math:`W_{0}` with

.. math::
   \mathbf{R}^{(B/I)}(t_{0})=\mathbf{R}_{z}(W_{0})\mathbf{R}_{x}(\pi/2-\delta)\mathbf{R}_{z}(\pi/2+\alpha)


Parameters
----------
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
initial_orientation : numpy.ndarray[numpy.float64[3, 3]]
    Orientation of target frame in base frame at initial time.
initial_time : float
    Initial time (reference epoch for rotation matrices).
rotation_rate : float
    Constant rotation rate [rad/s] about rotational axis.
Returns
-------
SimpleRotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` for Earth,
using a simple rotation model with constant orientation of the rotation axis (body-fixed z-axis), and constant rotation rate about this axis:

.. code-block:: python 
  
  # Set parameters describing the rotation between the two frames 
  initial_orientation = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]]) 
  initial_time = 12345 # [sec since J2000] 
  rotation_rate = 2e-5 # [rad/s] 
  original_frame = "J2000" 
  target_frame = "Earth_Fixed_Simplified" 
  # Create the rotation model settings and assign to body settings of "Earth" 
  body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.simple( 
      original_frame, 
      target_frame, 
      initial_orientation, 
      initial_time, 
      rotation_rate) 


    )";



    } else if(name == "simple_from_spice" && variant==0) {
        return R"(
        
Factory function for creating simple rotation model settings using initial orientation and rotation rates from Spice.

Factory function for settings object, defining a :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` rotation model with the added functionality that the initial orientation and rotation rate are extracted from Spice, as opposed to provided manually.
Note that `only` the initial orientation and rotation rate ( at the time defined by ``initial_time`` ) are extracted from Spice - for
the full Spice rotation model see :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.spice`.
Also note the distinction between the ``target_frame`` and ``target_frame_spice`` parameters.


Parameters
----------
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Target frame of rotation model - name of frame that Tudat assigns to the body-fixed frame
target_frame_spice : str
    Spice reference of target frame - name of the frame in Spice for which the initial orientation and rotation rate are extracted.
initial_time : float
    Initial time (reference epoch for rotation matrices).
Returns
-------
SimpleRotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class



Notes
-----
In order to create a :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` object which describes a synchronous rotation w.r.t. some ``central_body``,
we require an ``ephemeris_settings`` attribute to the :class:`~tudatpy.numerical_simulation.environment_setup.BodySettings` object of the ``central_body``.



Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` for Earth,
using a simple rotation model with constant orientation of the rotation axis (body-fixed z-axis), and constant rotation rate about this axis.
The initial orientation and rotation rate are extracted from Spice at the time defined by ``initial_time``:

.. code-block:: python 
   
   # set parameters for time at which initial data is extracted from spice 
   initial_time = 12345 
   # set parameters for defining the rotation between frames 
   original_frame = "J2000" 
   target_frame = "IAU_Earth_Simplified" 
   target_frame_spice = "IAU_Earth" 
   # create rotation model settings and assign to body settings of "Earth" 
   body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.simple_from_spice( 
   original_frame, target_frame, target_frame_spice, initial_time) 


    )";



    } else if(name == "synchronous" && variant==0) {
        return R"(
        
Factory function for creating synchronous rotational ephemeris settings.

Factory function for settings object, defining a synchronous rotation model where rotation of a body is defined from its relative orbit w.r.t. some central body. Specifically
- the body-fixed x-axis is *always* pointing towards the central body
- the body-fixed z-axis is *always* perpendicular to the orbital plane (along the direction of :math:`\mathbf{x}\times\mathbf{v}` )
- the body-fixed y-axis completes the right-handed reference frame

Such a model can be useful for, for instance, approximate rotation of tidally locked natural satellites or nadir-pointing spacecraft.


Parameters
----------
central_body_name : str
    Name of the base frame of rotation model.
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Spice reference of target frame.
Returns
-------
SynchronousRotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SynchronousRotationModelSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` for the martian moon Phobos,
We do so by assigning a synchronous rotation model to the rotation model settings of Phobos, using in this case ``"ECLIPJ2000"`` as the base frame,
and ``"Phobos_Fixed"`` as the target frame.

.. code-block:: python 
   
   # define parameters describing the synchronous rotation model
   central_body_name = "Mars"
   original_frame = "ECLIPJ2000"
   target_frame = "Phobos_Fixed"
   # create rotation model settings for target frame and assign to body settings of "Phobos" 
   body_settings.get( "Phobos" ).rotation_model_settings = environment_setup.rotation_model.synchronous(
   central_body_name, original_frame, target_frame)


    )";



    } else if(name == "spice" && variant==0) {
        return R"(
        
Factory function for creating rotation model settings from the Spice interface.

Factory function for settings object, defining a rotation model directly (and entirely) from Spice interface.


Parameters
----------
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
spice_frame_name : str, default = ""
    Name of the spice reference frame name that will be used to compute the rotation to the target frame. For instance, if target_frame is set to "IAU_Earth", and ``spice_frame_name`` is set to "IAU_Mars", Tudat will extract the rotation to the IAU_Mars frame from Spice, and assign this rotation to the "IAU_Earth" frame in Tudat. By default, this input is left empty, which corresponds to it being equal to the ``target_frame``.
Returns
-------
RotationModelSettings
    Instance of :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` class.





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` for Earth,
using full rotation model data from Spice:

.. code-block:: python 
   
   # define parameters describing the rotation between frames 
   original_frame = "J2000" 
   target_frame = "IAU_Earth" 
   # create rotation model settings and assign to body settings of "Earth" 
   body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.spice( 
   original_frame, target_frame) 


    )";



    } else if(name == "gcrs_to_itrs" && variant==0) {
        return R"(
        
Factory function for creating high-accuracy Earth rotation model settings.

Factory function for settings object, defining high-accuracy Earth rotation model according to the IERS 2010 Conventions.
This settings class has various options to deviate from the default settings, typical applications will use default.
Note that for this model the original frame must be J2000 or GCRS (in the case of the former, the frame bias between GCRS and J2000 is automatically corrected for). The target frame (e.g. body-fixed frame) name is ITRS.
The precession-nutation theory may be any member of :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.IAUConventions` (``iau_2000a`` / ``iau_2000b`` or ``iau_2006``).
Alternative options to modify the input (not shown here) include the EOP correction file, input time scale, short period UT1 and polar motion variations.
The target frame (e.g. body-fixed frame) name is ITRS.


Parameters
----------
precession_nutation_theory : IAUConventions, default=tba::iau_2006
    Setting theory for modelling Earth nutation.
base_frame : str, default='GCRS'
    Base frame of rotation model
Returns
-------
GcrsToItrsRotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.GcrsToItrsRotationModelSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` for Earth,
using a high-accuracy Earth rotation model as defined by IERS 2010 conventions:


.. code-block:: python 
   
   # define parameters describing the rotation between frames 
   precession_nutation_theory = environment_setup.rotation_model.IAUConventions.iau_2006 
   original_frame = "J2000" 
   # create rotation model settings and assign to body settings of "Earth" 
   body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.gcrs_to_itrs( 
   precession_nutation_theory, original_frame) 


    )";



    } else if(name == "constant_rotation_model" && variant==0) {
        return R"(
        
Factory function for creating simple rotation model settings for target-frames with constant orientation.

Factory function for settings object, defining simple rotation model setting objects with constant rotation matrix.
These model settings are for target frames which do not have a rotational rate in the base frame and are fully defined by their initial orientation.


Parameters
----------
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
initial_orientation : numpy.ndarray[numpy.float64[3, 3]]
    Rotation matrix from inertial to body-fixed (base to target) frame at initial time (constant throughout).
Returns
-------
SimpleRotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class.





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` for Earth,
using a constant rotation matrix between Earth-fixed and inertial frame:

.. code-block:: python 
  
  # define parameters describing the constant orientation between frames 
  original_frame = "ECLIPJ2000"  
  target_frame = "Earth_fixed"  
  constant_orientation = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])  
  # create rotation model settings and assign to body settings of "Earth" 
  body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.constant( 
      original_frame, 
      target_frame, 
      constant_orientation ) 


    )";



    } else if(name == "aerodynamic_angle_based" && variant==0) {
        return R"(
        
Factory function for creating rotation model settings based on custom aerodynamic angles (attack, sideslip, bank).

Factory function for creating rotation model settings based on custom aerodynamic angles:
angle of attack :math:`\alpha`, sideslip angle :math:`\beta` and bank angle :math:`\sigma`. The use of this function is typical for
simulating the dynamics of a (guided) re-entry vehicle. It calculates the rotation matrix from inertial frame to the body-fixed frame
of the current body B (typically a vehicle) w.r.t. the body-fixed frame of a central body C (e.g., the body at which the re-entry is taking place.
The full algorithm for :math:`R^{(I/B)}` is described by Mooij (1994), and is composed of:

*  The rotation from inertial frame to the body fixed frame of body C, using the existing rotation model of body C
*  The rotation from body-fixed frame of body C to the vehicle's vertical frame V. This rotation uses the current latitude and longitude angles.
*  The rotation of the vehicle's vertical frame V to its trajectory frame T. This rotation uses the current heading and flight path angles.
*  The rotation of the vehicle's trajectory frame T to its aerodynamic frame A. This rotation uses the current bank angle
*  The rotation of the vehicle's aerodynamic frame A to its body-fixed frame. This rotation uses the current angle of attack and sideslip angles

In the above algorithm, the latitude, longitude, heading and flight-path angles are computed from the vehicle's current translational state, in the body-fixed
frame of body C. The angle of attack, sideslip angle and bank angle are to be defined by the user, through a single custom function that is passed to
the ``angle_function`` argument of this functions


Parameters
----------
central_body : str
    Name of the central body C that is to be used.
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
angle_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]], default = None
    Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\alpha,\beta,\sigma]`. If this input is left empty, these angles are both fixed to 0.
Returns
-------
CustomRotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model.






    )";



    } else if(name == "zero_pitch_moment_aerodynamic_angle_based" && variant==0) {
        return R"(
        
Factory function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank.

Factory function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank. This function is
largely identical to the :func:`~aerodynamic_angle_based`, with the difference that the angle of attack :math:`\alpha` is not provided as a custom value by the user, but is
calculated from the body's aerodynamic moment coefficients, such that we have :math:`C_{m}=0`. This requires aerodynamic moment coefficients to be defined for the vehicle that
depend on (among others) the body's angle of attack


Parameters
----------
central_body : str
    Name of the central body C that is to be used.
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
angle_funcion : Callable[[float], numpy.ndarray[numpy.float64[2, 1]]], default = None
    Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\beta,\sigma]`. If this input is left empty, these angles are both fixed to 0.
Returns
-------
CustomRotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model.






    )";



    } else if(name == "custom_inertial_direction_based" && variant==0) {
        return R"(
        
Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction

Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction :math:`\hat{\mathbf{T}}_{I}`. Specifically, it ensures
that the rotation matrix from body-fixed to inertial frame is set up such that :math:`\hat{\mathbf{T}}_{I}=R^{(I/B)}\hat{\mathbf{i}}` (where :math:`\mathbf{i}` is the unit-vector in local x-direction).
The complete rotation matrix requires an additional angle :math:`\phi` (rotation of the body about its body-fixed x-axis), which is set to 0 by default.

The full rotation matrix is computed from a 3-2-1 Euler angle rotation
:math:`R^{(I/B)}=R_{z}(\psi)R_{y}(\theta)R_{x}(\phi)`, with :math:`\psi` and :math:`\theta` computed from the suitable decomposition of :math:`\hat{\mathbf{T}}_{I}`.
This function is typically used for simulating the (guided) dynamics of a spacecraft under thrust, where the thrust is provided in the x-direction of the body-fixed frame. By providing a suitable
``inertial_body_axis_direction``, this thrust can be defined to point in an arbitrary direction (typically defined by a guidance algorithm) in the inertial frame as a function of time.

NOTE: this function may be extended in the future to allow an arbitrary body-fixed direction to align with an arbitrary inertial direction. At present, its functionality is limited to imposing the inertial direction of the body-fixed x-axis.


Parameters
----------
inertial_body_axis_direction : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
    Custom function defined by the user, which imposes the inertial orientation of the body-fixed x-axis, by providing :math:`\hat{\mathbf{T}}_{I}(t)`.
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
free_rotation_angle_function : Callable[[float], float], default = None
    Custom function provided by the user, which returns a value for the free rotation angle :math:`\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
Returns
-------
BodyFixedDirectionBasedRotationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model.






    )";



    } else if(name == "orbital_state_direction_based" && variant==0) {
        return R"(
        
Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector.

Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector. This function is
similar to the :func:`~custom_inertial_direction_based` function, with the exception that the :math:`\hat{\mathbf{T}}_{I}` vector is not defined by thee user, but is defined by the
relative position vector :math:`\mathbf{r}_{C}` or velocity vector :math:`\mathbf{r}_{C}` of the vehicle w.r.t. some body C. The inputs to this function allow :math:`\hat{\mathbf{T}}_{I}` to
be set to :math:`\pm\mathbf{r}_{C}` or :math:`\pm\mathbf{v}_{C}`, for any body C. It is typically used for simplified or preliminary thrust analyses.


Parameters
----------
central_body : str
    Name of central body w.r.t. which the position/velocity vector is to be computed
is_colinear_with_velocity : bool
    Boolean defining whether :math:`\hat{\mathbf{T}}_{I}` is to be aligned with velocity (if true) or position (if false)
direction_is_opposite_to_vector : bool
    Boolean defining whether :math:`\hat{\mathbf{T}}_{I}` is to be in the same direction as position/velocity (if false), or in the opposite direction (if true).
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
free_rotation_angle_function : Callable[[float], float], default = None
    Custom function provided by the user, which returns a value for the free rotation angle :math:`\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
Returns
-------
BodyFixedDirectionBasedRotationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model.






    )";



    } else if(name == "mars_high_accuracy" && variant==0) {
        return R"(
        
Factory function for creating a high-accuracy Mars rotation model.

Factory function for creating a high-accuracy Mars rotation model, using the default parameters of `Konopliv et al. (2016) <https://www.sciencedirect.com/science/article/abs/pii/S0019103516001305>`_
and the mathematical model of ` Konopliv et al. (2006) <https://www.sciencedirect.com/science/article/pii/S0019103506000297>`_. The rotation matrix formulation is given in Eq. (13)-(19) of that paper.
Note that, at the moment, all parameters in this rotation model are hard coded, and cannot be adapted by the user (except by estimating a number of its constituent parameters, see :ref:`\`\`parameter\`\`` module )
As such, this model is at present applicable to Mars rotation only. If you require more fine-grained control of the parameters, please contact the Tudat support team      


Parameters
----------
base_frame : str, default = "ECLIPJ2000"
    Name of the base frame of rotation model.
target_frame : str, default = "Mars_Fixed"
    Name of the target frame of rotation model.
Returns
-------
RotationModelSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.PlanetaryRotationModelSettings` class, which defines the required settings for the rotation model.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace shape {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "BodyShapeSettings") {
         return R"(

        Base class for providing settings for body shape model.

        Functional (base) class for settings of body shape models that require no information in addition to their type.
        Body shape model settings requiring additional information must be defined using an object derived from this class.





     )";





    } else if(name == "SphericalBodyShapeSettings") {
         return R"(

        Class for defining model settings of a strictly spherical body shape.

        `BodyShapeSettings` derived class for strictly spherical body shape model settings.




     )";


    } else if(name == "SphericalBodyShapeSettings.radius") {
         return R"(

        **read-only**

        Radius specifying spherical body shape.

        :type: float
     )";





    } else if(name == "OblateSphericalBodyShapeSettings") {
         return R"(

        Class for defining model settings of a oblate spherical body shape.

        `BodyShapeSettings` derived class for oblate spherical body shape model settings.




     )";


    } else if(name == "OblateSphericalBodyShapeSettings.equatorial_radius") {
         return R"(

        **read-only**

        Equatorial radius of the oblate spherical body shape.

        :type: float
     )";


    } else if(name == "OblateSphericalBodyShapeSettings.flattening") {
         return R"(

        **read-only**

        Flattening of spheroid shape model.

        :type: float
     )";





    } else if(name == "PolyhedronBodyShapeSettings") {
         return R"(

        Class for defining model settings of a polyhedron body shape.

        `BodyShapeSettings` derived class for polyhedron body shape model settings.




     )";


    } else if(name == "PolyhedronBodyShapeSettings.vertices_coordinates") {
         return R"(

        Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
        row per vertex, 3 columns). 


        :type: numpy.ndarray
     )";


    } else if(name == "PolyhedronBodyShapeSettings.vertices_defining_each_facet") {
         return R"(

        Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
        the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
        when seen from the outside of the polyhedron. 


        :type: numpy.ndarray
     )";


    } else if(name == "PolyhedronBodyShapeSettings.compute_altitude_with_sign") {
         return R"(

        Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
        having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is
        recommended, as it reduces the CPU time for computing the altitude. 


        :type: bool, default=True
     )";


    } else if(name == "PolyhedronBodyShapeSettings.just_compute_distance_to_vertices") {
         return R"(

        Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
        is set to *true*) or to all polyhedron features (vertices, facets and edges; happens if flag is set to
        *false*). Depending on the application, it might be useful to set the flag to *true* for medium to high
        altitudes, as it allows significantly reducing the CPU time (the resulting altitude errors depend on the
        resolution of the used polyhedron and altitude itself). 


        :type: bool, default=False
     )";





    } else if(name == "HybridBodyShapeSettings") {
         return R"(

        Class for defining model settings of a hybrid body shape.

        `BodyShapeSettings` derived class for hybrid body shape model settings.




     )";


    } else if(name == "HybridBodyShapeSettings.low_resolution_body_shape_settings") {
         return R"(

        Settings of the shape model that is to be used to compute the altitude at high altitudes (above the switchover
        altitude). 


        :type: BodyShapeSettings
     )";


    } else if(name == "HybridBodyShapeSettings.high_resolution_body_shape_settings") {
         return R"(

        Settings of the shape model that is to be used to compute the altitude at low altitudes (below the switchover
        altitude). 


        :type: BodyShapeSettings
     )";


    } else if(name == "HybridBodyShapeSettings.switchover_altitude") {
         return R"(

        Altitude at which the model used to compute the altitude is changed. The high-resolution model is used for
        altitudes below the switchover altitude, the low-resolution model for altitudes above it. 


        :type: float
     )";






    } else if(name == "spherical" && variant==0) {
        return R"(
        
Factory function for creating spherical body shape model settings.

Factory function for settings object, defining strictly spherical body shape model entirely from single radius parameter.


Parameters
----------
radius : float
    Radius specifying spherical body shape.
Returns
-------
SphericalBodyShapeSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.SphericalBodyShapeSettings` class





Examples
--------
In this example, we create a :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` using a perfectly spherical shape model:

.. code-block:: python 
   
   # define parameters describing perfectly spherical model 
   body_radius = 6378.0E3 
   # create shape model settings  
   body_settings.get( "Earth" ).shape_settings = environment_setup.shape.spherical( body_radius ) 


    )";



    } else if(name == "spherical_spice" && variant==0) {
        return R"(
        
Factory function for creating spherical body shape model settings entirely from spice.

Factory function for settings object, defining spherical body shape model entirely from spice parameters.

Returns
-------
BodyShapeSettings
    Instance of :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` class





Examples
--------
In this example, we create a :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` using a perfectly spherical shape model and data from Spice:

.. code-block:: python 
   
   # create shape model settings  
   body_settings.get( "Earth" ).shape_settings = environment_setup.shape.spherical_spice( )


    )";



    } else if(name == "oblate_spherical" && variant==0) {
        return R"(
        
Factory function for creating oblate spherical body shape model settings.

Factory function for settings object, defining oblate spherical body shape model from equatorial radius and flattening parameter.


Parameters
----------
equatorial_radius : float
    Equatorial radius specifying oblate spherical body shape.
flattening : float
    Flattening parameter specifying oblate spherical body shape.
Returns
-------
OblateSphericalBodyShapeSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.OblateSphericalBodyShapeSettings` class





Examples
--------
In this example, we create a :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` using a perfectly oblate spherical shape model:

.. code-block:: python 
   
   # define parameters describing oblate spherical model 
   body_radius = 6378.0E3 
   body_flattening = 1.0 / 300.0 
   # create shape model settings  
   body_settings.get( "Earth" ).shape_settings = environment_setup.shape.oblate_spherical( body_radius, body_flattening ) 


    )";



    } else if(name == "polyhedron" && variant==0) {
        return R"(
        
Factory function for creating a polyhedron body shape model settings.

Factory function for settings object, defining a polyhedron shape model.

Note 1: The evaluation of the altitude with a polyhedron model tends to be computationally expensive. To reduce the
computational time, it might be useful to instead define a hybrid shape model (see
:func:`~tudatpy.numerical_simulation.environment_setup.shape.hybrid`), which allows using a high-resolution
polyhedron (with a large number of facets) at low altitudes and a low-resolution one (with smaller number of facets)
at high-altitudes.

Note 2: If the goal of using the shape model is only to detect collisions with the surface and not to explicitly
obtain the altitude, it is instead recommended to use the Laplacian of the gravitational potential (see
:func:`~tudatpy.numerical_simulation.propagation_setup.dependent_variable.gravity_field_laplacian_of_potential`).
This allows reducing the computational time, but is only valid if the same polyhedron model that is used to define
the gravitational acceleration should also be used to detect the impacts.


Parameters
----------
vertices_coordinates : numpy.ndarray
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
    row per vertex, 3 columns).

vertices_defining_each_facet : numpy.ndarray
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
    when seen from the outside of the polyhedron.

compute_altitude_with_sign : bool, default=True
    Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
    having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is
    recommended, as it reduces the CPU time.

just_compute_distance_to_vertices : bool, default=False
    Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
    is set to *true*) or to all polyhedron features (vertices, facets and edges; happens if flag is set to
    *false*). Depending on the application, it might be useful to set the flag to *true* for medium to high
    altitudes, as it allows significantly reducing the CPU time (the resulting altitude errors depend on the
    resolution of the used polyhedron and altitude itself).

Returns
-------
PolyhedronBodyShapeSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived
    :class:`~tudatpy.numerical_simulation.environment_setup.shape.PolyhedronBodyShapeSettings` class







    )";



    } else if(name == "hybrid" && variant==0) {
        return R"(
        
Factory function for creating hybrid body shape model settings.

Factory function for settings object, defining a hybrid shape model.

The hybrid shape model is constituded by two shape models: a low-resolution model which is used at high altitudes
(above the switchover altitude) and a high-resolution model used at low altitudes (below the switchover altitude).
In each computation of the altitude, the altitude is first computed with the low-resolution model. The
low-resolution altitude is then compared to the switchover altitude to decide whether to compute the high-resolution
altitude.

The hybrid shape model is useful when the evaluation of the high-resolution model is computationally expensive
(e.g. polyhedron model).


Parameters
----------
low_resolution_body_shape_settings : BodyShapeSettings
    Settings of the shape model that is to be used to compute the altitude at high altitudes (above the switchover
    altitude).

high_resolution_body_shape_settings : BodyShapeSettings
    Settings of the shape model that is to be used to compute the altitude at low altitudes (below the switchover
    altitude).

switchover_altitude : float
    Altitude at which the model used to compute the altitude is changed. The high-resolution model is used for
    altitudes below the switchover altitude, the low-resolution model for altitudes above it.

Returns
-------
HybridBodyShapeSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived
    :class:`~tudatpy.numerical_simulation.environment_setup.shape.HybridBodyShapeSettings` class







    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace shape_deformation {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "BodyDeformationSettings") {
         return R"(

        Base class for providing settings for body shape deformation model.

        Functional (base) class for settings of body shape deformation models that require no information in addition to their type.
        Body shape deformation model settings requiring additional information must be defined using an object derived from this class.





     )";





    } else if(name == "BasicSolidBodyDeformationSettings") {
         return R"(

        Class for defining model settings for simple tidal solid-body shape deformation.

        `BodyDeformationSettings` derived class for simple tidal solid-body shape deformation.




     )";


    } else if(name == "BasicSolidBodyDeformationSettings.radius") {
         return R"(

        **read-only**

        Radius specifying spherical body shape.

        :type: float
     )";






    } else if(name == "basic_solid_body_tidal" && variant==0) {
        return R"(
        
Factory function for creating basic tidal solid-body shape deformation

Factory function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{m}` and :math:`l_{m}` (with only :math:`m=2,3` currently supported). This function implements equations (7.5) and (7.6) of the `IERS 2010 Conventions <https://iers-conventions.obspm.fr/conventions_material.php>`_.


Parameters
----------
tide_raising_bodies : list[ string ]
    List of bodies that raise a tide that induces the shape variation.
displacement_love_numbers : dict[ int, [float,float] ]
    Dictionary of pairs. The dictionary key the spherical harmonic degree :math:`l` of the tidal deformation (2 or 3 are currenty supported). The dictionary value is comprised of a pair of floats representing the :math:`h_{2}` and :math:`l_{2}` deformation Love numbers
reference_radius : float, default = NaN
    Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
Returns
-------
BasicSolidBodyDeformationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class





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
          tide_raising_bodies, love_numbers ) )


    )";



    } else if(name == "degree_two_basic_solid_body_tidal" && variant==0) {
        return R"(
        
Factory function for creating degree 2 basic tidal solid-body shape deformation

Factory function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{2}` and :math:`l_{2}`. This function implements equations (7.5) of the IERS 2010 Conventions, and provides a simplified interface (for degree 2 only) of :func:`basic_solid_body_tidal`.


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
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class





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
          tide_raising_bodies, h2_love_number, l2_shida_number ) )


    )";



    } else if(name == "iers_2010_solid_body_tidal" && variant==0) {
        return R"(
        
Factory function for creating full IERS 2010 shape deformation model

Factory function for creating full IERS 2010 shape deformation model, computing the tidal shape variation due to the full model defined in Section 7.1.1 of the 2010 IERS conventions, implementing Eqs. (7.5)-(7.13), including all terms from Tables 7.3a and 7.3b. At present, none of the input parameters of the model can be varied.

Returns
-------
BodyDeformationSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape_deformation.BodyDeformationSettings` defining the IERS 2010 settings






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace vehicle_systems {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "BodyPanelGeometrySettings") {
         return R"(

        Base class for defining the geometrical properties of a single panel on the vehicle's exterior





     )";





    } else if(name == "BodyPanelSettings") {
         return R"(

        Class for defining the complete properties of a single panel on the vehicle's exterior





     )";





    } else if(name == "FullPanelledBodySettings") {
         return R"(

        Class for providing the complete settings for a panelled body exterior





     )";






    } else if(name == "frame_fixed_panel_geometry" && variant==0) {
        return R"(
        
Factory function for creating settings for a vehicle exterior panel that is fixed to a given frame.

Factory function for creating settings for a vehicle exterior panel that is fixed to a given frame, meaning
that the orientation of the panel is fully defined by the rotation model(s) defined in the vehicle.
The constant surface normal :math:`\hat{\mathbf{n}}^{\mathcal{F}}` in frame :math:`\mathcal{F}` is provided by the user.      
If the ``frame_orientation`` of this function is left empty, the panel is fixed to the body-frame, and 
:math:`\mathcal{F}` is the  body-fixed frame :math:`\mathcal{B}`.

Alternatively, the ``frame_orientation`` may be defined as the identifier of the frame fixed to one of the
vehicle parts (solar array, antenna, etc.). See :func:`~full_panelled_body_settings` for the definition
of rotation models of vehicle parts.

Note that this panel model does not contain information on panel location or shape, only its area and surface normal,
and is therefore not suitable for computation of panel shadowing of torque computations.


Parameters
----------
surface_normal : np.array
    Panel outward surface normal vector (in specified frame)
area : float
    Panel surface area
frame_orientation : str, default = ""
    Identifier of the frame to which the panel is fixed (if body-fixed frame, this can be left empty)
Returns
-------
BodyPanelGeometrySettings
    Object defining settings for panel geometry






    )";



    } else if(name == "time_varying_panel_geometry" && variant==0) {
        return R"(
        
Factory function for creating settings for a vehicle exterior panel that has time-variable orientation in a given frame.

As :func:`~frame_fixed_panel_geometry`, but with a time-variable outward surface normal :math:`\hat{\mathbf{n}}^{\mathcal{F}}(t)`


Parameters
----------
surface_normal_function : np.array
    Panel outward surface normal vector (in specified frame)
area : float
    Panel surface area
frame_orientation : str, default = ""
    Identifier of the frame to which the panel is fixed (if body-fixed frame, this can be left empty)
Returns
-------
BodyPanelGeometrySettings
    Object defining settings for panel geometry






    )";



    } else if(name == "body_tracking_panel_geometry" && variant==0) {
        return R"(
        
Factory function for creating settings for a vehicle exterior panel where the surface normal tracks a given body.

Factory function for creating settings for a vehicle exterior panel where the surface normal tracks a given body, for instance
to define the surface normal of a solar array to always point towards the Sun, or an antenna to always point towards the Earth.
When using this option, the panel surface normal :math:`\hat{\mathbf{n}}` is computed in an inertial frame based on the tracked
body, and then (if necessary) rotated to the body-fixed frame. 
Note that this panel model does not contain information on panel location or shape, only its area and surface normal,
and is therefore not suitable for computation of panel shadowing of torque computations.


Parameters
----------
body_to_track : str
    Name of the body towards (or away from) which the panel surface normal is to point
towards_tracked_body : bool
    Boolean defining whether the normal vector points towards (if true) or away from (if false) the tracked body
area : float
    Panel surface area
frame_orientation : str, default = ""
    Identifier of the frame in which the panel is defined (with time-variable orientation, defined by tracked body). Note that this option is typically only relevant for internal  book-keeping, and can be left empty
Returns
-------
BodyPanelGeometrySettings
    Object defining settings for panel geometry






    )";



    } else if(name == "body_panel_settings" && variant==0) {
        return R"(
        
Factory function for creating settings for a full panel

Factory function for creating settings for a full panel (presently only geometry and reflection properties). The panel
can also be endowed with an identifier to specify the type of the panel. This has no direct consequences for the model,
but may be useful in estimation, to for instance estimate the reflection properties of all panels specified with identified "MLI"
as a single parameter


Parameters
----------
panel_geometry : BodyPanelGeometrySettings
    Geometric properties of the panel (size and orientation, at least)
panel_reflection_law : BodyPanelReflectionLawSettings, default = None
    Reflection law settings of the panel (default none)
panel_type_id : str, default = ""
    Optional identifier for panel type
Returns
-------
BodyPanelSettings
    Object defining settings for a panel






    )";



    } else if(name == "full_panelled_body_settings" && variant==0) {
        return R"(
        
Factory function for creating settings for a full panelled vehicle exterior

Factory function for creating settings for a full panelled vehicle exterior, taking a list of panel settings,
and (optionally) a list of rotation model settings for vehicle parts. The identifiers for the rotation models
are used to specify the names of part-fixed frames, which are used by the ``frame_orientation`` inputs to factory
functions creating settings for :class:`~BodyPanelGeometrySettings`. For instance, assigning a rotation model
to frame ``LRO_SolarArray`` (dict key for ``part_rotation_model_settings``) allows panels defined in the frame
with this same frame orientation to be defined. The associated rotation model defines rotations from body-fixed
frame :math:`\mathcal{B}` to part-fixed frame :math:`\mathcal{F}_{j}` (for part :math:`j`). The rotation from part-fixed
(where the surface normal is defined) to inertial frame is then computed from 
:math:`\mathbf{R}^{I/\mathcal{F}_{j}}=\mathbf{R}^{I/\mathcal{B}}\mathbf{R}^{\mathcal{B}/\mathcal{F}_{j}}`, where :math:`\mathbf{R}^{I/\mathcal{B}}` 
defines the body's orientation, and :math:`\mathbf{R}^{\mathcal{B}/\mathcal{F}_{j}}` the part orientation (w.r.t. a body-fixed frame)


Parameters
----------
panel_settings : list[BodyPanelSettings]
    List of settings for body panels.
part_rotation_model_settings : dict[str,RotationModelSettings], default = dict()
    Rotation model settings per vehicle part (default empty, indicating no part-fixed frames are defined)
Returns
-------
FullPanelledBodySettings
    Object defining full panelled vehicle exterior






    )";



    } else if(name == "box_wing_panelled_body_settings" && variant==0) {
        return R"(
        
Factory function for creating a simple box-wing spacecraft exterior shape with reflection law settings.

This function creates a :func:`~full_panelled_body_settings` with ``panel_settings`` generated from simple box-wing 
settings. The assumptions behind the box-wing model are:

* The spacecraft shape is defined by a rectangular box (cuboid) and solar array
* The box has its faces parallel to the xy-, xz- and yz-planes
* The solar array surface normal always points towards the Sun
* Each box face has identical reflection law settings, defined by :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.specular_diffuse_body_panel_reflection` settings.
* The solar array has reflection law settings, defined by :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.specular_diffuse_body_panel_reflection` settings.


Parameters
----------
length : float
    Box length (size in body-fixed x-direction).
width : float
    Box width (size in body-fixed y-direction).
height : float
    Box height (size in body-fixed z-direction).
solar_array_area : float
    Surface area of the solar array.
box_specular_reflectivity : float
    Box secular reflectivity :math:`\rho`.
box_diffuse_reflectivity : float
    Box secular reflectivity :math:`\delta`.
solar_array_specular_reflectivity : float
    Solar array secular reflectivity :math:`\rho`.
solar_array_diffuse_reflectivity : float
    Solar array secular reflectivity :math:`\delta`.
box_instantaneous_reradiation : bool
    Boolean denoting whether absorbed radiation is instantaneously retransmitted from box (yes, if true).
solar_array_instantaneous_reradiation : bool
    Boolean denoting whether absorbed radiation is instantaneously retransmitted from solar array (yes, if true).
Returns
-------
FullPanelledBodySettings
    Object defining full panelled vehicle exterior






    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace propagation_setup {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "create_acceleration_models" && variant==0) {
        return R"(
        
Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.

Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
bodies and central bodies are provided as two separate lists with the same order.


Parameters
----------
body_system : SystemOfBodies
    System of bodies to be used in the propagation.
selected_acceleration_per_body : Dict[str, Dict[str, List[AccelerationSettings]]]
    Key-value container, with key denoting the body undergoing the acceleration, and the value containing an additional key-value container, with the body exerting acceleration, and list of acceleration settings exerted by this body.
bodies_to_propagate : list
    List of bodies to propagate.
central_bodies : list
    List of central bodies, each referred to each propagated body in the same order.
Returns
-------
AccelerationMap
    Set of accelerations acting on the bodies to propagate, provided as dual key-value container, similar to the acceleration settings input, but now with ``AccelerationModel`` lists as inner value





Examples
--------
In this example, the acceleration model includes a spherical harmonic (degree and order 5) gravitational acceleration
and aerodynamic acceleration exerted by the Earth, as well as point-mass gravity exerted by the Sun and the Moon.
The variable ``accelerations_settings_vehicle`` denotes the list of bodies exerting accelerations and the types of
accelerations, while the variable ``acceleration_settings`` associates this list with the body undergoing the
acceleration (``"Vehicle"``).

.. code-block:: python 
   
   # Define bodies that are propagated. 
   bodies_to_propagate = ["Vehicle"] 
   
   # Define central bodies. 
   central_bodies = ["Earth"] 
   
   # Define accelerations acting on Vehicle 
   accelerations_settings_vehicle = dict( 
       Sun=[propagation_setup.acceleration.point_mass_gravity()], 
       Moon=[propagation_setup.acceleration.point_mass_gravity()], 
       Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5), 
              propagation_setup.acceleration.aerodynamic()] 
   ) 
   
   # Create global accelerations settings dictionary 
   acceleration_settings = {"Vehicle": accelerations_settings_vehicle} 
   
   # Create acceleration models 
   acceleration_models = propagation_setup.create_acceleration_models( 
       bodies, acceleration_settings,  bodies_to_propagate, central_bodies) 


    )";



    } else if(name == "create_torque_models" && variant==0) {
        return R"(
        
Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.

Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
bodies is provided as a list.


Parameters
----------
body_system : SystemOfBodies
    System of bodies to be used in the propagation.
selected_torque_per_body : Dict[str, Dict[str, List[TorqueSettings]]]
    Key-value container, with key denoting the body undergoing the torque, and the value containing an additional key-value container, with the body exerting torque, and list of torque settings exerted by this body.
bodies_to_propagate : list
    List of bodies to propagate.
Returns
-------
TorqueModelMap
    Set of torques acting on the bodies to propagate, provided as dual key-value container, similar to the torque settings input, but now with ``TorqueModel`` lists as inner value





Examples
--------

In this example, the following torques are exerted on the vehicle: spherical-harmonic gravitational torque
(up to order 4 and degree 4) and aerodynamic torque exerted by the Earth, second-degree gravitational torque
exerted by the Sun and the Moon.
The variable ``torques_settings_vehicle`` denotes the list of bodies exerting torques and the types of
torques, while the variable ``torque_settings`` associates this list with the body undergoing the
torque.

.. code-block:: python 
  
  # Define bodies that are propagated 
  bodies_to_propagate = ["Vehicle"] 
  
  # Define torques per each exerting body 
  torque_settings_vehicle = dict( 
      Sun=[propagation_setup.torque.second_degree_gravitational()], 
      Moon=[propagation_setup.torque.second_degree_gravitational()], 
      Earth=[propagation_setup.torque.spherical_harmonic_gravitational(4, 4), 
             propagation_setup.torque.aerodynamic()] 
  ) 
  
  # Create global torque settings dictionary 
  torque_settings = {"Vehicle": torque_settings_vehicle} 
  
  # Create torque models 
  torque_models = propagation_setup.create_torque_models( 
      bodies, torque_settings,  bodies_to_propagate ) 


    )";



    } else if(name == "create_mass_rate_models" && variant==0) {
        return R"(
        
Function to create a set of mass-rate models from associated settings.

Function to create a set of mass-rate models from a map of bodies and mass-rate model types.
If the mass-rate depends on any acceleration models (e.g. thrust), the acceleration
models must be provided as an input.


Parameters
----------
body_system : SystemOfBodies
    System of bodies to be used in the propagation.
selected_mass_rates_per_body : Dict[str, List[MassRateModelSettings]]
    Key-value container, with key denoting the body with changing mass, and the value containing a list of mass rate settings (in most cases, this list will have only a single entry)
acceleration_models : AccelerationMap
    Sorted list of acceleration models, as created by :func:`create_acceleration_models`
Returns
-------
MassRateModelMap
    Set of mass-rate models, as key-value container, same as the settings input, with the difference that the rate settings objects have been processed into the associated objects calculating the actual mass-rate changes.






    )";



    } else {
        return "No documentation found.";
    }

}


    
namespace acceleration {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "AvailableAcceleration") {
         return R"(

        Enumeration of available acceleration types.

        Enumeration of acceleration types supported by tudat.





     )";


    } else if(name == "AvailableAcceleration.undefined_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.point_mass_gravity_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.aerodynamic_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.cannon_ball_radiation_pressure_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.spherical_harmonic_gravity_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.mutual_spherical_harmonic_gravity_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.polyhedron_gravity_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.ring_gravity_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.thrust_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.relativistic_correction_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.empirical_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.direct_tidal_dissipation_in_central_body_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.direct_tidal_dissipation_in_orbiting_body_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.panelled_radiation_pressure_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.quasi_impulsive_shots_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.solar_sail_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.custom_acceleration_type") {
         return R"(
     )";


    } else if(name == "AvailableAcceleration.radiation_pressure_type") {
         return R"(
     )";




    } else if(name == "AccelerationSettings") {
         return R"(

        Functional base class to define settings for accelerations.

        Class for providing settings for acceleration model. This class is a functional (base) class for
        settings of acceleration models that  require no information in addition to their type.
        Classes defining settings for acceleration models requiring additional information must be derived from this class.
        Bodies exerting and undergoing acceleration are set externally from this class.
        This class can be used for the easy setup of acceleration models
        (see createAccelerationModels.h), but users may also chose to do so manually.
        (Derived) Class members are all public, for ease of access and modification.





     )";





    } else if(name == "SphericalHarmonicAccelerationSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for the spherical harmonic acceleration.

        Class for providing settings for spherical harmonics acceleration model,
        including the maximum degree and order up to which the field is to be expanded. Note that
        the minimum degree and order are currently always set to zero.





     )";





    } else if(name == "MutualSphericalHarmonicAccelerationSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for the mutual spherical harmonic acceleration.

        Class for providing settings for the mutual spherical harmonics acceleration model,
        including the maximum degree and order up to which the fields of the bodies are to be expanded. Note that
        the minimum degree and order are currently always set to zero.





     )";





    } else if(name == "RelativisticAccelerationCorrectionSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for the relativistic acceleration correction.

        Class to provide settings for typical relativistic corrections to the dynamics of an orbiter: the
        Schwarzschild, Lense-Thirring and de Sitter terms (see 'General relativity and Space Geodesy' by L. Combrinck,
        2012).





     )";





    } else if(name == "EmpiricalAccelerationSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for the empirical acceleration.

        Class to provide settings for empirical accelerations. These are expressed in the
        RSW frame, for which the magnitude is determined empirically (typically during an orbit determination process).
        The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
        a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
        the RSW frame.





     )";





    } else if(name == "CustomAccelerationSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for custom acceleration.

        Class to provide settings for custom accelerations. This is done by means of a function and, if necessary,
        an associated scaling function.





     )";





    } else if(name == "DirectTidalDissipationAccelerationSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for direct tidal dissipation acceleration.

        Class to provide settings for direct tidal dissipation accelerations. Creates settings for tidal accelerations.
        The direct of tidal effects in a satellite system is applied directly as an acceleration
        (as opposed to a modification of spherical harmonic coefficients).





     )";





    } else if(name == "ThrustAccelerationSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for thrust acceleration, listing the engine models that are to be used

        Class to provide settings for thrust acceleration, listing the engine models that are to be used





     )";





    } else if(name == "MomentumWheelDesaturationAccelerationSettings") {
         return R"(

        `AccelerationSettings`-derived class to define settings for momentum wheel desaturation acceleration.

        Class to provide settings for momentum wheel desaturation acceleration. Settings for the direction and magnitude
        of the thrust are included.





     )";






    } else if(name == "point_mass_gravity" && variant==0) {
        return R"(
        
Creates settings for the point-mass gravity acceleration.

Creates settings for the point-mass gravity acceleration. The direct acceleration (acceleration w.r.t. an inertial frame) is computed from:

.. math::
   \mathbf{a}=\frac{\mu}{{r}^{2}}\hat{\mathbf{r}}

with :math:`\mathbf{r}` the position vector measured from the center of mass of the body exerting the acceleration. 

The body exerting the acceleration needs to have a gravity field model (:ref:`\`\`gravity_field\`\`` module) defined to use this acceleration. 

Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation \math:`C`, choosing this option may create a direct point-mass attraction (:math:`\mu=\mu_{B}`), a central point-mass attraction (:math:`\mu=\mu_{B}+\mu_{A}`) or a third-body point-mass attraction (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/third_body_acceleration.html>`_ for more details).

Returns
-------
AccelerationSettings
    Acceleration settings object.





Examples
--------
In this example, we define the point mass gravity acceleration exerted by the Earth on the vehicle:

.. code-block:: python 
   
   # Create acceleration dict 
   accelerations_acting_on_vehicle = dict() 
   # Add aerodynamic acceleration exerted by Earth 
   accelerations_acting_on_vehicle["Earth"] = [propagation_setup.acceleration.point_mass_gravity()] 


    )";



    } else if(name == "aerodynamic" && variant==0) {
        return R"(
        
Creates settings for the aerodynamic acceleration.

Creates settings for the aerodynamic acceleration. The acceleration is computed from:

.. math::
   \mathbf{a}=-\frac{1}{m}\mathbf{R}^{(I/\text{Aero})}\left(\frac{1}{2}\rho v_{\text{air}}^{2}S_{ref}\begin{pmatrix} C_{D} \\ C_{S} \\ C_{L}\end{pmatrix}\right)

with :math:`\mathbf{R}^{(I/\text{Aero})}` the rotation matrix from the aerodynamic frame of the body undergoing acceleration to the inertial frame (computed from the body's current state, and the rotation of the body exerting the acceleration), :math:`\rho` the local freesream atmospheric density, :math:`v_{\text{air}}` the airspeed,  :math:`C_{D,S,L}` the drag, side and lift coefficients (which may depend on any number of properties of the body/environment) with reference area :math:`S_{ref}`, and :math:`m` the mass of the body undergoing acceleration
The body exerting the acceleration needs to have an
atmosphere (:ref:`\`\`atmosphere\`\`` module), shape (:ref:`\`\`shape\`\`` module) and rotation model (:ref:`\`\`rotation_model\`\`` module) defined. The body undergoing the acceleration needs to have aerodynamic coefficients (:ref:`\`\`aerodynamic_coefficients\`\`` module) defined.

Returns
-------
AccelerationSettings
    Acceleration settings object.





Examples
--------
In this example, we define the aerodynamic acceleration exerted by the Earth on the vehicle:

.. code-block:: python 
   
   # Create acceleration dict 
   accelerations_acting_on_vehicle = dict() 
   # Add aerodynamic acceleration exerted by Earth 
   accelerations_acting_on_vehicle["Earth"] = [propagation_setup.acceleration.aerodynamic()] 


    )";



    } else if(name == "radiation_pressure" && variant==0) {
        return R"(
        
Creates settings for the radiation pressure acceleration.

Creates settings for the radiation pressure acceleration. This function takes the source model of the body exerting the acceleration, and the target model of the 
body undergoing the acceleration, and links these models to set up the specific acceleration model, regardless of how the source and target models are defined. 
For more extensive details on how this is done, check out the 
`radiation pressure acceleration page <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/radiation_pressure_acceleration.html>`_.

Returns
-------
AccelerationSettings
    Acceleration settings object.






    )";



    } else if(name == "cannonball_radiation_pressure" && variant==0) {
        return R"(
        
Creates settings for the cannonball radiation pressure acceleration.

Creates settings for the radiation pressure acceleration, for which a cannonball model is used. The acceleration is computed from:

.. math::

   \mathbf{a}=\left(\frac{P}{4\pi c}\right)\left(\frac{C_{r}S_{ref}}{m}\right)\frac{\hat{\mathbf{r}}}{r^{2}} 

with :math:`P` the total emitted radiation power for the body exerting the acceleration, :math:`C_{r}` the radiation pressure coefficient with reference area :math:`S_{ref}`, :math:`\mathbf{r}` the vector from the body exerting the acceleration to the body undergoing the acceleration, and :math:`m` the mass of the body undergoing acceleration

In this model,
the effective acceleration is colinear with the vector connecting the source of radiation and the target.
The body undergoing the acceleration needs to have a radiation pressure model defined, while the body emitting
radiation needs to have radiative properties defined.

Returns
-------
AccelerationSettings
    Acceleration settings object.





Examples
--------
In this example, we define the aerodynamic acceleration exerted by the Sun on the vehicle:

.. code-block:: python 
   
   # Create acceleration dict 
   accelerations_acting_on_vehicle = dict() 
   # Add cannonball radiation pressure acceleration exerted by Sun 
   accelerations_acting_on_vehicle["Sun"] = [propagation_setup.acceleration.cannonball_radiation_pressure()] 


    )";



    } else if(name == "spherical_harmonic_gravity" && variant==0) {
        return R"(
        
Creates settings for the spherical harmonic gravity acceleration.

Creates settings for the spherical harmonic gravity acceleration, accounting for a finite (given as input) number
of degrees and orders. The direct acceleration (acceleration w.r.t. an inertial origin) is computed from:

.. math::
   \mathbf{a}=\mathbf{R}^{(I/B)}\nabla^{(B)}U(\mathbf{r})

with :math:`\mathbf{r}` the position vector measured from the center of mass of the body exerting the acceleration, :math:`\mathbf{R}^{(I/B)}` the rotation matrix from body-fixed to inertial frame, and :math:`\nabla^{(B)}` the gradient operator in a body-fixed frame, and :math:`U` the spherical harmonic gravitational potential, expanded up to the provided ``maximum_degree`` and ``maximum_order``. 

The body exerting the acceleration needs to have a spherical harmonic gravity field model (see :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`) and a rotation model (:ref:`\`\`rotation_model\`\`` module) defined. 

Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation :math:`C`, choosing this option may create a direct spherical harmonic attraction (:math:`\mu=\mu_{B}`), a central spherical harmonic attraction (:math:`\mu=\mu_{B}+\mu_{A}`) or a third-body spherical harmonic attraction (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/third_body_acceleration.html>`_ for more details).


Parameters
----------
maximum_degree : int
    Maximum degree of the spherical harmonic expansion.
maximum_order : int
    Maximum order of the spherical harmonic expansion.
Returns
-------
SphericalHarmonicAccelerationSettings
    Spherical harmonic acceleration settings object.





Examples
--------
In this example, we define the spherical harmonic gravity acceleration (where the gravity field is expanded
up to degree 12 and order 6) exerted by the Earth on the vehicle:

.. code-block:: python 
   
   # Define the maximum degree and order 
   maximum_degree = 12 
   maximum_order = 6 
   # Create acceleration dict 
   accelerations_acting_on_vehicle = dict() 
   # Add aerodynamic acceleration exerted by Earth 
   accelerations_acting_on_vehicle["Earth"] = [propagation_setup.acceleration.spherical_harmonic_gravity( 
        maximum_degree,  
        maximum_order)] 


    )";



    } else if(name == "mutual_spherical_harmonic_gravity" && variant==0) {
        return R"(
        
Creates settings for the mutual spherical harmonic gravity acceleration.

Creates settings for the mutual spherical harmonic gravity acceleration. This model computes the total spherical harmonic acceleration exerted by a body :math:`B` on a body :math:`A`, where the influence of the gravity field coefficients of body :math:`A` itself has been included. The model includes couplings between the mass of each body, and the gravity field coefficients of the other body. It does not include the 'figure-figure' interactions (coupling between the two-bodies' gravity field coefficients). It corresponds to the model presented by Lainey et al. (2004); Dirkx et al. (2016).  
The model combines the spherical harmonic accelerations of the two bodies (see :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic`) on each other. The direct acceleration (acceleration w.r.t. an inertial origin) is computed from:

.. math::

   \mathbf{a}={-\frac{\mu_{_{B}}}{{r}^{2}}\hat{\mathbf{r}}}+{\mathbf{R}^{(I/B)}\nabla^{(B)}U_{\hat{B}}(\mathbf{r})}-{\frac{\mu_{_{B}}}{\mu_{_{A}}}\mathbf{R}^{(I/A)}\nabla^{(A)}U_{\hat{A}}(-\mathbf{r})}

where :math:`U_{\hat{B}}` and :math:`U_{\hat{A}}` denote the spherical harmonic gravity fields a degree :math:`>=1` of bodies :math:`B` and :math:`A`, respectively.
Both the body exerting the acceleration and the body undergoing it need to
have spherical harmonic gravity field and rotation models defined.

Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation :math:`C`, choosing this option may create a direct spherical harmonic attraction (as above), a central spherical harmonic attraction (:math:`\mu_{B}\rightarrow\mu_{B}+\mu_{A}`, in the above equation and in :math:`U_{\hat{B}}`) or a third-body spherical harmonic attraction (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/third_body_acceleration.html>`_ for more details).

For the case where a third-body mutual spherical harmonic acceleration,
additional parameters have to be provided that denote the expansion degree/order of the central body (``maximum_degree_central_body`` and ``maximum_order_central_body``)


Parameters
----------
maximum_degree_body_exerting : int
    Maximum degree of the spherical harmonic expansion for the body exerting the acceleration.
maximum_order_body_exerting : int
    Maximum order of the spherical harmonic expansion for the body exerting the acceleration.
maximum_degree_body_undergoing : int
    Maximum degree of the spherical harmonic expansion for the body undergoing the acceleration.
maximum_order_body_undergoing : int
    Maximum order of the spherical harmonic expansion for the body undergoing the acceleration.
maximum_degree_central_body : int, default=0
    Maximum degree of the spherical harmonic expansion for the central body, if needed.
maximum_order_central_body : int, default=0
    Maximum order of the spherical harmonic expansion for the central body, if needed.
Returns
-------
MutualSphericalHarmonicAccelerationSettings
    Spherical harmonic acceleration settings object.





Examples
--------
In this example, we define the spherical harmonic gravity accelerations exerted on Io by Jupiter and vice-versa.

.. code-block:: python 
   
   # Define the maximum degree and order for both bodies 
   maximum_degree_of_io = 12 
   maximum_order_of_io = 12 
   maximum_degree_of_jupiter = 4 
   maximum_order_of_jupiter = 4 
   # Create acceleration dict 
   acceleration_settings_on_io = dict() 
   # Add aerodynamic acceleration exerted by Earth 
   acceleration_settings_on_io["Jupiter"] = [propagation_setup.acceleration.mutual_spherical_harmonic_gravity( 
        maximum_degree_of_jupiter, 
        maximum_order_of_jupiter, 
        maximum_degree_of_io, 
        maximum_order_of_io)] 

For the case where the mutual spherical harmonic acceleration is a third-body acceleration,
additional parameters have to be provided to denote the expansion of the spherical harmonics of the central body.
In the following example, we consider the spherical harmonic gravity acceleration mutually exerted between
Ganymede and Io when propagating w.r.t. Jupiter:

.. code-block:: python 
   
   # Define the maximum degree and order for both bodies 
   maximum_degree_of_jupiter = 4 
   maximum_order_of_jupiter = 4 
   maximum_degree_of_ganymede = 4 
   maximum_order_of_ganymede = 4 
   maximum_degree_of_io = 12 
   maximum_order_of_io = 12 
   # Create acceleration dict 
   acceleration_settings_on_io = dict() 
   # Add the acceleration to the dict 
   acceleration_settings_on_io["Jupiter"] = [propagation_setup.acceleration.mutual_spherical_harmonic_gravity( 
        maximum_degree_of_jupiter, 
        maximum_order_of_jupiter, 
        maximum_degree_of_ganymede, 
        maximum_order_of_ganymede, 
        maximum_degree_of_io, 
        maximum_order_of_io)] 


    )";



    } else if(name == "polyhedron_gravity" && variant==0) {
        return R"(
        
Creates settings for the polyhedron gravity acceleration.

Creates settings for the polyhedron gravity acceleration, which follows from defining a body to have polyhedral gravity.

Returns
-------
AccelerationSettings
    Acceleration settings object.






    )";



    } else if(name == "ring_gravity" && variant==0) {
        return R"(
        
Creates settings for the ring gravity acceleration.

Creates settings for the ring gravity acceleration, which follows from defining a body to have ring gravity.

Returns
-------
AccelerationSettings
    Acceleration settings object.






    )";



    } else if(name == "relativistic_correction" && variant==0) {
        return R"(
        
Creates settings for the relativistic acceleration correction.

Creates settings for typical relativistic acceleration corrections: the Schwarzschild, Lense-Thirring and de
Sitter terms, where each of the three terms can be toggled on or of (see 'General relativity and Space Geodesy' by L. Combrinck, 2012). It implements the model of
2010 Conventions (chapter 10, section 3). Here, the ‘primary body’ for a planetary orbiter should always be set
as the Sun (only relevant for de Sitter correction). The angular momentum vector of the orbited body is only
relevant for Lense-Thirring correction.


Parameters
----------
use_schwarzschild : bool, default=False
    Boolean defining whether or not to use the Schwarzschild contribution to the acceleration correction
use_lense_thirring : bool
    Boolean defining whether or not to use the Lense-Thirring contribution to the acceleration correction
use_de_sitter : bool, default=False
    Boolean defining whether or not to use the de Sitter contribution to the acceleration correction
de_sitter_central_body : str, default=""
    Body used as 'perturbed' in the calculation of the de Sitter acceleration. For the case of an Earth-orbiting satellite, this would be the Sun
lense_thirring_angular_momentum : numpy.ndarray, default=numpy.array([0, 0, 0])
    Angular momentum vector per unit mass (in global frame) that is to be used for the calculation of the Lense-Thirring acceleration
Returns
-------
RelativisticAccelerationCorrectionSettings
    Relativistic acceleration correction settings object.





Examples
--------
In this example, we define the relativistic correction acceleration for a Mars orbiter:

.. code-block:: python 
   
   # Select terms to be used 
   use_schwarzschild = True 
   use_lense_thirring = True 
   use_de_sitter = True 
   # Define central body for De-Sitter term 
   de_sitter_central_body = "Sun" 
   # Define angular momentum for the Lense-Thirring term 
   lense_thirring_angular_momentum = ... # numpy.ndarray 3D vector 
   # Create acceleration dict 
   acceleration_settings_on_vehicle = dict() 
   # Add the acceleration to the dict 
   acceleration_settings_on_vehicle["Mars"] = [propagation_setup.acceleration.relativistic_correction( 
      use_schwarzschild, 
      use_lense_thirring,  
      use_de_sitter,  
      de_sitter_central_body,  
      lense_thirring_angular_momentum)] 


    )";



    } else if(name == "empirical" && variant==0) {
        return R"(
        
Creates settings for empirical acceleration.

Creates settings for empirical accelerations. These are expressed in the
RSW frame, for which the magnitude is determined empirically (typically during an orbit determination process).
The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
the RSW frame. The empirical acceleration is calculated as:

 .. math::

    \mathbf{a}=R^{I/RSW}\left(\mathbf{a}_{\text{const.}}+\mathbf{a}_{\sin}\sin\theta+\mathbf{a}_{\cos}\cos\theta \right)

Here, :math:`R^{I/RSW}` is the rotation matrix from the RSW frame (of the body undergoing the acceleration w.r.t. the
body exerting the acceleration), :math:`theta` is the true anomaly, and the three constituent acceleration vectors are
the inputs provided in the above code block. The body 'exerting' the acceleration is considered to be the
central body, w.r.t. which the true anomaly is calculated.


Parameters
----------
constant_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
    Constant term, defined in the RSW frame.
sine_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
    Sine term (function of the true anomaly), defined in the RSW frame..
cosine_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
    Cosine term (function of the true anomaly), defined in the RSW frame..
Returns
-------
EmpiricalAccelerationSettings
    Empirical acceleration settings object.





Examples
--------
In this example, we define the empirical acceleration acting on the vehicle. The body that 'exerts' the acceleration
is here used to determine the body w.r.t. which the true anomaly has o be calculated when determining the sine/cosine
contributions. This central body must be endowed with a gravity field (so that it possesses a gravitational parameter
for the Cartesian to Keplerian conversion)

.. code-block:: python 
   
   # Define the nine terms (constant, sine and cosine) 
   constant_acceleration = ... # 3D numpy.ndarray vector 
   sine_acceleration = ... # 3D numpy.ndarray vector 
   cosine_acceleration = ... # 3D numpy.ndarray vector 
   # Create acceleration dict 
   acceleration_settings_on_vehicle = dict() 
   # Add the acceleration to the dict 
   acceleration_settings_on_vehicle["Sun"] = [propagation_setup.acceleration.empirical( 
       constant_acceleration,  
       sine_acceleration,  
       cosine_acceleration)] 


    )";



    } else if(name == "custom_acceleration" && variant==0) {
        return R"(
        
Creates settings for custom acceleration.

Creates settings for a custom accelerations, this acceleration must be parameterized as a function of time,
and expressed with an inertial orientation.


Parameters
----------
acceleration_function : callable[[float], list]
    Custom acceleration function with time as an independent variable, returning the acceleration in an inertial frame (*e.g.* with global frame orientation) as a function of time.
Returns
-------
CustomAccelerationSettings
    Custom acceleration settings object.





Examples
--------
In this example, we define a simple, standalone, custom acceleration (depending only on time),
with the following (arbitrary, for the sake of example) mathematical definition:

 .. math::

    \mathbf{a}=\begin{pmatrix}C\\0\\0 \end{pmatrix}\sin\left(\frac{t-t_{0}}{T}\right)

with :math:`C=10^{-8}`, :math:`t_{0}=0` and :math:`T=86400`.
More complex custom acceleration functions can be created by, for instance, extracting the
custom function from a user-defined class, which may in turn have other simulation objects
(*e.g.* :class:`SystemOfBodies`) as members, allowing the custom acceleration to depend on
the current simulation state/environment

.. code-block:: python 
   
   # Define custom function 
   def custom_function( current_time ): 
       period = 86400.0 
       reference_time = 0.0 
       phase = np.pi * ( current_time - reference_time ) / period 
       return np.array([1.0E-8, 0.0, 0.0]) * np.sin(phase) 
   
   acceleration_settings_on_vehicle["Vehicle"] = [propagation_setup.acceleration.custom_acceleration( 
       custom_function)] 


    )";



    } else if(name == "direct_tidal_dissipation_acceleration" && variant==0) {
        return R"(
        
Creates settings for tidal acceleration.

Creates settings for tidal accelerations. The direct of tidal effects in a satellite system is applied directly as
an acceleration (as opposed to a modification of spherical harmonic coefficients).
The model is based on Lainey et al. (2007, 2012). It can compute the acceleration due to tides, and in
particular tidal dissipation, on a planetary satellite. The acceleration computed can account for either the
effect of tide raised on the satellite by the planet or on the planet by the satellite. The satellite is assumed
to be tidally locked to the planet.


Parameters
----------
k2_love_number : float
    Value of the k2 Love number.
time_lag : float
    Value of the tidal time lag.
include_direct_radial_component : bool, default=True
    It denotes whether the term independent of the time lag is to be computed.
use_tide_raised_on_planet : bool, default=True
    It denotes whether the tide raised on the planet is to be modelled (if true) or the tide raised on the satellite (if false).
Returns
-------
DirectTidalDissipationAccelerationSettings
    Direct tidal dissipation acceleration settings object.





Examples
--------
In this example, we define the tidal dissipation exerted by Jupiter on Io directly, instead of computing it
through the spherical harmonic gravity:

.. code-block:: python 
   
   # Define parameters 
   love_number = 0.1 
   time_lag = 100.0 
   # Add entry to acceleration settings dict 
   acceleration_settings_on_io["Jupiter"] = [propagation_setup.acceleration.direct_tidal_dissipation( 
      love_number,  
      time_lag,  
      False,  
      False)] 


    )";



    } else if(name == "quasi_impulsive_shots_acceleration" && variant==0) {
        return R"(
        
Creates settings for incorporating quasi-impulsive shots into the acceleration.

The acceleration model is purpose-built to represent short bursts of thrust, such as a momentum wheel desaturation.
A typical use case is precise orbit determination, but the functionality can be used just as well in propagation
(for instance to model an impulsive manuever in a continuous manner when going from preliminary modelling to
full modelling). The thrust is modelled similarly to Fig. 3 of Alessi et al. (2012), with the main difference
being that a third-order polynomial to go from zero acceleration to the maximum acceleration level is employed.
By using a 3rd-order polynomial and imposing continuity in the value and first derivative of the acceleration,
defining the rise time (time it takes acceleration to go from 0 to its maximum level), the total time where
there is non-zero thrust (total maneuver time), and the total Delta V exerted by a single maneuver,
the acceleration profile is fully defined.


Parameters
----------
thrust_mid_times : list[float]
    Set of middle point in times in the maneuver denoting the epoch of each maneuver.
delta_v_values : list[numpy.ndarray]
    Set of delta V, one for each maneuver.
total_maneuver_time : float
    Total duration of every maneuver.
maneuver_rise_time : float
    Time taken by the acceleration to go from zero to its maximum level.
Returns
-------
MomentumWheelDesaturationAccelerationSettings
    Momentum wheel desaturation acceleration settings object.





Examples
--------
In this example, we define an acceleration model to represent two quasi-impulsive shots, with a total duration of
30 seconds, and a rise time of 5 seconds. The maneuvers are to be done at :math:`t=86400` and :math:`t=2*86400`.
The first maneuver is exclusively in :math:`x`-direction, and the second one exclusively in :math:`y`-direction,
with both maneuvers having a magnitude of 1 mm/s

.. code-block:: python 
   
   # Define the quasi-impulsive shot settings 
   mid_times = [ 86400.0, 2.0 * 86400.0] 
   delta_v_values = [ np.array([1.0E-3, 0.0, 0.0]), np.array([0.0, 1.0E-3, 0.0]) ] 
   maneuver_duration = 30 
   maneuver_duration = 5 
   
   # Create acceleration dict 
   acceleration_settings_on_spacecraft = dict() 
   
   # Add quasi-impulsive acceleration exerted by Spacecraft itself (!) 
   acceleration_settings_on_spacecraft["Spacecraft"] = [propagation_setup.acceleration.quasi_impulsive_shots_acceleration( 
        mid_times, 
        delta_v_values, 
        maneuver_duration, 
        maneuver_duration)] 


    )";



    } else if(name == "thrust_from_engines" && variant==0) {
        return R"(
        
Creates settings for thrust acceleration using a list of engine models.

Creates settings for thrust acceleration using a list of engine models. See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_
for more details on the definition of a thrust model in Tudat.


Parameters
----------
engine_names : List[str]
    List of engine names to use when computing thrust.
Returns
-------
ThrustAccelerationSettings
    Thrust acceleration settings object.






    )";



    } else if(name == "thrust_from_engine" && variant==0) {
        return R"(
        
Creates settings for thrust acceleration using a single engine models.

Creates settings for thrust acceleration using a single engine models. See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_
for more details on the definition of a thrust model in Tudat.


Parameters
----------
engine_name : str
    Name of engine to use when computing thrust.
Returns
-------
ThrustAccelerationSettings
    Thrust acceleration settings object.






    )";



    } else if(name == "thrust_from_all_engines" && variant==0) {
        return R"(
        
Creates settings for thrust acceleration using a single engine models.

Creates settings for thrust acceleration by combining thrust from all engines defined in the body. See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_
for more details on the definition of a thrust model in Tudat.

Returns
-------
ThrustAccelerationSettings
    Thrust acceleration settings object.






    )";



    } else if(name == "yarkovsky" && variant==0) {
        return R"(
        
Creates settings for the Yarkovsky acceleration.

Creates settings for the Yarkovsky acceleration, which is calculated based on Pérez-Hernández & Benet (2022). The acceleration is only
considered in the tangential direction and is proportional to:

.. math::

   \mathbf{a}=A_{2} \cdot (\frac{r_{0}}{r_{S}})^{2}

where :math:`A_{2}` is the Yarkovsky parameter, :math:`r_{0} = 1` AU and :math:`r_{S}` is the heliocentric distance in AU.


Parameters
----------
yarkovsky_parameter : float
    Value of the Yarkovsky parameter.
Returns
-------
AccelerationSettings
    Acceleration settings object.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace dependent_variable {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "PropagationDependentVariables") {
         return R"(

        Enumeration of available propagation dependent variables.

        Enumeration of propagation dependent variables supported by tudat.





     )";


    } else if(name == "PropagationDependentVariables.mach_number_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.altitude_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.airspeed_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.local_density_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.relative_speed_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.relative_position_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.relative_distance_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.relative_velocity_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.radiation_pressure_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.total_acceleration_norm_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.single_acceleration_norm_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.total_acceleration_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.single_acceleration_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.aerodynamic_force_coefficients_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.aerodynamic_moment_coefficients_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.rotation_matrix_to_body_fixed_frame_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.intermediate_aerodynamic_rotation_matrix_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.relative_body_aerodynamic_orientation_angle_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.body_fixed_airspeed_based_velocity_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.total_aerodynamic_g_load_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.local_temperature_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.geodetic_latitude_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.control_surface_deflection_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.total_mass_rate_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.tnw_to_inertial_frame_rotation_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.periapsis_altitude_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.apoapsis_altitude_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.total_torque_norm_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.single_torque_norm_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.total_torque_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.single_torque_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.body_fixed_groundspeed_based_velocity_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.keplerian_state_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.modified_equinoctial_state_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.spherical_harmonic_acceleration_terms_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.spherical_harmonic_acceleration_norm_terms_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.body_fixed_relative_cartesian_position_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.body_fixed_relative_spherical_position_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.total_gravity_field_variation_acceleration_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.single_gravity_field_variation_acceleration_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.single_gravity_field_variation_acceleration_terms_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.acceleration_partial_wrt_body_translational_state_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.local_dynamic_pressure_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.euler_angles_to_body_fixed_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.current_body_mass_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.radiation_pressure_coefficient_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.gravity_field_potential_type") {
         return R"(
     )";


    } else if(name == "PropagationDependentVariables.gravity_field_laplacian_of_potential_type") {
         return R"(
     )";




    } else if(name == "VariableSettings") {
         return R"(

        Functional base class to define settings for variables.

        This class is a functional base class for defining settings for variables.
        Any variable that requires additional information in addition to what can be provided here, should be defined by a
        dedicated derived class.





     )";





    } else if(name == "SingleDependentVariableSaveSettings") {
         return R"(

        `VariableSettings`-derived class to define settings for dependent variables that are to be saved during propagation.

        Functional base class for defining settings for dependent variables that are to be computed and saved during propagation.
        Any dependent variable that requires additional information in addition to what can be provided here, should be
        defined by a dedicated derived class.





     )";





    } else if(name == "SingleAccelerationDependentVariableSaveSettings") {
         return R"(

        `SingleDependentVariableSaveSettings`-derived class to save a single acceleration (norm or vector) during propagation.

        Class to define settings for saving a single acceleration (norm or vector) during propagation. Note: this acceleration is returned in the inertial frame!




     )";






    } else if(name == "mach_number" && variant==0) {
        return R"(
        
Function to add the Mach number to the dependent variables to save.

Function to add the Mach number to the dependent variables to save. The calculation of the altitude uses the atmosphere model of the central body and the current state of the body for which the Mach number is to be calculated.

Parameters
----------
body : str
    Body whose Mach number is to be saved.
central_body : str
    Body with atmosphere with respect to which the Mach number is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.





Examples
--------

To create settings for saving of a Mach number of a body name 'Spacecraft'
w.r.t. the atmosphere of body 'Earth', use:

.. code-block:: python
   :emphasize-lines: 18

   # Define save settings for Mach number
   propagation_setup.dependent_variable.mach_number( "Spacecraft", "Earth" )


    )";



    } else if(name == "altitude" && variant==0) {
        return R"(
        
Function to add the altitude to the dependent variables to save.

Function to add the altitude to the dependent variables to save. The calculation of the altitude uses the shape model of the central body and the current state of the body for which the altitude is to be calculated.

Parameters
----------
body : str
    Body whose altitude is to be saved.
central_body : str
    Body with respect to which the altitude is computed (requires this body to have a shape model defined).
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "airspeed" && variant==0) {
        return R"(
        
Function to add the airspeed to the dependent variables to save.

Function to add the airspeed to the dependent variables to save. The calculation of the airspeed uses the rotation and wind models of the central body (to determine the motion of the atmosphere in inertial space), and the current state of the body for which the airspeed is to be calculated.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with atmosphere with respect to which the airspeed is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "body_fixed_airspeed_velocity" && variant==0) {
        return R"(
        
Function to add the airspeed velocity vector to the dependent variables to save.

Function to add the airspeed velocity vector to the dependent variables to save. The airspeed velocity vector is *not provided in an inertial frame*, but instead a frame centered on, and fixed to, the central body. It defines the velocity vector of a body w.r.t. the relative atmosphere It requires the central body to have an atmosphere.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the airspeed is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "body_fixed_groundspeed_velocity" && variant==0) {
        return R"(
        
Function to add the groundspeed velocity vector to the dependent variables to save.

Function to add the groundspeed velocity vector to the dependent variables to save. The groundspeed velocity vector is *not provided in an inertial frame*, but instead a frame centered on, and fixed to, the central body. It defines the velocity vector of a body w.r.t. 'the ground' or (alternatively and identically) the relative atmosphere in the case the atmosphere would be perfectly co-rotating with the central body.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the groundspeed is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "density" && variant==0) {
        return R"(
        
Function to add the local freestream density to the dependent variables to save.

Function to add the freestream density (at a body's position) to the dependent variables to save. The calculation of the density uses the atmosphere model of the central body, and the current state of the body for which the density is to be calculated.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
body_with_atmosphere : str
    Body with atmosphere with respect to which the density is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "temperature" && variant==0) {
        return R"(
        
Function to add the local freestream temperature to the dependent variables to save.

Function to add the freestream temperature (at a body's position) to the dependent variables to save. The calculation of the temperature uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
body_with_atmosphere : str
    Body with atmosphere with respect to which the temperature is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "dynamic_pressure" && variant==0) {
        return R"(
        
Function to add the local freestream dynamic pressure to the dependent variables to save.

Function to add the freestream dynamic pressure (at a body's position) to the dependent variables to save. The calculation of the temperature uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
body_with_atmosphere : str
    Body with atmosphere with respect to which the dynamic pressure is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "local_aerodynamic_g_load" && variant==0) {
        return R"(
        
Function to add the total aerodynamic G-load to the dependent variables to save.

Function to add the total aerodynamic G-load of a body to the dependent variables to save. The calculation uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
body_with_atmosphere : str
    Body with atmosphere exerting the aerodynamic acceleration for which the g-load is to be computed
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "relative_position" && variant==0) {
        return R"(
        
Function to add the relative position vector to the dependent variables to save.

Function to add a body's relative position vector with respect to a second body to the dependent variables to save. The relative position is computed between the bodies' centers of mass.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
relative_body : str
    Body with respect to which the relative position is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "relative_distance" && variant==0) {
        return R"(
        
Function to add the relative distance to the dependent variables to save.

Function to add a body's relative distance (norm of the position vector) with respect to a second body to the dependent variables to save. The relative distance is computed between the bodies' centers of mass.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
relative_body : str
    Body with respect to which the relative distance is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "relative_velocity" && variant==0) {
        return R"(
        
Function to add the relative velocity vector to the dependent variables to save.

Function to add a body's relative velocity vector with respect to a second body to the dependent variables to save. The relative velocity is computed between the bodies' centers of mass.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
relative_body : str
    Body with respect to which the relative velocity is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "relative_speed" && variant==0) {
        return R"(
        
Function to add the relative speed to the dependent variables to save.

Function to add a body's relative speed (norm of the relative velocity vector) with respect to a second body to the dependent variables to save. The relative speed is computed between the bodies' centers of mass.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
relative_body : str
    Body with respect to which the relative speed is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "keplerian_state" && variant==0) {
        return R"(
        
Function to add the Keplerian state to the dependent variables to save.

Function to add the Keplerian state to the dependent variables to save. The Keplerian state is returned in this order: 1: Semi-major Axis. 2: Eccentricity. 3: Inclination. 4: Argument of Periapsis. 5. Right Ascension of the Ascending Node. 6: True Anomaly.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the Keplerian state is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "modified_equinoctial_state" && variant==0) {
        return R"(
        
Function to add the modified equinoctial state to the dependent variables to save.

Function to add the modified equinoctial state to the dependent variables to save. The value of the parameter I is automatically chosen as +1 or -1, depending on whether the inclination is smaller or larger than 90 degrees. The elements are returned in the order :math:`p`, :math:`f`, :math:`g`, :math:`h`, :math:`k`, :math:`L`

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the modified equinoctial state is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "single_acceleration" && variant==0) {
        return R"(
        
Function to add a single acceleration to the dependent variables to save.

Function to add a single acceleration vector to the dependent variables to save. The requested acceleration is defined by its type, and the bodies undergoing and exerting the acceleration. This acceleration vector represents the acceleration in 3D in the inertial reference frame. NOTE: When requesting a third-body perturbation be saved, you may use either the direct acceleration type, or the third body type. For instance, for saving a point-mass third-body perturbation, you may specify either ``point_mass_gravity_type`` or ``third_body_point_mass_gravity_type`` as acceleration type.

Parameters
----------
acceleration_type : AvailableAcceleration
    Acceleration type to be saved.
body_undergoing_acceleration : str
    Body undergoing acceleration.
body_exerting_acceleration : str
    Body exerting acceleration.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.





Examples
--------

To create settings for saving a point mass acceleration acting on body called 'Spacecraft',
exerted by a body named 'Earth', use:

.. code-block:: python
   :emphasize-lines: 18

   # Define save settings for point-mass acceleration on Spacecraft by Earth
   propagation_setup.dependent_variable.single_acceleration(
           propagation_setup.acceleration.point_mass_gravity_type, 'Spacecraft', 'Earth' )


    )";



    } else if(name == "single_acceleration_norm" && variant==0) {
        return R"(
        
Function to add a single scalar acceleration to the dependent variables to save.

Function to add a single scalar acceleration (norm of the acceleration vector) to the dependent variables to save. The requested acceleration is defined by its type, and the bodies undergoing and exerting the acceleration. NOTE: When requesting a third-body perturbation be saved, you may use either the direct acceleration type, or the third body type. For instance, for saving a point-mass third-body perturbation, you may specify either ``point_mass_gravity_type`` or ``third_body_point_mass_gravity_type`` as acceleration type.

Parameters
----------
acceleration_type : AvailableAcceleration
    Acceleration type to be saved
body_undergoing_acceleration : str
    Body undergoing acceleration.
body_exerting_acceleration : str
    Body exerting acceleration.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.





Examples
--------

To create settings for saving norm of a point mass acceleration acting on body called 'Spacecraft',
exerted by a body named 'Earth', use:

.. code-block:: python
   :emphasize-lines: 18

   # Define save settings for point-mass acceleration on Spacecraft by Earth
   propagation_setup.dependent_variable.single_acceleration_norm(
           propagation_setup.acceleration.point_mass_gravity_type, 'Spacecraft', 'Earth' )


    )";



    } else if(name == "total_acceleration_norm" && variant==0) {
        return R"(
        
Function to add the total scalar acceleration (norm of the vector) acting on a body to the dependent variables to save.


Parameters
----------
body : str
    Body undergoing acceleration.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "total_acceleration" && variant==0) {
        return R"(
        
Function to add the total acceleration vector acting on a body to the dependent variables to save.


Parameters
----------
body : str
    Body undergoing acceleration.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "single_torque_norm" && variant==0) {
        return R"(
        
Function to add a single torque (norm of the torque vector) to the dependent variables to save.


Parameters
----------
torque_type : AvailableTorque
    Torque type to be saved.
body_undergoing_torque : str
    Body undergoing torque.
body_exerting_torque : str
    Body exerting torque.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "single_torque" && variant==0) {
        return R"(
        
Function to add a single torque vector to the dependent variables to save.


Parameters
----------
torque_type : AvailableTorque
    Torque type to be saved.
body_undergoing_torque : str
    Body undergoing torque.
body_exerting_torque : str
    Body exerting torque.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "total_torque_norm" && variant==0) {
        return R"(
        
Function to add the total torque (norm of the torque vector) to the dependent variables to save.


Parameters
----------
body : str
    Body whose dependent variable should be saved.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "total_torque" && variant==0) {
        return R"(
        
Function to add the total torque vector to the dependent variables to save.


Parameters
----------
body : str
    Body whose dependent variable should be saved.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "spherical_harmonic_terms_acceleration" && variant==0) {
        return R"(
        
Function to add single degree/order contributions of a spherical harmonic acceleration vector to the dependent variables to save.

Function to add single degree/order contributions of a spherical harmonic acceleration vector to the dependent variables to save. The spherical harmonic acceleration consists of a (truncated) summation of contributions at degree :math:`l` and order :math:`m`. Using this function, you can save the contributions of separate :math:`l,m` entries to the total acceleration. For instance, when requesting dependent variables for :math:`l,m=2,2`, the contribution due to the combined influence of :math:`ar{C}_{22}` and `ar{S}_{22}` are provided 

Parameters
----------
body_undergoing_acceleration : str
    Body undergoing acceleration.
body_exerting_acceleration : str
    Body exerting acceleration.
component_indices : list[tuple]
    Tuples of (degree, order) indicating the terms to save.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.





Examples
--------

To create settings for saving spherical harmonic acceleration contributions of degree/order
2/0, 2/1 and 2/2, acting on a body names 'Spacecraft', exerted by a body named 'Earth',
use the following for the acceleration. The resulting dependent variable will contain
nine entries (three acceleration components for 2/0, 2/1 and 2/2, respectively).

.. code-block:: python
   :emphasize-lines: 18

   # Define degree/order combinations for which to save acceleration contributions
   spherical_harmonic_terms = [ (2,0), (2,1), (2,2) ]

   # Define save settings for separate spherical harmonic contributions
   propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration( "Spacecraft", "Earth", spherical_harmonic_terms )


    )";



    } else if(name == "spherical_harmonic_terms_acceleration_norm" && variant==0) {
        return R"(
        
Function to add a single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.

Function to add single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.

Parameters
----------
body_undergoing_acceleration : str
    Body undergoing acceleration.
body_exerting_acceleration : str
    Body exerting acceleration.
component_indices : list[tuple]
    Tuples of (degree, order) indicating the terms to save.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.





Examples
--------

To create settings for saving spherical harmonic acceleration contributions of degree/order
2/0, 2/1 and 2/2, acting on a body names 'Spacecraft', exerted by a body named 'Earth',
use the following for the acceleration. The resulting dependent variable will contain
three entries (one acceleration norm for 2/0, 2/1 and 2/2, respectively).

.. code-block:: python
   :emphasize-lines: 18


   # Define degree/order combinations for which to save acceleration contributions
   spherical_harmonic_terms = [ (2,0), (2,1), (2,2) ]

   # Define save settings for separate spherical harmonic contributions
   propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration_norm( "Spacecraft", "Earth", spherical_harmonic_terms )


    )";



    } else if(name == "aerodynamic_force_coefficients" && variant==0) {
        return R"(
        
Function to add the aerodynamic force coefficients to the dependent variables to save.

Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic coefficient interface to be defined for the vehicle. The coefficients are returned in the following order: C_D, C_S, C_l (if coefficient interface defined in aerodynamic frame), or C_X, C_Y, C_Z (if coefficient interface defined in body frame).

Parameters
----------
body : str
    Body undergoing acceleration.
central_body : str
    Body exerting acceleration (e.g. body with atmosphere).
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "aerodynamic_moment_coefficients" && variant==0) {
        return R"(
        
Function to add the aerodynamic moment coefficients to the dependent variables to save.

Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic coefficient interface to be defined for the vehicle. The coefficients are returned in the following order: C_l, C_m, C_n , respectively about the X, Y, Z axes of the body-fixed frame, see (see Mooij, 1994 [1]_)

Parameters
----------
body : str
    Body undergoing acceleration.
central_body : str
    Body exerting acceleration (e.g. body with atmosphere).
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "latitude" && variant==0) {
        return R"(
        
Function to add the latitude to the dependent variables to save.

Function to add the latitude of a body, in the body-fixed frame of a central body, to the dependent variables to save.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the latitude is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "geodetic_latitude" && variant==0) {
        return R"(
        
Function to add the geodetic latitude to the dependent variables to save.

Function to add the geodetic latitude, in the body-fixed frame of a central body, to the dependent variables to save. If the central body has a spherical shape model, this value is identical to the latitude. If the central body has an oblate spheroid shape model, the calculation of the geodetic latitude uses the flattening of the this shape model to determine the geodetic latitude

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the geodetic latitude is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "longitude" && variant==0) {
        return R"(
        
Function to add the longitude to the dependent variables to save.

Function to add the longitude of a body, in the body-fixed frame of a central body, to the dependent variables to save.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the longitude is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "heading_angle" && variant==0) {
        return R"(
        
Function to add the heading angle to the dependent variables to save.

Function to add the heading angle to the dependent variables to save, as defined by Mooij, 1994 [1]_ .

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the heading angle is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "flight_path_angle" && variant==0) {
        return R"(
        
Function to add the flight path angle to the dependent variables to save.

Function to add the flight path angle to the dependent variables to save, as defined by Mooij, 1994 [1]_ .

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the flight path angle is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "angle_of_attack" && variant==0) {
        return R"(
        
Function to add the angle of attack to the dependent variables to save.

Function to add the angle of attack angle to the dependent variables to save, as defined by Mooij, 1994 [1]_ .

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the angle of attack is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "sideslip_angle" && variant==0) {
        return R"(
        
Function to add the sideslip angle to the dependent variables to save, as defined by Mooij, 1994 [1]_ .


Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the sideslip angle is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "bank_angle" && variant==0) {
        return R"(
        
Function to add the bank angle to the dependent variables to save, as defined by Mooij, 1994 [1]_ .


Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the bank angle is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "radiation_pressure" && variant==0) {
        return R"(
        
Function to add the radiation pressure to the dependent variables to save.

Function to add the local radiation pressure, in N/m^2, to the dependent variables to save. It requires a radiation source model to be defined for the radiating body.

Parameters
----------
target_body : str
    Name of body at the location of which the radiation pressure is to be returned
source_body : str
    Name of body from which the radiation originates that causes the radiation pressure at the target body
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "total_gravity_field_variation_acceleration" && variant==0) {
        return R"(
        
Function to add the acceleration induced by the total time-variability of a gravity field to the dependent variables to save.

Function to add the acceleration induced by the total time-variability of a gravity field to the dependent variables to save. This function does not distinguish between different sources of variations of the gravity field, and takes the full time-variation when computing the contribution to the acceleration. To select only one contribution, use the :func:`single_gravity_field_variation_acceleration` function.

Parameters
----------
body_undergoing_acceleration : str
    Body whose dependent variable should be saved.
body_exerting_acceleration : str
    Body exerting the acceleration.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "single_gravity_field_variation_acceleration" && variant==0) {
        return R"(
        
Function to add the acceleration induced by a single time-variability of a gravity field to the dependent variables to save.

Function to add the acceleration induced by a single time-variability of a gravity field to the dependent variables to save. The user specifies the type of variability for which the induced acceleration is to be saved.

Parameters
----------
body_undergoing_acceleration : str
    Body whose dependent variable should be saved.
body_exerting_acceleration : str
    Body exerting the acceleration.
deformation_type : BodyDeformationTypes
    Type of gravity field variation for which the acceleration contribution is to be saved
identifier : str, default=""
    Identifier for the deformation type. To be used in case multiple realizations of a single variation type are present in the given body. Otherwise, this entry can be left empty
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "single_per_term_gravity_field_variation_acceleration" && variant==0) {
        return R"(
        
Function to add the acceleration induced by a single time-variability of a gravity field, at a given list of degrees/orders, to the dependent variables to save. This combines the functionality of the :func:`single_gravity_field_variation_acceleration` and :func:`spherical_harmonic_terms_acceleration` variables


Parameters
----------
body_undergoing_acceleration : str
    Body whose dependent variable should be saved.
body_exerting_acceleration : str
    Body exerting the acceleration.
component_indices : list[tuple]
    Tuples of (degree, order) indicating the terms to save.
deformation_type : BodyDeformationTypes
    Type of gravity field variation for which the acceleration contribution is to be saved
identifier : str, default=""
    Identifier for the deformation type. To be used in case multiple realizations of a single variation type are present in the given body. Otherwise, this entry can be left empty
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "inertial_to_body_fixed_rotation_frame" && variant==0) {
        return R"(
        
Function to add the rotation matrix from inertial to body-fixed frame to the dependent variables to save.

Function to add the rotation matrix from inertial to body-fixed frame to the dependent variables to save. This requires the rotation of the body to be defined (either in the environment or the state vector). NOTE: a rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`),

Parameters
----------
body : str
    Body for which the rotation matrix is to be saved.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "tnw_to_inertial_rotation_matrix" && variant==0) {
        return R"(
        
Function to add the rotation matrix from the TNW to the inertial frame to the dependent variables to save.

Function to add the rotation matrix from the TNW to the inertial frame to the dependent variables to save. It has the x-axis pointing along the velocity vector, the z-axis along the orbital angular momentum vector, and the y-axis completing the right-handed system. NOTE: a rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`),

Parameters
----------
body : str
    Body for which the rotation matrix is to be saved.
central_body : str
    Body with respect to which the TNW frame is determined.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "rsw_to_inertial_rotation_matrix" && variant==0) {
        return R"(
        
Function to add the rotation matrix from the RSW to the inertial frame to the dependent variables to save.

Function to add the rotation matrix from the RSW to the inertial frame to the dependent variables to save. It has the x-axis pointing along the position vector (away from the central body), the z-axis along the orbital angular momentum vector, and the y-axis completing the right-handed system. NOTE: a rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`),

Parameters
----------
body : str
    Body for which the rotation matrix is to be saved.
central_body : str
    Body with respect to which the TNW frame is determined.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "inertial_to_body_fixed_313_euler_angles" && variant==0) {
        return R"(
        
Function to add the 3-1-3 Euler angles for the rotation from inertial to body-fixed frame to the dependent variables to save.

Function to add the 3-1-3 Euler angles for the rotation from inertial to body-fixed frame to the dependent variables to save. This requires the rotation of the body to be defined (either in the environment or the state vector).

Parameters
----------
body : str
    Body for which the rotation angles are to be saved.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "intermediate_aerodynamic_rotation_matrix_variable" && variant==0) {
        return R"(
        
Function to add the rotation matrix between any two reference frames used in aerodynamic calculations.

Function to add the rotation matrix between any two reference frames used in aerodynamic calculations. The list of available frames is defined by the :class:`AerodynamicsReferenceFrames` enum. NOTE: a rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`),

Parameters
----------
body : str
    Body whose dependent variable should be saved.
base_frame : AerodynamicsReferenceFrames
    Base reference frame for the rotation.
target_frame : AerodynamicsReferenceFrames
    Target reference frame for the rotation.
central_body : str
    Central body w.r.t. which the state of the body is considered.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "periapsis_altitude" && variant==0) {
        return R"(
        
Function to add the altitude of periapsis to the dependent variables to save.

Function to add the periapsis altitude of the current osculating orbit to the dependent variables to save. The altitude depends on the shape of the central body. This function takes the current (osculating) orbit of the body w.r.t. the central body, and uses this Kepler orbit to extract the position/altitude of periapsis.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the altitude of periapsis is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "apoapsis_altitude" && variant==0) {
        return R"(
        
Function to add the altitude of apoapsis to the dependent variables to save.

Function to add the apoapsis altitude of the current osculating orbit to the dependent variables to save. The altitude depends on the shape of the central body. This function takes the current (osculating) orbit of the body w.r.t. the central body, and uses this Kepler orbit to extract the position/altitude of apoapsis.

Parameters
----------
body : str
    Body whose dependent variable should be saved.
central_body : str
    Body with respect to which the altitude of apoapsis is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "central_body_fixed_spherical_position" && variant==0) {
        return R"(
        
Function to add the spherical, body-fixed position to the dependent variables to save.

Function to add the spherical position to the dependent variables to save. The spherical position is return as the radius, latitude, longitude, defined in the body-fixed frame of the central body

Parameters
----------
body : str
    Body whose spherical position is to be saved.
central_body : str
    Body with respect to which the spherical, body-fixed is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "central_body_fixed_cartesian_position" && variant==0) {
        return R"(
        
Function to add the relative Cartesian position, in the central body's fixed frame, to the dependent variables to save.


Parameters
----------
body : str
    Body whose relative cartesian position is to be saved.
central_body : str
    Body with respect to which the cartesian, body-fixed is computed.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "body_mass" && variant==0) {
        return R"(
        
Function to add the current body mass to the dependent variables to save.


Parameters
----------
body : str
    Body whose mass should be saved.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "radiation_pressure_coefficient" && variant==0) {
        return R"(
        
Function to add the current radiation pressure coefficient to the dependent variables to save.


Parameters
----------
body : str
    Body whose dependent variable should be saved.
emitting_body : str
    Emitting body.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "total_mass_rate" && variant==0) {
        return R"(
        
Function to add the total mass rate to the dependent variables to save.

Function to add the total mass rate to the dependent variables to save. It requires the body mass to be numerically propagated.

Parameters
----------
body : str
    Body whose mass rate should be saved.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "gravity_field_potential" && variant==0) {
        return R"(
        
Function to add the gravitational potential to the dependent variables to save.

Function to add the gravitational potential to the dependent variables to save. The gravitational potential is defined by the bodies undergoing and exerting the acceleration.

Parameters
----------
body_undergoing_acceleration : str
    Body whose dependent variable should be saved.
body_exerting_acceleration : str
    Body exerting acceleration.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "gravity_field_laplacian_of_potential" && variant==0) {
        return R"(
        
Function to add the laplacian of the gravitational potential to the dependent variables to save.

Function to add the laplacian of the gravitational potential to the dependent variables to save. The laplacian is defined by the bodies undergoing and exerting the acceleration.

Parameters
----------
body_undergoing_acceleration : str
    Body whose dependent variable should be saved.
body_exerting_acceleration : str
    Body exerting acceleration.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "minimum_body_distance" && variant==0) {
        return R"(
        
Function to compute the minimum distance between a given body, and a set of other bodies.

Function to compute the minimum distance between a given body, and a set of other bodies. This function takes the instantaneous position of body ``body_name``, and each body in the list ``bodies_to_check``, and computes the body from this list closest to ``body_name``. In this calculation, the positions of the bodies are evaluated at the current propagation time, and therefore **light time is ignored**. In addition, this functions does not consider visbility requirements (e.g. is a planet between two bodies). The dependent variable is of size 2, and consists of: (0) The distance between the body, and the closest other body; (1) The index from ``bodies_to_check`` for which the distance (given by the first index) is closest to ``body`` Typically, this function is used to compute the closest body in a constellation of satellites. 

Parameters
----------
body_name : str
    Body for which the distance to other bodies is to be computed.
bodies_to_check : list[ str ]
    List of bodies for which it is to be checked which of these bodies is closest to ``body_name``.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "minimum_visible_station_body_distances" && variant==0) {
        return R"(
        
Function to compute the minimum distance between a ground station, and a set of other bodies visible from that station.

Function to compute the minimum distance between a ground station, and a set of other bodies visible from that station This function takes the instantaneous position of the ground station ``station_name`` on ``body_name``, and each body in the list ``bodies_to_check``, and computes the body from this list closest to this ground station, only taking into account those bodies from this list which are visible from teh ground station. For this function, visibility is defined by a single elevation angle cutoff (at the ground station) below which a body is deemed to not be visible. In this calculation, the positions of the bodies are evaluated at the current propagation time, and therefore **light time is ignored**. The dependent variable is of size 3, and consists of: (0) The distance between the ground station, and the closest visible body; (1) The index from ``bodies_to_check`` for which the distance (given by the first index) is closest to thee ground station, and the body is visible. (2) Elevation angle for closest body. In case, no body is visible from the station, this function returns [NaN, -1, NaN]. Typically, this function is used to compute the closest body between a ground staion and a constellation of satellites. 

Parameters
----------
body_name : str
    Body on which ground station is located, for which the distance to other bodies is to be computed.
station_name : str
    Name of ground station, for which the distance to other bodies is to be computed.
bodies_to_check : list[ str ]
    List of bodies for which it is to be checked which of these bodies is closest to ``station_name`` on ``body_name``.
minimum_elevation_angle : float
    Minimum elevation angle (at ground station) below which the distance to the ``bodies_to_check`` is not considered.
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "custom_dependent_variable" && variant==0) {
        return R"(
        
Function to compute a custom dependent variable.

Function to compute a custom dependent variable, which can be implemented by the user as a Python function. The custom dependent variable is typically dependent on the current properties of the environment (e.g. bodies in the environment) or a user-defined guidance class (or similar)

Parameters
----------
custom_function : Callable[[], numpy.ndarray].
    Function taking no input, and returning the custom dependent variable (as a numpy Nx1 array).
variable_size : int
    Size N of the array returned by the ``custom_function``
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "received_irradiance" && variant==0) {
        return R"(
        
Function to save the received irradiance from a give source.

Function to save the received irradiance (in W/m^{2}) from a given source, as used in the computation of radiation pressure

Parameters
----------
target_body : str
    Name of body at the location of which the irradiance is to be returned
source_body : str
    Name of body from which the radiation originates that causes the irradiance at the target body
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else if(name == "received_irradiance_shadow_function" && variant==0) {
        return R"(
        
Function to save the shadow function that reduces the received irradiance from a given source.

Function to save the shadow function that reduces the received irradiance from a given source, as a result of occultation of radiation between the source and the target. The shadow function is by definition between 0 (no irradiance is received, total occultation) and 1 (all irradiance is received, no occultation)

Parameters
----------
target_body : str
    Name of body at the location of which the irradiance is to be returned
source_body : str
    Name of body from which the radiation originates that causes the irradiance at the target body
Returns
-------
SingleDependentVariableSaveSettings
    Dependent variable settings object.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace integrator {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "AvailableIntegrators") {
         return R"(

        Enumeration of integrators available with tudat.





     )";


    } else if(name == "AvailableIntegrators.runge_kutta_fixed_step_size_type") {
         return R"(
     )";


    } else if(name == "AvailableIntegrators.runge_kutta_variable_step_size_type") {
         return R"(
     )";


    } else if(name == "AvailableIntegrators.bulirsch_stoer_type") {
         return R"(
     )";


    } else if(name == "AvailableIntegrators.adams_bashforth_moulton_type") {
         return R"(
     )";



    } else if(name == "CoefficientSets") {
         return R"(

        Coefficient sets for Runge-Kutta-type integrators.

        Coefficient sets for Runge-Kutta-type integrators. The coefficients are defined
        in a Butcher Tableau, with an coefficient set yielding an x(y) method yielding an integrator
        with global truncation error of :math:`O(\Delta t^{x})`. Some of these coefficients also contain an embedded integrator of :math:`O(\Delta t^{y})`
        for step size control.





     )";


    } else if(name == "CoefficientSets.euler_forward") {
         return R"(
Coefficients for the classic forward Euler method
     )";


    } else if(name == "CoefficientSets.rk_4") {
         return R"(
Coefficients for the original Runge-Kutta method of order 4
     )";


    } else if(name == "CoefficientSets.explicit_mid_point") {
         return R"(
Coefficients for the explicit midpoint method
     )";


    } else if(name == "CoefficientSets.explicit_trapezoid_rule") {
         return R"(
Coefficients for the explicit trapezoid rule, also called Heun's method or improved Euler's method
     )";


    } else if(name == "CoefficientSets.ralston") {
         return R"(
Coefficients for Ralston's method
     )";


    } else if(name == "CoefficientSets.rk_3") {
         return R"(
Coefficients for the Runge-Kutta method of order 3
     )";


    } else if(name == "CoefficientSets.ralston_3") {
         return R"(
Coefficients for Ralston's third-order method
     )";


    } else if(name == "CoefficientSets.SSPRK3") {
         return R"(
Coefficients for the Strong Stability Preserving Runge-Kutta third-order method
     )";


    } else if(name == "CoefficientSets.ralston_4") {
         return R"(
Coefficients for Ralston's fourth-order method
     )";


    } else if(name == "CoefficientSets.three_eight_rule_rk_4") {
         return R"(
Coefficients for the classic Runge Kutta 3/8-rule fourth-order method
     )";


    } else if(name == "CoefficientSets.rkf_12") {
         return R"(
Coefficients for the Runge-Kutta-Fehlberg method of order 2 with an embedded 1st order
     )";


    } else if(name == "CoefficientSets.heun_euler") {
         return R"(
Coefficients for the Heun's method of order 2 with an embedded Euler method of order 1
     )";


    } else if(name == "CoefficientSets.rkf_45") {
         return R"(
Coefficients for the Runge-Kutta-Fehlberg method of order 5 with an embedded 4th order
     )";


    } else if(name == "CoefficientSets.rkf_56") {
         return R"(
Coefficients for the Runge-Kutta-Fehlberg method of order 6 with an embedded 5th order
     )";


    } else if(name == "CoefficientSets.rkf_78") {
         return R"(
Coefficients for the Runge-Kutta-Fehlberg method of order 8 with an embedded 7th order
     )";


    } else if(name == "CoefficientSets.rkdp_87") {
         return R"(
Coefficients for the Dormand-Prince method of order 7 with an embedded 8th order
     )";


    } else if(name == "CoefficientSets.rkf_89") {
         return R"(
Coefficients for the Runge-Kutta-Fehlberg method of order 9 with an embedded 8th order
     )";


    } else if(name == "CoefficientSets.rkv_89") {
         return R"(
Coefficients for the Runge-Kutta-Verner method of order 9 with an embedded 8th order
     )";


    } else if(name == "CoefficientSets.rkf_108") {
         return R"(
Coefficients for the Runge-Kutta-Feagin method of order 8 with an embedded 10th order
     )";


    } else if(name == "CoefficientSets.rkf_1210") {
         return R"(
Coefficients for the Runge-Kutta-Feagin method of order 10 with an embedded 12ve order
     )";


    } else if(name == "CoefficientSets.rkf_1412") {
         return R"(
Coefficients for the Runge-Kutta-Feagin method of order 12 with an embedded 14th order
     )";



    } else if(name == "OrderToIntegrate") {
         return R"(

        Enumeration defining Which integrator order needs to be integrated, only used for coefficient sets with an embedded order.





     )";


    } else if(name == "OrderToIntegrate.lower") {
         return R"(
For a method of order :math:`p`, with embedded method of order :math:`q`, the step is taken using the method with order :math:`\min(p,q)`
     )";


    } else if(name == "OrderToIntegrate.higher") {
         return R"(
For a method of order :math:`p`, with embedded method of order :math:`q`, the step is taken using the method with order :math:`\max(p,q)`
     )";



    } else if(name == "ExtrapolationMethodStepSequences") {
         return R"(

        Enumeration of available extrapolation method substep sequences, with :math:`n_{j}` defining the number of substeps in iteration :math:`j`.





     )";


    } else if(name == "ExtrapolationMethodStepSequences.bulirsch_stoer_sequence") {
         return R"(
Sequence for which :math:`n_{j}=2n_{j-2}` (2, 4, 6, 8, 12, 16, 24, ....)
     )";


    } else if(name == "ExtrapolationMethodStepSequences.deufelhard_sequence") {
         return R"(
Sequence for which :math:`n_{j}=2(j+1)` (2, 4, 6, 8, 10, 12, 14, ....)
     )";



    } else if(name == "MinimumIntegrationTimeStepHandling") {
         return R"(

        Enumeration defining possible behaviours when :math:`\Delta t_{rec}<\Delta t_{\min}`. in step-size control (e.g. recommended time step is smaller than minimum time step)





     )";


    } else if(name == "MinimumIntegrationTimeStepHandling.throw_exception_below_minimum") {
         return R"(
Exception is throw, and propagation is terminated
     )";


    } else if(name == "MinimumIntegrationTimeStepHandling.set_to_minimum_step_silently") {
         return R"(
The final time step is set to :math:`\Delta t=\Delta t_{\min}`, violating requirements of step-size control algorithm, without any message to user"
     )";


    } else if(name == "MinimumIntegrationTimeStepHandling.set_to_minimum_step_single_warning") {
         return R"(
The final time step is set to :math:`\Delta t=\Delta t_{\min}`, violating requirements of step-size control algorithm, a warning is printed to the terminal the first time this happens during a propagation"
     )";


    } else if(name == "MinimumIntegrationTimeStepHandling.set_to_minimum_step_every_time_warning") {
         return R"(
The final time step is set to :math:`\Delta t=\Delta t_{\min}`, violating requirements of step-size control algorithm, a warning is printed to the terminal every time this happens during a propagation"
     )";




    } else if(name == "IntegratorStepSizeControlSettings") {
         return R"(

        Base class to define settings for step-size control algorithm.

        Base class to define settings for step-size control algorithm, typically created by one of the factory functions provided in this module





     )";





    } else if(name == "IntegratorStepSizeValidationSettings") {
         return R"(

        Base class to define settings for step-size validation algorithm.

        Base class to define settings for step-size validation algorithm, typically created by one of the factory functions provided in this module





     )";





    } else if(name == "IntegratorSettings") {
         return R"(

        Functional base class to define settings for integrators.

        Class to define settings for numerical integrators, for instance for use in numerical integration of equations of motion/
        variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that
        require more settings to define have their own derived class.





     )";





    } else if(name == "RungeKuttaFixedStepSizeSettings") {
         return R"(

        `IntegratorSettings`-derived class to define settings for Runge Kutta integrators with a fixed step size





     )";





    } else if(name == "BulirschStoerIntegratorSettings") {
         return R"(

        `IntegratorSettings`-derived class to define settings for Bulirsch-Stoer integrator settings.





     )";





    } else if(name == "AdamsBashforthMoultonSettings") {
         return R"(

        `IntegratorSettings`-derived class to define settings for Adams-Bashforth-Moulton integrator settings.





     )";






    } else if(name == "step_size_validation" && variant==0) {
        return R"(
        
Creates settings step size validation in a variable step-size integrator.

Factory function to create settings step size validation in a variable step-size integrator. The validation
model takes the proposed new step size  :math:`\Delta t_{rec}` as input, and checks if it meets predefined conditions, specifically
whether the proposed time step falls in a given predefined range :math:`[\Delta t_{\min}, \Delta t_{\max}]`.
This function also provides the option of handling recommended step sizes below :math:`\Delta t_{\min}` in various ways,
and control on how to deal with recommend Inf/NaN step sizes.


Parameters
----------
minimum_step : float
    Value of minimum permitted time step :math:`\Delta t_{\min}`.
maximum_step : float
    Value of maximum permitted time step :math:`\Delta t_{\max}`.
minimum_step_size_handling : MinimumIntegrationTimeStepHandling, default = throw_exception_below_minimum
    Entry defining the behaviour when :math:`\Delta t_{rec}<\Delta t_{\min}`.
accept_infinity_step : bool, default = False
    Entry defining whether to accept a step size of infinity (if False, exception is throw in such cases)
accept_nan_step : bool, default = False
    Entry defining whether to accept a step size of NaN (if False, exception is throw in such cases)
Returns
-------
IntegratorStepSizeValidationSettings
    Object containing settings for step-size validation.






    )";



    } else if(name == "step_size_control_elementwise_scalar_tolerance" && variant==0) {
        return R"(
        
Creates settings for integrator step-size control, using element-wise analysis for the propagated states.

Function to create settings for integrator step-size control, using element-wise analysis for the propagated states. For a propagated
state :math:`\mathbf{x}` with entries :math:`x_{i}`, and an estimate :math:`\boldsymbol{\epsilon}` for the current local error, the following
algorithm is performed per element :math:`i` to calculate the required error :math:`\epsilon_{i,req}` on this element:

.. math::

   \epsilon_{i,req}=\epsilon_{r}x_{i}+\epsilon_{a}
   
A proposed modification to the step size is then computed, using the most constraining of all state elements

.. math::
   
   \bar{\Delta t_{rec.}}&=\Delta t\left(\min_{i}\left(\frac{\epsilon_{i,req}}{\epsilon_{i}}\right)\right)^{p}\\
   \Delta t_{rec.}&=K\bar{\Delta t_{rec.}}

with :math:`p` the order of the local truncation error of the method for which step-size control is being applied, 
:math:`\Delta t_{rec.}` the new, recommended step size, and :math:`\Delta t` the current step size. The factor :math:`K` is a safety factor
used make the time step slightly smaller than strictly required.

A minimum and maximum change in time step may be provided by the user, such that if :math:`\Delta t_{rec.}/\Delta t` is too large or too small,
the proposed increase/decrease to the step size is constrained to this limit value. That is, if :math:`\Delta t_{rec.}/\Delta t` proposed by the algorithm is
1000, and the ``maximum_factor_increase`` input is equal to 20, the algorithm will use :math:`\Delta t_{rec.}/\Delta t=20` in what follows. 

For cases where :math:`\bar{\Delta t_{rec.}}/\Delta t < 1`, the step is recommended to be recomputed with the new proposed step size (e.g. the current step 
is not accepted, and will be re-attempted with a smaller step size). For cases where :math:`\bar{\Delta t_{rec.}}/\Delta t > 1`, the step is accepted, and 
the next step will be performed with the new, higher, step size.


Parameters
----------
relative_error_tolerance : float
    Value of relative error tolerance :math:`\epsilon_{r}`.
absolute_error_tolerance : float
    Value of absolute error tolerance :math:`\epsilon_{a}`.
safety_factor : float, default = 0.8
    Safety factor :math:`K` for step size control
minimum_factor_increase : float, default = 0.1
    Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
maximum_factor_increase : float, default = 4.0
    Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
Returns
-------
IntegratorStepSizeControlSettings
    Object containing settings for per-element step-size control.






    )";



    } else if(name == "step_size_control_elementwise_matrix_tolerance" && variant==0) {
        return R"(
        
Creates settings for integrator step-size control, using element-wise analysis for the propagated states.

Function to create settings for integrator step-size control, using element-wise analysis for the propagated states. This function
is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`,
with the differences that the tolerances are provided as a vector/matrix (which must be of equal size as the propagates state), such that
different tolerances can be provided for each state element. The behaviour of the algorithm is then such that 
:math:`\epsilon_{r}\rightarrow\epsilon_{r,i}` and :math:`\epsilon_{a}\rightarrow\epsilon_{a,i}`. 

If the size of the tolerances used as input differ from one another, or differ from the size of the state vector, an exception is thrown


Parameters
----------
relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
    Values of relative error tolerance :math:`\boldsymbol{\epsilon}_{r}`.
absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
    Values of absolute error tolerance :math:`\boldsymbol{\epsilon}_{a}`.
safety_factor : float, default = 0.8
    Safety factor :math:`K` for step size control
minimum_factor_increase : float, default = 0.1
    Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
maximum_factor_increase : float, default = 4.0
    Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
Returns
-------
IntegratorStepSizeControlSettings
    Object containing settings for per-element step-size control.






    )";



    } else if(name == "step_size_control_blockwise_scalar_tolerance" && variant==0) {
        return R"(
        
Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`,
with the difference that the error estimation :math:`\boldsymbol{\epsilon}` is not used on an element-by-element basis, but using the norms
of user defined matrix blocks. This is for instance very useful when propagating Cartesian states, where the tolerances are then
typically applied twice: once to the norm of the position error, and once to the norm of the velocity error.

The algorithm is then run, using the modification 
that :math:`\epsilon_{i}\rightarrow||\boldsymbol{\epsilon_{[i,k],[j,l]}}||`. Where the indices on the right-hand side denote start row :math:`i`,
start column :math:`j`, number of rows :math:`k` and number of columns :math:`l`. over
which the state error norm is to be taken. For a single Cartesian state vector, the norm is taken on blocks :math:`[0,3],[0,1]` and :math:`[3,3],[0,1]`
            


Parameters
----------
block_indices : list[tuple[int,int,int,int]]
    List of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order.
relative_error_tolerance : float
    Value of relative error tolerance :math:`\epsilon_{r}`.
absolute_error_tolerance : float
    Value of absolute error tolerance :math:`\epsilon_{a}`.
safety_factor : float, default = 0.8
    Safety factor :math:`K` for step size control
minimum_factor_increase : float, default = 0.1
    Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
maximum_factor_increase : float, default = 4.0
    Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
Returns
-------
IntegratorStepSizeControlSettings
    Object containing settings for per-element step-size control.






    )";



    } else if(name == "step_size_control_blockwise_matrix_tolerance" && variant==0) {
        return R"(
        
Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance`,
with the differences that the tolerances are provided as a list (which must be of equal size as the number of state blocks used), such that
different tolerances can be provided for each state block. 

If the size of the tolerances used as input differ from one another, or differ from the number of blocks, an exception is thrown
            


Parameters
----------
block_indices : list[tuple[int,int,int,int]]
    List of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order.
relative_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
    Values of relative error tolerance :math:`\boldsymbol{\epsilon}_{r}`.
absolute_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
    Values of absolute error tolerance :math:`\boldsymbol{\epsilon}_{a}`.
safety_factor : float, default = 0.8
    Safety factor :math:`K` for step size control
minimum_factor_increase : float, default = 0.1
    Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
maximum_factor_increase : float, default = 4.0
    Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
Returns
-------
IntegratorStepSizeControlSettings
    Object containing settings for per-element step-size control.






    )";



    } else if(name == "step_size_control_custom_blockwise_scalar_tolerance" && variant==0) {
        return R"(
        
Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance`,
but rather than providing the ``block_indices`` directly, a function to determine the block indices, based on the size of the 
propagated state, is provided. For instance, the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.standard_cartesian_state_element_blocks`
can be provided to this function (as ``block_indices_function``), which will adapt the block indices depending on the size of the propagated state
(e.g. regardless of how many bodies are propagated, step size control will always be done on position and velocity element blocks)  

            


Parameters
----------
block_indices_function : Callable[[int,int],list[tuple[int,int,int,int]]]
    Function returning list of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order, with number of rows and columns of propagated state as input.
relative_error_tolerance : float
    Value of relative error tolerance :math:`\epsilon_{r}`.
absolute_error_tolerance : float
    Value of absolute error tolerance :math:`\epsilon_{a}`.
safety_factor : float, default = 0.8
    Safety factor :math:`K` for step size control
minimum_factor_increase : float, default = 0.1
    Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
maximum_factor_increase : float, default = 4.0
    Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
Returns
-------
IntegratorStepSizeControlSettings
    Object containing settings for per-element step-size control.






    )";



    } else if(name == "step_size_control_custom_blockwise_matrix_tolerance" && variant==0) {
        return R"(
        
Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance`,
but uses blockwise tolerances (as in :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_blockwise_matrix_tolerance`)     
           


Parameters
----------
block_indices_function : Callable[[int,int],list[tuple[int,int,int,int]]]
    Function returning list of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order, with number of rows and columns of propagated state as input.
relative_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
    Values of relative error tolerance :math:`\boldsymbol{\epsilon}_{r}`.
absolute_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
    Values of absolute error tolerance :math:`\boldsymbol{\epsilon}_{a}`.
safety_factor : float, default = 0.8
    Safety factor :math:`K` for step size control
minimum_factor_increase : float, default = 0.1
    Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
maximum_factor_increase : float, default = 4.0
    Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
Returns
-------
IntegratorStepSizeControlSettings
    Object containing settings for per-element step-size control.






    )";



    } else if(name == "standard_cartesian_state_element_blocks" && variant==0) {
        return R"(
        
Function to generate step size control blocks on position and velocity elements for numerical integration

Function to generate step size control blocks on position and velocity elements for numerical integration, typically provided
to the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance` or
:func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_matrix_tolerance` function.
By providing this function to one of these step-size control functions, the final column of the state vector is taken (such that
it works  both for state-only, and variational equations and state propagation) and combined into :math:`N` blocks of size 3.
The step-size control is then done on each of these blocks, which will represent the position and velocity blocks.


Parameters
----------
number_of_rows : int
    Number of rows in state vector
number_of_columns : int
    Number of columns in state vector
Returns
-------
list[tuple[int,int,int,int]]
    List of matrix blocks over which the step size control is to be done (see ``block_indices_function`` input to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance`)






    )";



    } else if(name == "runge_kutta_fixed_step" && variant==0) {
        return R"(
        
Creates the settings for the Runge-Kutta fixed step size integrator.

Factory function to create settings for the Runge-Kutta integrator with a constant step size.
Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


Parameters
----------
time_step : float
    Initial time step to be used.
coefficient_set : CoefficientSets
    Coefficient set (Butcher's tableau) to be used in the integration.
order_to_use : OrderToIntegrate, default=OrderToIntegrate.lower
    If the coefficient set is supposed to be for variable step sizes (with an embedded method of a different order),
    this parameter can be used to set the order that will be used.

assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
    integrator (true) or only at the end of each integration step (false).

Returns
-------
IntegratorSettings
    Object containing settings for the integrator.





Examples
--------
In this example, settings for the classical RK4 integrator with 30 second time step are created

.. code-block:: python

  # Create RK4 settings
  integrator_settings = integrator.runge_kutta_fixed_step(
      time_step = 30.0,
      coefficient_set = integrator.rk_4 )

In this example, settings for fixed-step integration using the higher-order (8th-order) of the two
embedded propagators of the RKF7(8) method are created, with a time-step of 120 seconds.

.. code-block:: python

  # Create 8th order RKF settings
  integrator_settings = integrator.runge_kutta_fixed_step(
      time_step = 120.0,
      coefficient_set = integrator.rkf_78,
      order_to_use = integrator.higher )      
  


    )";



    } else if(name == "runge_kutta_variable_step" && variant==0) {
        return R"(
        
Creates the settings for the Runge-Kutta fixed step size integrator.

Factory function to create settings for the Runge-Kutta variable step size integrator.
Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).
The step-size control algorithm is defined by a :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeControlSettings` and
:class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeValidationSettings` object, created using one of the factory functions
listed above.


Parameters
----------
initial_time_step : float
    Initial time step to be used.
coefficient_set : CoefficientSets
    Coefficient set (Butcher's tableau) to be used in the integration.
step_size_control_settings : IntegratorStepSizeControlSettings
    Object used to control the step size, by computing a new step size :math:`\Delta t_{rec.}`, from the embedded Runge-Kutta integrator pair,
    and recommending whether the steps is to be accepted, or recomputed with a different time step.

step_size_validation_settings : IntegratorStepSizeValidationSettings
    Object used to validate whether the :math:`\Delta t_{rec.}` provided by model defined by the ``step_size_control_settings`` meets with user-defined
    criteria (minimum, maximum values, etc.)

assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
    integrator (true) or only at the end of each integration step (false).

Returns
-------
IntegratorSettings
    Object containing settings for the integrator.





Examples
--------
In this example, settings for the varianle step RK4(5) integrator are created, with the same tolerances (:math:`10^{-10}`)
applied element-wise on the propagated state. The minimum and maximum time steps are set to 0.001 and 1000 seconds,
and the initial step is set to 30 seconds. All other inputs are left on their defaults

.. code-block:: python

  # Create RK4(5) settings
  control_settings = integator.step_size_control_elementwise_scalar_tolerance( 1.0E-10, 1.0E-10 )
  validation_settings = integrator.step_size_validation( 0.001, 1000.0 )
  integrator_settings = integrator.runge_kutta_variable_step(
      initial_time_step = 30.0,
      coefficient_set = integrator.rkf_45,
      step_size_control_settings = control_settings,
      step_size_validation_settings = validation_settings )

In this example, the above is modified such that step-size control is applied on position and velocity
element blocks.

.. code-block:: python

  # Create RK4(5) settings
  control_settings = integator.step_size_control_custom_blockwise_scalar_tolerance( 
      integrator.standard_cartesian_state_element_blocks
      1.0E-10, 1.0E-10 )
  validation_settings = integrator.step_size_validation( 0.001, 1000.0 )
  integrator_settings = integrator.runge_kutta_variable_step(
      initial_time_step = 30.0,
      coefficient_set = integrator.rkf_45,
      step_size_control_settings = control_settings,
      step_size_validation_settings = validation_settings )


    )";



    } else if(name == "bulirsch_stoer_fixed_step" && variant==0) {
        return R"(
        
Creates the settings for the fixed time-step Bulirsch-Stoer integrator.

Factory function to create settings for the fixed time-step Bulirsch-Stoer integrator. The
underlying method is the same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.bulirsch_stoer_variable_step`, 
but using a fixed, user-defined, time step.


Parameters
----------
time_step : float
    Time step to be used.
extrapolation_sequence : ExtrapolationMethodStepSequences
    Extrapolation sequence to be used for the integration (defining the number of substeps in iteration :math:`i`).
maximum_number_of_steps : int
    Number of entries from the sequence to be used (e.g., total number of iterations used for a single extrapolation and time step).
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.





Examples
--------
In this example, settings for the Bulirsch-Stoer integrator with 600 second time step are created, using the typical
sequence, using 6 iterations of the same step. By using the Bulirsch-Stoer sequence, this means that the same step is done
using 2, 4, 6, 8, 12 and 16 substeps     

.. code-block:: python

  # Create BS settings
  integrator_settings = integrator.bulirsch_stoer_fixed_step(
      time_step = 300.0,
      extrapolation_sequence = integrator.bulirsch_stoer_sequence,
      maximum_number_of_steps = 6 )


    )";



    } else if(name == "bulirsch_stoer_variable_step" && variant==0) {
        return R"(
        
Creates the settings for the variable time-step Bulirsch-Stoer integrator.

Factory function to create settings for the variable time-step Bulirsch-Stoer integrator. This integrator
works by performing the same (typically very large) step multiple times, using an ever increasing number of substeps.
Each substep is performed using the modified midpoint method. The succesive integrations from :math:`t_{i}` to :math:`t_{i}+\Delta t`
are (in principle) done using ever increasing accuracy, as the size of the substep decreases. This integrator works
by extrapolating the behaviour to a substep length of 0 (e.g. an infinite number of substeps), at which the solution should be perfect.
The number of substeps on the :math:`i^{t}` iteration are done using the number of substeps defined by  entry :math:`i` of the
``extrapolation_sequence`` input. The number of iterations for a single step is defined by the ``maximum_number_of_steps`` entry.
For instance, using the ``bulirsch_stoer_sequence`` sequence, and 5 iterations, the same step is done using 2, 4, 6, 8 and 12 substeps,
and the results are then extrapolated to an infinite number of steps. Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).

The step-size control algorithm is defined by a :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeControlSettings` and
:class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeValidationSettings` object, created using one of the factory functions
listed above. The time step control uses the result from the final, and second to final iteration to generate an error estimate of the current step.


Parameters
----------
time_step : float
    Initial time step to be used.
extrapolation_sequence : ExtrapolationMethodStepSequences
    Extrapolation sequence to be used for the integration (defining the number of substeps in iteration :math:`i`).
maximum_number_of_steps : int
    Number of entries from the sequence to be used (e.g., total number of iterations used for a single extrapolation and time step).
step_size_control_settings : IntegratorStepSizeControlSettings
    Object used to control the step size, by computing a new step size :math:`\Delta t_{rec.}`, from the embedded Runge-Kutta integrator pair,
    and recommending whether the steps is to be accepted, or recomputed with a different time step.

step_size_validation_settings : IntegratorStepSizeValidationSettings
    Object used to validate whether the :math:`\Delta t_{rec.}` provided by model defined by the ``step_size_control_settings`` meets with user-defined
    criteria (minimum, maximum values, etc.)

assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.





Examples
--------
In this example, settings for the Bulirsch-Stoer integrator with 600 second initial time step are created, using the typical
sequence, using 6 iterations of the same step. By using the Bulirsch-Stoer sequence, this means that the same step is done
using 2, 4, 6, 8, 12 and 16 substeps. The same tolerances (:math:`10^{-10}`)
applied element-wise on the propagated state. The minimum and maximum time steps are set to 0.1 and 10000 seconds,
and the initial step is set to 600 seconds. All other inputs are left on their defaults

.. code-block:: python

  # Create BS settings
  control_settings = integator.step_size_control_elementwise_scalar_tolerance( 1.0E-10, 1.0E-10 )
  validation_settings = integrator.step_size_validation( 0.1, 10000.0 )
  integrator_settings = integrator.bulirsch_stoer_variable_step(
      initial_time_step = 600.0,
      extrapolation_sequence = integrator.bulirsch_stoer_sequence,
      maximum_number_of_steps = 6 
      step_size_control_settings = control_settings,
      step_size_validation_settings = validation_settings )


    )";



    } else if(name == "adams_bashforth_moulton" && variant==0) {
        return R"(
        
Creates the settings for the Adams-Bashforth-Moulton integrator.

Factory function to create settings for the Adams-Bashforth-Moulton multistep integrator.
For this integrator, the step size and order are both according to a control algorithm
similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`.      
The integrator is initialized using an RKF7(8) integrator.

NOTE: this integrator's step-size and order control algorithm work in a method that is overly simplistic,
when increasing/decreasing the order, existing function evaluations are re-used, without any recomputations.
Similarly, when halving or doubling the time-step, the existing interpolating polynomial is evaluated at the relevant points.
This can lead to unwanted behaviour, where the time-step reduces to unrealistically low values. It is strongly
recommended that a reasonable minimum step is provided to this function, to partially mitigate this behaviour.


Parameters
----------
initial_time_step : float
    Initial time step to be used.
minimum_step_size : float
    Minimum time step to be used during the integration.
maximum_step_size : float
    Maximum time step to be used during the integration.
relative_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
absolute_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
minimum_order
    Minimum order of the integrator.
maximum_order
    Maximum order of the integrator.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
bandwidth : float, default=200.0
    Maximum error factor for doubling the step size.
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.






    )";



    } else if(name == "adams_bashforth_moulton_fixed_order" && variant==0) {
        return R"(
        
Creates the settings for the Adams-Bashforth-Moulton integrator of fixed order.

Same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.adams_bashforth_moulton`, but
with fixed order and variable step


Parameters
----------
initial_time_step : float
    Initial time step to be used.
minimum_step_size : float
    Minimum time step to be used during the integration.
maximum_step_size : float
    Maximum time step to be used during the integration.
relative_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
absolute_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
order
    Order of the integrator.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
bandwidth : float, default=200.0
    Maximum error factor for doubling the step size.
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.






    )";



    } else if(name == "adams_bashforth_moulton_fixed_step" && variant==0) {
        return R"(
        
Creates the settings for the Adams-Bashforth-Moulton fixed-step integrator.

Same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.adams_bashforth_moulton`, but
with fixed step and variable order


Parameters
----------
time_step : float
    Initial time step to be used.
relative_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
absolute_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
minimum_order
    Minimum order of the integrator.
maximum_order
    Maximum order of the integrator.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
bandwidth : float, default=200.0
    Maximum error factor for doubling the step size.
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.






    )";



    } else if(name == "adams_bashforth_moulton_fixed_step_fixed_order" && variant==0) {
        return R"(
        
Creates the settings for the Adams-Bashforth-Moulton fixed-step, fixed-order integrator.

Same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.adams_bashforth_moulton`, but
with fixed step and fixed order


Parameters
----------
time_step : float
    Initial time step to be used.
order
    Order of the integrator.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.






    )";



    } else if(name == "print_butcher_tableau" && variant==0) {
        return R"(
        
Print the Butcher tableau of a given coefficient set.


Parameters
----------
coefficient_set : CoefficientSets
    Coefficient set of which the Butcher tableau will be printed.





    )";



    } else if(name == "runge_kutta_variable_step_size" && variant==0) {
        return R"(
        
Creates the settings for the Runge-Kutta variable step size integrator with scalar tolerances.

NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.runge_kutta_variable_step` INTERFACE INSTEAD
    
Factory function to create settings for the Runge-Kutta variable step size integrator with scalar tolerances.
For this integrator, the step size is varied based on the tolerances and safety factor provided.
The tolerance is composed of an absolute and a relative part.
Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


Parameters
----------
initial_time_step : float
    Initial time step to be used.
coefficient_set : CoefficientSets
    Coefficient set (Butcher's tableau) to be used in the integration.
minimum_step_size : float
    Minimum time step to be used during the integration.
maximum_step_size : float
    Maximum time step to be used during the integration.
relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
    Relative vector tolerance to adjust the time step.
absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
    Absolute vector tolerance to adjust the time step.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
    integrator (true) or only at the end of each integration step (false).

safety_factor : float, default=0.8
    Safety factor used in the step size control.
maximum_factor_increase : float, default=4.0
    Maximum increase between consecutive time steps, expressed as the factor between new and old step size.

minimum_factor_increase : float, default=0.1
    Minimum increase between consecutive time steps, expressed as the factor between new and old step size.

throw_exception_if_minimum_step_exceeded : bool, default=true
    If set to false, the variable step integrator will use the minimum step size specified when the algorithm
    computes the optimum one to be lower, instead of throwing an exception.

Returns
-------
RungeKuttaVariableStepSettingsScalarTolerances
    RungeKuttaVariableStepSettingsScalarTolerances object.






    )";



    } else if(name == "runge_kutta_variable_step_size_vector_tolerances" && variant==0) {
        return R"(
        
Creates the settings for the Runge-Kutta variable step size integrator with vector tolerances.

NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.runge_kutta_variable_step` INTERFACE INSTEAD

Factory function to create settings for the Runge-Kutta variable step size integrator with vector tolerances.
For this integrator, the step size is varied based on the tolerances and safety factor provided.
The tolerance is composed of an absolute and a relative part.
Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


Parameters
----------
initial_time_step : float
    Initial time step to be used.
coefficient_set : CoefficientSets
    Coefficient set (Butcher's tableau) to be used in the integration.
minimum_step_size : float
    Minimum time step to be used during the integration.
maximum_step_size : float
    Maximum time step to be used during the integration.
relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
    Relative vector tolerance to adjust the time step.
absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
    Absolute vector tolerance to adjust the time step.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
    integrator (true) or only at the end of each integration step (false).

safety_factor : float, default=0.8
    Safety factor used in the step size control.
maximum_factor_increase : float, default=4.0
    Maximum increase between consecutive time steps, expressed as the factor between new and old step size.

minimum_factor_increase : float, default=0.1
    Minimum increase between consecutive time steps, expressed as the factor between new and old step size.

throw_exception_if_minimum_step_exceeded : bool, default=true
    If set to false, the variable step integrator will use the minimum step size specified when the algorithm
    computes the optimum one to be lower, instead of throwing an exception.

Returns
-------
RungeKuttaVariableStepSizeSettingsVectorTolerances
    RungeKuttaVariableStepSizeSettingsVectorTolerances object.






    )";



    } else if(name == "bulirsch_stoer" && variant==0) {
        return R"(
        
Creates the settings for the Bulirsch-Stoer integrator.


NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.bulirsch_stoer_variable_step` INTERFACE INSTEAD

Factory function to create settings for the Bulirsch-Stoer integrator.
For this integrator, the step size is varied based on the tolerances and safety factor provided.
The tolerance is composed of an absolute and a relative part.
Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).


Parameters
----------
initial_time_step : float
    Initial time step to be used.
extrapolation_sequence : ExtrapolationMethodStepSequences
    Extrapolation sequence to be used in the integration.
maximum_number_of_steps : int
    Number of entries in the sequence (e.g., number of integrations used for a single extrapolation).
minimum_step_size : float
    Minimum time step to be used during the integration.
maximum_step_size : float
    Maximum time step to be used during the integration.
relative_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
absolute_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
safety_factor : float, default=0.7
    Safety factor used in the step size control.
maximum_factor_increase : float, default=10.0
    Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
minimum_factor_increase : float, default=0.1
    Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
Returns
-------
BulirschStoerIntegratorSettings
    BulirschStoerIntegratorSettings object.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace mass_rate {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "AvailableMassRateModels") {
         return R"(

        Enumeration of available mass rate models.

        Enumeration of mass rate models supported by tudat.





     )";


    } else if(name == "AvailableMassRateModels.undefined_mass_rate_type") {
         return R"(
     )";


    } else if(name == "AvailableMassRateModels.custom_mass_rate_type") {
         return R"(
     )";


    } else if(name == "AvailableMassRateModels.from_thrust_mass_rate_type") {
         return R"(
     )";




    } else if(name == "MassRateModelSettings") {
         return R"(

        Functional base class to define settings for mass rates.

        Base class for providing settings for a mass rate model, that defines the models to be used to numerically propagate the
        mass of a body during a simulation. If any additional information (besides the type of the mass rate model) is required,
        these must be implemented in a derived class (dedicated for each mass rate model type).





     )";





    } else if(name == "FromThrustMassRateSettings") {
         return R"(

        `MassRateModelSettings`-derived class to define settings for a mass rate model derived from a thrust model.

        `MassRateModelSettings`-derived class to define settings for a mass rate model derived from a thrust model.





     )";





    } else if(name == "CustomMassRateSettings") {
         return R"(

        `MassRateModelSettings`-derived class to define settings for a custom mass rate model.

        `MassRateModelSettings`-derived class to define settings for a custom mass rate model.





     )";






    } else if(name == "from_thrust" && variant==0) {
        return R"(
        
Creates the settings for a mass rate model defined from a thrust model.

Creates the settings for a mass rate model defined from a thrust model. The mass rate model is derived from
the associated body's engine model. It is possible to consider only a specific engine or all engines.


Parameters
----------
use_all_thrust_models : bool, default=true
    Denotes whether all engines of the associated body are to be combined into a single thrust model.
associated_thrust_source : str, default=""
    Name of engine model from which thrust is to be derived (must be empty if the first argument is set to true).
Returns
-------
FromThrustMassRateSettings
    From thrust mass rate settings object.






    )";



    } else if(name == "custom_mass_rate" && variant==0) {
        return R"(
        
Creates the settings for a mass rate model defined from a thrust model.

Creates the settings for a custom mass rate model defined through a mass rate function. The function must have
time as an independent variable.


Parameters
----------
mass_rate_function : callable[[float], float]
    Function of time defining the custom mass rate.
Returns
-------
CustomMassRateSettings
    Custom mass rate settings object.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace propagator {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "TranslationalPropagatorType") {
         return R"(

        Enumeration of available translational propagator types.





     )";


    } else if(name == "TranslationalPropagatorType.undefined_translational_propagator") {
         return R"(
     )";


    } else if(name == "TranslationalPropagatorType.cowell") {
         return R"(
Propagation of Cartesian elements (state vector size 6), without any transformations
     )";


    } else if(name == "TranslationalPropagatorType.encke") {
         return R"(
Propagation of the difference in Cartesian elements of the orbit w.r.t. an unperturbed reference orbit. The reference orbit is generated from the initial state/central body, and not updated during the propagation (see Wakker, 2015 [2]_)
     )";


    } else if(name == "TranslationalPropagatorType.gauss_keplerian") {
         return R"(
Propagation of Keplerian elements (state vector size 6), with true anomaly as the 'fast' element  (see Vallado, 2001 [4]_)
     )";


    } else if(name == "TranslationalPropagatorType.gauss_modified_equinoctial") {
         return R"(
Propagation of Modified equinoctial elements (state vector size 6), with the element :math:`I` defining the location of the singularity based on the initial condition (see Hintz, 2008 [3]_)
     )";


    } else if(name == "TranslationalPropagatorType.unified_state_model_quaternions") {
         return R"(
Propagation of Unified state model using quaternions (state vector size 7, see Vittaldev et al., 2012 [1]_)
     )";


    } else if(name == "TranslationalPropagatorType.unified_state_model_modified_rodrigues_parameters") {
         return R"(
Propagation of Unified state model using modified Rodrigues parameters (state vector size 7, last element represents shadow parameter, see Vittaldev et al., 2012 [1]_)
     )";


    } else if(name == "TranslationalPropagatorType.unified_state_model_exponential_map") {
         return R"(
Propagation of Unified state model using exponential map (state vector size 7, last element represents shadow parameter, see Vittaldev et al., 2012 [1]_)
     )";



    } else if(name == "RotationalPropagatorType") {
         return R"(

        Enumeration of available rotational propagator types.





     )";


    } else if(name == "RotationalPropagatorType.undefined_rotational_propagator") {
         return R"(
     )";


    } else if(name == "RotationalPropagatorType.quaternions") {
         return R"(
Entries 1-4: The quaternion defining the rotation from inertial to body-fixed frame (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#definition-of-rotational-state>`_) Entries 5-7: The body's angular velocity vector, expressed in its body-fixed frame.
     )";


    } else if(name == "RotationalPropagatorType.modified_rodrigues_parameters") {
         return R"(
Entries 1-4: The modified Rodrigues parameters defining the rotation from inertial to body-fixed frame (with entry four the shadow parameter) Entries 5-7: The body's angular velocity vector, expressed in its body-fixed frame.
     )";


    } else if(name == "RotationalPropagatorType.exponential_map") {
         return R"(
Entries 1-4: The exponential map defining the rotation from inertial to body-fixed frame (with entry four the shadow parameter) Entries 5-7: The body's angular velocity vector, expressed in its body-fixed frame.
     )";



    } else if(name == "StateType") {
         return R"(

        Enumeration of available integrated state types.





     )";


    } else if(name == "StateType.hybrid_type") {
         return R"(
     )";


    } else if(name == "StateType.translational_type") {
         return R"(
     )";


    } else if(name == "StateType.rotational_type") {
         return R"(
     )";


    } else if(name == "StateType.body_mass_type") {
         return R"(
     )";


    } else if(name == "StateType.custom_type") {
         return R"(
     )";



    } else if(name == "PropagationTerminationTypes") {
         return R"(

        Enumeration of possible propagation termination types





     )";


    } else if(name == "PropagationTerminationTypes.time_stopping_condition") {
         return R"(
     )";


    } else if(name == "PropagationTerminationTypes.cpu_time_stopping_condition") {
         return R"(
     )";


    } else if(name == "PropagationTerminationTypes.dependent_variable_stopping_condition") {
         return R"(
     )";


    } else if(name == "PropagationTerminationTypes.hybrid_stopping_condition") {
         return R"(
     )";


    } else if(name == "PropagationTerminationTypes.custom_stopping_condition") {
         return R"(
     )";




    } else if(name == "PropagatorSettings") {
         return R"(

        Functional base class to define settings for propagators.

        Base class to define settings for propagators. Derived classes are split into settings for single- and multi-arc dynamics.





     )";




    } else if(name == "PropagatorSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user because this is a base class.





    )";



    } else if(name == "PropagatorSettings.reset_initial_states" && variant==0) {
            return R"(

        Function to reset the initial state used as input for numerical integration.


        Parameters
        ----------
        initial_states : numpy.ndarray
            Initial states to be reset for the numerical propagation.




    )";




    } else if(name == "MultiArcPropagatorSettings") {
         return R"(

        `PropagatorSettings`-derived class to define settings for multi-arc dynamics.





     )";




    } else if(name == "MultiArcPropagatorSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "HybridArcPropagatorSettings") {
         return R"(

        `PropagatorSettings`-derived class to define settings for hybrid-arc dynamics.





     )";




    } else if(name == "HybridArcPropagatorSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "SingleArcPropagatorSettings") {
         return R"(

        `PropagatorSettings`-derived class to define settings for single-arc dynamics.


        Attributes
        ----------
        termination_settings : PropagationTerminationSettings
            Settings for creating the object that checks whether the propagation is finished.




     )";




    } else if(name == "SingleArcPropagatorSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "TranslationalStatePropagatorSettings") {
         return R"(

        `SingleArcPropagatorSettings`-derived class to define settings for single-arc translational dynamics.


        Attributes
        ----------
        acceleration_settings : SelectedAccelerationMap
            Settings for retrieving the accelerations acting on the body during propagation.




     )";




    } else if(name == "TranslationalStatePropagatorSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";



    } else if(name == "TranslationalStatePropagatorSettings.reset_initial_states" && variant==0) {
            return R"(

        Function to reset the initial state used as input for numerical integration.


        Parameters
        ----------
        initial_states : numpy.ndarray
            Initial states to be reset for the numerical propagation.




    )";



    } else if(name == "TranslationalStatePropagatorSettings.recreate_state_derivative_models" && variant==0) {
            return R"(

        Function to (re)create the integrated state models (e.g. acceleration/torque/mass models).

        Function to create the integrated state models (e.g. acceleration/torque/mass models) for
        each fo the propagators state types contained in `propagatorSettingsMap_`.


        Parameters
        ----------
        bodies : SystemOfBodies
            System of bodies used in the propagation.




    )";



    } else if(name == "TranslationalStatePropagatorSettings.single_type_settings" && variant==0) {
            return R"(

        Function to retrieve a single type of propagator.

        Function to retrieve a single type of propagator (translational, rotational or mass). This function is
        often used in multi-type propagation.


        Parameters
        ----------
        state_type : IntegratedStateType
            State type to be retrieved.




    )";




    } else if(name == "RotationalStatePropagatorSettings") {
         return R"(

        `SingleArcPropagatorSettings`-derived class to define settings for single-arc rotational state propagation.





     )";




    } else if(name == "RotationalStatePropagatorSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "MultiTypePropagatorSettings") {
         return R"(

        `SingleArcPropagatorSettings`-derived class to define settings for propagation of multiple quantities.





     )";


    } else if(name == "MultiTypePropagatorSettings.propagator_settings_per_type") {
         return R"(

        None

        :type: dict[IntegratedStateType, list[SingleArcPropagatorSettings]]
     )";




    } else if(name == "MultiTypePropagatorSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";



    } else if(name == "MultiTypePropagatorSettings.reset_initial_states" && variant==0) {
            return R"(

        Function to reset the initial states used as input for numerical integration.


        Parameters
        ----------
        initial_states : numpy.ndarray
            Initial states to be reset for the numerical propagation.




    )";



    } else if(name == "MultiTypePropagatorSettings.recreate_state_derivative_models" && variant==0) {
            return R"(

        Function to (re)create the integrated state models (e.g. acceleration/torque/mass models).

        Function to create the integrated state models (e.g. acceleration/torque/mass models) for
        each of the propagators state types contained in `propagatorSettingsMap_`.


        Parameters
        ----------
        bodies : SystemOfBodies
            System of bodies used in the propagation.




    )";



    } else if(name == "MultiTypePropagatorSettings.single_type_settings" && variant==0) {
            return R"(

        Function to retrieve a single type of propagator.

        Function to retrieve a single type of propagator (translational, rotational or mass). This function is
        often used in multi-type propagation.


        Parameters
        ----------
        state_type : IntegratedStateType
            State type to be retrieved.




    )";




    } else if(name == "PropagationTerminationSettings") {
         return R"(

        Functional base class to define termination settings for the propagation.





     )";




    } else if(name == "PropagationTerminationSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user because this is a base class.





    )";




    } else if(name == "PropagationDependentVariableTerminationSettings") {
         return R"(

        `PropagationTerminationSettings`-derived class to define termination settings for the propagation from dependent variables.





     )";




    } else if(name == "PropagationDependentVariableTerminationSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "PropagationTimeTerminationSettings") {
         return R"(

        `PropagationTerminationSettings`-derived class to define termination settings for the propagation from propagation time.





     )";




    } else if(name == "PropagationTimeTerminationSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "PropagationCPUTimeTerminationSettings") {
         return R"(

        `PropagationTerminationSettings`-derived class to define termination settings for the propagation from CPU time.





     )";




    } else if(name == "PropagationCPUTimeTerminationSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "PropagationCustomTerminationSettings") {
         return R"(

        `PropagationTerminationSettings`-derived class to define custom termination settings for the propagation.





     )";




    } else if(name == "PropagationCustomTerminationSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "PropagationHybridTerminationSettings") {
         return R"(

        `PropagationTerminationSettings`-derived class to define hybrid termination settings for the propagation.





     )";




    } else if(name == "PropagationHybridTerminationSettings.ctor" && variant==0) {
            return R"(

        Constructor.

        Instances of this class are typically not generated by the user. Settings objects for integrators should be
        instantiated through the factory functions of a derived class.





    )";




    } else if(name == "PropagationPrintSettings") {
         return R"(

        Class to save settings on what is to be written to the console during the propagation of a single arc.





     )";


    } else if(name == "PropagationPrintSettings.print_number_of_function_evaluations") {
         return R"(

        Boolean defining whether the number of function evaluations that
        were performed is to be printed to the console (after propagation).


        :type: bool
     )";


    } else if(name == "PropagationPrintSettings.print_propagation_clock_time") {
         return R"(

        Boolean defining whether the total clock time taken for the propagation
        is to be printed to the console (after propagation).


        :type: bool
     )";


    } else if(name == "PropagationPrintSettings.print_termination_reason") {
         return R"(

        Boolean defining whether the reason for propagation termination
        is to be printed to the console (after propagation).


        :type: bool
     )";


    } else if(name == "PropagationPrintSettings.print_initial_and_final_conditions") {
         return R"(

        Boolean defining whether the initial and final conditions (state and time)
        are to be printed to the console (beforee and after propagation, respectively).


        :type: bool
     )";


    } else if(name == "PropagationPrintSettings.results_print_frequency_in_seconds") {
         return R"(

        Variable indicating how often (in seconds of simulation time)
        the current state and time are to be printed to the console (by default, set to NaN - they are never printed).
        In case this setting is active (e.g. not NaN), and the ``results_print_frequency_in_steps`` setting is active,
        the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


        :type: Float
     )";


    } else if(name == "PropagationPrintSettings.results_print_frequency_in_steps") {
         return R"(

        Variable indicating how often (in number of full integration steps)
        the current state and time are to be printed to the console (by default, set to 0 - they are never printed).
        In case this setting is active (e.g. not 0), and the ``results_print_frequency_in_seconds`` setting is active,
        the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


        :type: int
     )";


    } else if(name == "PropagationPrintSettings.print_dependent_variables_during_propagation") {
         return R"(

        Boolean defining whether the dependent variables are to be printed during the propagation along with the state,
        at steps/epochs define by the ``results_print_frequency_in_seconds`` and/or ``results_print_frequency_in_steps`` inputs.


        :type: float
     )";


    } else if(name == "PropagationPrintSettings.print_state_indices") {
         return R"(

        Boolean defining whether the meaning and indices of the
        entries of the state vector are to be printed to
        the console (before the propagation).

        .. note:: The same information can be retrieved from the
                  :py:attr:`SingleArcSimulationResults.propagated_state_ids`
                  attribute.


        :type: bool
     )";


    } else if(name == "PropagationPrintSettings.print_processed_state_indices") {
         return R"(

        Boolean defining whether the meaning and indices of the
        entries of the processed state vector are to be printed to
        the console (after the propagation). The distinction between the 
        propagated and processed (or conventuional) state representation is described in
        detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/processed_propagated_elements.html>`_.
        Summarizing: the processed state is the 'typical' formulation of the state (for translational dynamics: Cartesian states).

        .. note:: The same information can be retrieved from the
                  :py:attr:`SingleArcSimulationResults.processed_state_ids`
                  attribute.
                  


        :type: bool
     )";


    } else if(name == "PropagationPrintSettings.print_dependent_variable_indices") {
         return R"(

        Boolean defining whether the meaning and indices of the
        entries of the dependent variable data are to be printed to
        the console (before the propagation).

        .. note:: The same information can be retrieved from the
                  :py:attr:`SingleArcSimulationResults.dependent_variable_ids`
                  attribute.


        :type: bool
     )";




    } else if(name == "PropagationPrintSettings.enable_all_boolean_printing" && variant==0) {
            return R"(

        Function enabling all True/False printing (e.g. sets all boolean attributes to True)       






    )";



    } else if(name == "PropagationPrintSettings.enable_all_printing" && variant==0) {
            return R"(

        Function enabling all True/False printing (e.g. sets all boolean attributes to True), and setting the non-boolean 
        attributes to values defined here. 
              



        Parameters
        ----------
        results_print_frequency_in_seconds : float
            See ``results_print_frequency_in_seconds`` class attribute
        results_print_frequency_in_steps : int
            See ``results_print_frequency_in_steps`` class attribute




    )";



    } else if(name == "PropagationPrintSettings.disable_all_printing" && variant==0) {
            return R"(

        Function enabling all printing (e.g. sets all boolean attributes to False, and disables all other output as well)






    )";




    } else if(name == "PropagatorProcessingSettings") {
         return R"(

        Base class to define settings on how the numerical results are to be used

        Base class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
        Instances of this class are typically not created by the user. Settings objects for derived class of single-, multi- and hybrid arc propagation are 
        instantiated through the factory functions to define propagator settings (such as :func:`~translational` or :func:`~multi_arc`) in this module
            





     )";


    } else if(name == "PropagatorProcessingSettings.clear_numerical_solution") {
         return R"(

        Boolean defining whether the propagation results should be
        deleted after the propagation is terminated. If this is
        done, the :py:attr:`~state_history`,
        :py:attr:`~unprocessed_state_history` and
        :py:attr:`~dependent_variable_history` will not be
        available in the :py:class:`~tudatpy.numerical_simulation.propagator.SingleArcSimulationResults` class. Putting this setting to True (deleting the
        results) is only sensible when the
        :py:attr:`~set_integrated_result` is set to True. In that
        case, the propagated states are *not* accessible directly
        from this objects, but the results are used to update the
        environment, *e.g.* update the ephemeris of the propagated
        body with the numerical results.


        :type: bool
     )";


    } else if(name == "PropagatorProcessingSettings.set_integrated_result") {
         return R"(

        Boolean defining whether the propagation results are to
        be used to update the environment. If this variable is set
        to False, the numerical propagation results can be
        retrieved from this object (provided the
        :py:attr:`~clear_numerical_solution` is set to False),
        but the (for instance) Ephemeris of the propagated body
        is not updated with the propagation results. If this
        variable is set to True, the properties of the propagated
        :class:`~tudatpy.numerical_simulation.environment.Body`
        object will be updated as per the numerical results
        (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/printing_processing_results.html#automatic-processing>`_ for details).


        :type: bool
     )";





    } else if(name == "SingleArcPropagatorProcessingSettings") {
         return R"(

        Class to define settings on how the numerical results are to be used for single-arc propagations 

        Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
        derived from :class:`PropagatorProcessingSettings`.
        Instances of this class are typically not created by the user. A settings object is
        instantiated through the factory functions to define single-arc propagator settings (such as :func:`~translational` or :func:`~rotational`) in this module





     )";


    } else if(name == "SingleArcPropagatorProcessingSettings.print_settings") {
         return R"(

        Settings object defining which quantities should be printed to the console before, during and after the propagation. By default, this
        object is instantiated to print nothing.


        :type: PropagationPrintSettings
     )";


    } else if(name == "SingleArcPropagatorProcessingSettings.results_save_frequency_in_steps") {
         return R"(

        Variable indicating how often (in number of integrator steps)
        the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results 
        (by default, set to 1 - they are saved every time step). If this setting is set to 0, the data is never saved based on number of steps.
        In case this setting is active (e.g. not 0), and the ``results_save_frequency_in_seconds`` setting is active,
        the data is saved as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


        :type: int
     )";


    } else if(name == "SingleArcPropagatorProcessingSettings.results_save_frequency_in_seconds") {
         return R"(

        Variable indicating how often (in seconds of simulation time)
        the propagated time, state, dependent variables, etc. are to be saved to data structures containing the results 
        (by default, set to NaN - they are not saved based on propagation time; see below and ``results_save_frequency_in_steps`` attribute ).
        In case this setting is active, and set to :math:`\Delta t`, the data are saved as soon as the current time step is :math:`\ge \Delta t` after the
        last step at which data was saved.
        In case this setting is active (e.g. not NaN), and the ``results_save_frequency_in_steps`` setting is active,
        the current state is printed as soon as *one* of the two conditions (number of seconds, or number of steps) is met.


        :type: float
     )";





    } else if(name == "MultiArcPropagatorProcessingSettings") {
         return R"(

        Class to define settings on how the numerical results are to be used for multi-arc propagations 

        Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
        derived from :class:`PropagatorProcessingSettings`.
        Instances of this class are typically not created by the user. A settings object is
        instantiated through the factory function :func:`~multi_arc` to define multi-arc propagator settings.
        This object contains a list of :class:`SingleArcPropagatorProcessingSettings` objects, containing the processing settings for each constituent arc.





     )";


    } else if(name == "MultiArcPropagatorProcessingSettings.print_output_on_first_arc_only") {
         return R"(

        **read-only**

        Variable defining whether the ``set_print_settings_for_all_arcs`` function has been used to define identical print settings for each arc.


        :type: bool
     )";


    } else if(name == "MultiArcPropagatorProcessingSettings.print_output_on_first_arc_only") {
         return R"(

        **read-only**

        Variable defining whether the ``set_print_settings_for_all_arcs`` function has been used to define identical print settings for each arc.


        :type: bool
     )";


    } else if(name == "MultiArcPropagatorProcessingSettings.single_arc_settings") {
         return R"(

        **read-only**

        List containing the processing settings for each constituent arc


        :type: list[SingleArcPropagatorProcessingSettings]
     )";




    } else if(name == "MultiArcPropagatorProcessingSettings.set_print_settings_for_all_arcs" && variant==0) {
            return R"(

        Function that sets the same print settings for each arc in the multi-arc propagation.



        Parameters
        ----------
        single_arc_print_settings : PropagationPrintSettings
            Propagation print settings that are applied to each constituent single-arc settings, overriding any existing settings.




    )";




    } else if(name == "HybridArcPropagatorProcessingSettings") {
         return R"(

        Class to define settings on how the numerical results are to be used for hybrid-arc propagations 

        Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment),
        derived from :class:`PropagatorProcessingSettings`.
        Instances of this class are typically not created by the user. A settings object is
        instantiated through the factory function :func:`~hybrid_arc` to define hybrid-arc propagator settings.
        This object contains a :class:`SingleArcPropagatorProcessingSettings` object and a :class:`MultuArcPropagatorProcessingSettings` , 
        containing the processing settings for the constituents of the hybrid-arc propagatioon.





     )";


    } else if(name == "HybridArcPropagatorProcessingSettings.single_arc_settings") {
         return R"(

        Processing settings for the single-arc component of the hybrid-arc propagation.


        :type: SingleArcPropagatorProcessingSettings
     )";


    } else if(name == "HybridArcPropagatorProcessingSettings.multi_arc_settings") {
         return R"(

        Processing settings for the single-arc component of the multi-arc propagation.


        :type: MultiArcPropagatorProcessingSettings
     )";




    } else if(name == "HybridArcPropagatorProcessingSettings.set_print_settings_for_all_arcs" && variant==0) {
            return R"(

        Function that sets the same print settings for each arc in the multi-arc propagation, and the single-arc propagation.



        Parameters
        ----------
        single_arc_print_settings : PropagationPrintSettings
            Propagation print settings that are applied to each arc in the multi-arc propagation, and the single-arc propagation, overriding any existing settings.




    )";





    } else if(name == "translational" && variant==0) {
        return R"(
        
Factory function to create translational state propagator settings with stopping condition at given final time.

Factory function to create translational state propagator settings for N bodies.
The propagated state vector is defined by the combination of integrated bodies, and their central body, the combination
of which define the relative translational states for which a differential equation is to be solved. The propagator
input defines the formulation in which the differential equations are set up
The dynamical models are defined by an ``AccelerationMap``, as created by :func:`~create_acceleration_models` function.
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational.html>`_


Parameters
----------
central_bodies : list[str]
    List of central bodies with respect to which the bodies to be integrated are propagated.
acceleration_models : AccelerationMap
    Set of accelerations acting on the bodies to propagate, provided as acceleration models.
bodies_to_integrate : list[str]
    List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
initial_states : numpy.ndarray
    Initial states of the bodies to integrate (one initial state for each body, concatenated into a single array), provided in the same order as the bodies to integrate. The initial states must be expressed in Cartesian elements, w.r.t. the central body of each integrated body. The states must be defined with the same frame orientation as the global frame orientation of the environment (specified when creating a system of bodies, see for instance :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings` and :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`). Consequently, for N integrated bodies, this input is a vector with size size 6N.
initial_time : float
    Initial epoch of the numerical propagation
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
propagator : TranslationalPropagatorType, default=cowell
    Type of translational propagator to be used (see `TranslationalPropagatorType` enum).
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
TranslationalStatePropagatorSettings
    Translational state propagator settings object.






    )";



    } else if(name == "rotational" && variant==0) {
        return R"(
        
Factory function to create rotational state propagator settings.

Factory function to create rotational state propagator settings for N bodies.
The propagated state vector is defined by the integrated bodies, which defines the bodies for which the
differential equation defining the evolution of the rotational state between an
inertial and body-fixed frame are to be solved. The propagator input defines the
formulation in which the differential equations are set up. The dynamical models are
defined by an ``TorqueModelMap``, as created by ``create_torque_models`` function.
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/rotational.html>`_


Parameters
----------
torque_models : TorqueModelMap
    Set of torques acting on the bodies to propagate, provided as torque models.
bodies_to_integrate : list[str]
    List of bodies to be numerically propagated, whose order reflects the order of the central bodies.
initial_states : numpy.ndarray
    Initial rotational states of the bodies to integrate (one initial state for each body), provided in the same order as the bodies to integrate.
    Regardless of the propagator that is selected, the initial rotational state is always defined as four quaternion entries, and the angular velocity of the body,
    as defined in more detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html#definition-of-rotational-state>`_.

initial_time : float
    Initial epoch of the numerical propagation
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
propagator : RotationalPropagatorType, default=quaternions
    Type of rotational propagator to be used (see `RotationalPropagatorType` enum).
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
RotationalStatePropagatorSettings
    Rotational state propagator settings object.






    )";



    } else if(name == "mass" && variant==0) {
        return R"(
        
Factory function to create mass propagator settings

Factory function to create mass propagator settings 
It works by providing a key-value mass rate container, containing the list of mass rate settings objects associated to
each body. In this function, the dependent variables to save are provided
as a list of SingleDependentVariableSaveSettings objects. In this function, the termination conditions are set
through the termination settings object provided.
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/mass.html>`_


Parameters
----------
bodies_with_mass_to_propagate : list[str]
    List of bodies whose mass should be numerically propagated.
mass_rate_settings : SelectedMassRateModelMap
    Mass rates associated to each body, provided as a mass rate settings object.
initial_body_masses : numpy.ndarray
    Initial masses of the bodies to integrate (one initial mass for each body), provided in the same order as the bodies to integrate.
initial_time : float
    Initial epoch of the numerical propagation
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
MassPropagatorSettings
    Mass propagator settings object.






    )";



    } else if(name == "multitype" && variant==0) {
        return R"(
        
Factory function to create multitype propagator settings.

Factory function to create multitype propagator settings.
It works by providing a list of SingleArcPropagatorSettings objects. When using this function,
only the termination and output settings provided here are used, any such settings in the
constituent propagator settings are ignored
Details on the usage of this function are discussed in more detail in the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/multi_type.html>`_

.. note:: The propagated state contains the state types in the following order: Translational ( **C** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
          When propagating two bodies, an example of what the output state would look like is for instance:
          [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ]


Parameters
----------
propagator_settings_list : list[SingleArcPropagatorSettings]
    List of SingleArcPropagatorSettings objects to use.
integrator_settings : IntegratorSettings
    Settings defining the numerical integrator that is to be used for the propagation

    .. note:: The sign of the initial time step in the integrator settings defines whether the propagation will be forward or backward in time

initial_time : float
    Initial epoch of the numerical propagation
termination_settings : PropagationTerminationSettings
    Generic termination settings object to check whether the propagation should be ended.
output_variables : list[SingleDependentVariableSaveSettings], default=[]
    Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
MassPropagatorSettings
    Mass propagator settings object.






    )";



    } else if(name == "multi_arc" && variant==0) {
        return R"(
        
Factory function to create multi-arc propagator settings.

Factory function to create multi-arc propagator settings. It works by providing separate settings for
each arc in a list.


Parameters
----------
single_arc_settings : list[SingleArcPropagatorSettings]
    List of SingleArcPropagatorSettings objects to use, one for each arc.
transfer_state_to_next_arc : bool, default=False
    Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
MultiArcPropagatorSettings
    Multi-arc propagator settings object.






    )";



    } else if(name == "hybrid_arc" && variant==0) {
        return R"(
        
Factory function to create hybrid-arc propagator settings.

Factory function to create hybrid-arc propagator settings (i.e., a combination of single- and multi-arc dynamics).


Parameters
----------
single_arc_settings : SingleArcPropagatorSettings
    SingleArcPropagatorSettings object to use for the propagation.
multi_arc_settings : MultiArcPropagatorSettings
    Class to define settings on how the numerical results are to be used, both during the propagation (printing to console) and after propagation (resetting environment)
Returns
-------
HybridArcPropagatorSettings
    Hybrid-arc propagator settings object.






    )";



    } else if(name == "time_termination" && variant==0) {
        return R"(
        
Factory function to create time termination settings for the propagation.

Factory function to create time termination settings for the propagation.
The propagation is stopped when the final time provided is reached. Note that the termination time is set as the
absolute time (in seconds since J2000), not the time since the start of the propagation.
Depending on the sign of the time step of the numerical integrator, the termination time will be treated as an
upper bound (for positive time step) or lower bound (for negative time step).
The simulator will normally finish the final time-step, which may cause the termination time to be slightly exceeded.
This behaviour can be suppressed by providing the optional input argument
``terminate_exactly_on_final_condition=True``, in which case the final propagation step will be *exactly* on the
specified time.


Parameters
----------
termination_time : float
    Final time of the propagation.
terminate_exactly_on_final_condition : bool, default=False
    Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
Returns
-------
PropagationTimeTerminationSettings
    Time termination settings object.



Notes
-----
To reach *exactly* the final time, state derivative function evaluations beyond the final
time may be required by the propagator. Reaching the final condition exactly is an iterative process and
very minor deviations from the specified final condition can occur.



Examples
--------
In this example, we set the termination time of the propagation equal to one day (86400 $s$).

.. code-block:: python

  # Set termination time (in seconds since J2000)
  termination_time = simulation_start_epoch + 86400.0
  # Create time termination settings
  termination_settings = propagation_setup.propagator.time_termination( termination_time )


    )";



    } else if(name == "cpu_time_termination" && variant==0) {
        return R"(
        
Factory function to create CPU time termination settings for the propagation.

Factory function to create CPU time termination settings for the propagation.
The propagation is stopped when the final CPU time provided is reached.


Parameters
----------
cpu_termination_time : float
    Maximum CPU time for the propagation.
Returns
-------
PropagationCPUTimeTerminationSettings
    CPU time termination settings object.





Examples
--------
In this case, we set a CPU time termination setting so that the propagation stops once your computer has run it
for 120 seconds.

.. code-block:: python

  # Set CPU time to 120 seconds
  cpu_termination_time = 120.0
  # Create termination settings
  termination_settings = propagation_setup.propagator.cpu_time_termination( cpu_termination_time )


    )";



    } else if(name == "dependent_variable_termination" && variant==0) {
        return R"(
        
Factory function to create termination settings for the propagation based on a dependent variable.

Factory function to create termination settings for the propagation based on the value of a dependent variable.
The propagation is stopped when a provided upper or lower limit value is reached.
The simulator will normally finish the final time-step, which may cause the dependent variable to be slightly exceeded.
This behaviour can be suppressed by providing the optional input argument
``terminate_exactly_on_final_condition=True``, in which case the final propagation step will be *exactly* on the
specified dependent variable value.


Parameters
----------
dependent_variable_settings : SingleDependentVariableSaveSettings
    Dependent variable object to be used as termination setting.
limit_value : float
    Limit value of the dependent variable; if reached, the propagation is stopped.
use_as_lower_limit : bool, default=False
    Denotes whether the limit value should be used as lower or upper limit.
terminate_exactly_on_final_condition : bool, default=False
    Denotes whether the propagation is to terminate exactly on the final condition, or whether it is to terminate on the first step where it is violated.
termination_root_finder_settings : bool, default=None
    Settings object to create root finder used to converge on exact final condition.
Returns
-------
PropagationDependentVariableTerminationSettings
    Dependent variable termination settings object.



Notes
-----
To reach *exactly* the final dependent variable value, state derivative function evaluations beyond the final
time may be required by the propagator. Reaching the final condition exactly is an iterative process and
very minor deviations from the specified final condition can occur.



Examples
--------
Below, an example is shown for termination on a given vehicle altitude. The exact termination condition is defined
in the ``termination_settings``. The propagation is terminated once the *lower* limit of 25 km in altitude is
reached (as the ``use_as_lower_limit`` is set to ``True``). To use the above settings to terminate when an
*upper* limit of 25 km is reached, set this boolean to ``False``. In this example, we also want to stop exactly
at 25 km, so we set ``terminate_exactly_on_final_condition`` to ``True``, and we specify ``termination_root_finder_settings``.

.. code-block:: python

  # Set dependent variable to be checked as termination setting
  termination_variable = propagation_setup.dependent_variable.altitude( "Spacecraft", "Earth" )
  # Create termination settings
  termination_settings = propagation_setup.propagator.dependent_variable_termination(
    dependent_variable_settings = termination_variable,
    limit_value = 25.0E3,
    use_as_lower_limit = True,
    terminate_exactly_on_final_condition=True,
    termination_root_finder_settings=root_finders.secant(
        maximum_iteration=5,
        maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result)
    )
  )


    )";



    } else if(name == "custom_termination" && variant==0) {
        return R"(
        
Factory function to create custom termination settings for the propagation.

Factory function to create custom termination settings for the propagation.
The propagation is stopped when the condition provided is verified.
This custom function should take the current time as input and output a Boolean. It can use internal variables
and calculations, for example retrieved from the environment.


Parameters
----------
custom_condition : callable[[float], bool]
    Function of time (independent variable) which is called during the propagation and returns a boolean value denoting whether the termination condition is verified.
Returns
-------
PropagationCustomTerminationSettings
    Custom termination settings object.





Examples
--------

.. code-block:: python

  # Create custom function returning a bool
  def custom_termination_function(time: float):
      # Do something
      set_condition = ...
      # Return bool
      return set_condition

  # Create termination settings
  termination_settings = propagation_setup.propagator.custom_termination(
    custom_termination_function)


    )";



    } else if(name == "hybrid_termination" && variant==0) {
        return R"(
        
Factory function to create hybrid termination settings for the propagation.

Factory function to create hybrid termination settings for the propagation. This function can be used
to define that all conditions or a single condition of the conditions provided must be met to
stop the propagation. Each termination condition should be created according to each individual factory function
and then added to a list of termination conditions.

Note that, when using this option, the :attr:`~tudatpy.numerical_simulation.propagation.SingleArcSimulationResults.termination_details` of
the simulation results object (obtained from here after a propagation: `:attr:`~tudatpy.numerical_simulation.SingleArcSimulator.propagation_results`)
is of derived type :class:`~tudatpy.numerical_simulation.propagation.PropagationTerminationDetailsFromHybridCondition`, which
contains additional details on the hybrid termination (such as the specific conditions that were met).


Parameters
----------
termination_settings : list[PropagationTerminationSettings]
    List of single PropagationTerminationSettings objects to be checked during the propagation.
fulfill_single_condition : bool, default=False
    Whether only a single condition of those provided must be met to stop the propagation (true) or all of them simultaneously (false).
Returns
-------
PropagationHybridTerminationSettings
    Hybrid termination settings object.





Examples
--------
In the following example, the propagation will terminate once *one of the three* termination settings (simulated time, cpu time, altitude)
has reached the imposed limit value. The ``fulfill_single_condition`` variable determines whether the propagation
terminates once a *single* condition is met (if True, as above) or once *all* conditions must be met (False).

.. code-block:: python

  # Set simulation termination time
  termination_time = simulation_start_epoch + 86400.0
  # Create simulation time termination setting
  time_termination_settings = propagation_setup.propagator.time_termination( termination_time )

  # Set dependent variable termination setting
  termination_variable = propagation_setup.dependent_variable.altitude( "Spacecraft", "Earth" )
  # Create altitude-based termination setting
  altitude_termination_settings = propagation_setup.propagator.dependent_variable_termination(
    dependent_variable_settings = termination_variable,
    limit_value = 25.0E3,
    use_as_lower_limit = True)

  # Set cpu termination time
  cpu_termination_time = 120.0
  # Create cpu time termination setting
  cpu_termination_settings = propagation_setup.propagator.cpu_time_termination( cpu_termination_time )

  # Store termination setting objects in a list
  termination_settings_list = [time_termination_settings, altitude_termination_settings, cpu_termination_settings]

  # Create hybrid termination settings
  termination_settings = propagation_setup.propagator.hybrid_termination( termination_settings_list, fulfill_single_condition = True )


    )";



    } else if(name == "non_sequential_termination" && variant==0) {
        return R"(
        
Factory function to create non-sequential termination settings for the propagation.

Factory function to create non-sequential termination settings for the propagation. By using this setting,
the propagation of the dynamics along an arc is propagated starting from some point (initial time and state) along the arc, and then
propagating both forwards and backwards in time. This termination condition allows the user to specify 
termination conditions for the propagations forwards and backwards in time. These two propagations are then
internally performed separately, but the resulting propagation results provide the concatenated results
from the two and effectively constitute
By using this function, the propagation controlled by the ``forward_termination_settings`` automatically has a positive
time step, while the ``backward_termination_settings`` automatically has a negative time step. By definition,
both will start from the same initial time and state, which are provided in the propagation settings


Parameters
----------
forward_termination_settings : PropagationTerminationSettings
    Propagation termination setting for the forward-in-time propagation
backward_termination_settings : PropagationTerminationSettings
    Propagation termination setting for the backwards-in-time propagation
Returns
-------
PropagationTerminationSettings
    Termination settings object for forward- and backwards-in-time propagation.






    )";



    } else if(name == "add_dependent_variable_settings" && variant==0) {
        return R"(
        
Function to add dependent variables to existing propagator settings.

Function to add dependent variables to existing :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.SingleArcPropagatorSettings`
object. This function is added as an alternative to teh regular manner in which to defined dependent variables (use of input to factory
functions for single-arc propagator settings :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.translational`,
:func:`~tudatpy.numerical_simulation.propagation_setup.propagator.rotational`, :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.mass`,
:func:`~tudatpy.numerical_simulation.propagation_setup.propagator.multitype`). Typically, this function is used to modify
existing propagator settings in a loop when running multiple simulations


Parameters
----------
dependent_variable_settings : List[ SingleDependentVariableSaveSettings ]
    List of dependent variable settings that are to be added to propagator settings. Note that this function adds settings, and does not replace any existing settings (nor does it check for duplicate settings).
propagator_settings : SingleArcPropagatorSettings
    Propagator settings to which the additional dependent variables settings are to be added.
Returns
-------
None
    None






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace torque {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "AvailableTorque") {
         return R"(

        Enumeration of available torque types.

        Enumeration of torque types supported by tudat.





     )";


    } else if(name == "AvailableTorque.torque_free_type") {
         return R"(
     )";


    } else if(name == "AvailableTorque.undefined_torque_type") {
         return R"(
     )";


    } else if(name == "AvailableTorque.second_order_gravitational_torque_type") {
         return R"(
     )";


    } else if(name == "AvailableTorque.aerodynamic_torque_type") {
         return R"(
     )";


    } else if(name == "AvailableTorque.spherical_harmonic_gravitational_torque_type") {
         return R"(
     )";


    } else if(name == "AvailableTorque.inertial_torque_type") {
         return R"(
     )";


    } else if(name == "AvailableTorque.dissipative_torque_type") {
         return R"(
     )";


    } else if(name == "AvailableTorque.custom_torque_type") {
         return R"(
     )";




    } else if(name == "TorqueSettings") {
         return R"(

        Functional base class to define settings for torques.

        This is a functional base class to define settings for torques that require no information in addition to their type.
        Classes defining settings for torque models requiring additional information must be
        derived from this class.
        Bodies exerting and undergoing torque are set outside of this class.
        This class can be used for the easy setup of torque models
        (see createTorqueModels.h), but users may also chose to do so manually.
        (Derived) Class members are all public, for ease of access and modification.





     )";





    } else if(name == "SphericalHarmonicTorqueSettings") {
         return R"(

        `TorqueSettings`-derived class to define settings for torques caused by spherical harmonic gravity.

        `TorqueSettings`-derived class to define settings for torques caused by spherical harmonic gravity.





     )";






    } else if(name == "aerodynamic" && variant==0) {
        return R"(
        
Creates the settings for the aerodynamic torque.

Creates the settings for the aerodynamic torque exerted by a body with an atmosphere model and shape model on
another body. The body exerting the torque needs to have both an atmosphere model and a shape model defined.
Furthermore, the body undergoing the torque needs to have the aerodynamic coefficient interface and its moment
coefficients defined. In the case that the aerodynamic coefficients are defined as a function of the vehicle
orientation (e.g. angle of attack and sideslip angle), these angles can be manually or automatically defined.

Returns
-------
TorqueSettings
    Torque settings object.





Examples
--------

In this example, we define the aerodynamic torque exerted by the Earth on the vehicle.

.. code-block:: python

  # Create torque settings dict
  torque_settings_vehicle = {}
  # Add aerodynamic torque exerted by the Earth on the vehicle
  torque_settings_vehicle["Earth"] = [propagation_setup.torque.aerodynamic()]


    )";



    } else if(name == "spherical_harmonic_gravitational" && variant==0) {
        return R"(
        
Creates the settings for the spherical harmonic torque.

Torque exerted by a point mass on a body with an arbitrary degree/order spherical harmonics mass distribution.
The body exerting the torque only needs to have a gravitational model defined (point-mass or spherical harmonic),
while the body undergoing the torque needs to have a spherical harmonic gravity field defined.


Parameters
----------
maximum_degree : int
    Maximum degree of the spherical harmonic expansion.
maximum_order : int
    Maximum order of the spherical harmonic expansion.
Returns
-------
TorqueSettings
    Torque settings object.





Examples
--------

In this example, we define the spherical harmonic gravitational torque (up to degree 4 and order 4)
exerted by the Earth on the vehicle.

.. code-block:: python

  # Create torque settings dict
  torque_settings_vehicle = {}
  # Add aerodynamic torque exerted by the Earth on the vehicle
  torque_settings_vehicle["Earth"] = [propagation_setup.torque.spherical_harmonic_gravitational(4, 4)]


    )";



    } else if(name == "second_degree_gravitational" && variant==0) {
        return R"(
        
Creates the settings for the second-degree gravitational torque.

Torque exerted by a point mass on a body with a degree two spherical harmonics mass distribution.
A degree two spherical harmonics mass distribution can be represented by an inertia tensor; thus,
for this torque model, the body undergoing the torque needs to have an inertia tensor defined.
The body exerting the torque only needs to have a gravitational model defined (either point-mass or spherical
harmonics).

Returns
-------
TorqueSettings
    Torque settings object.





Examples
--------

In this example, we define the second degree gravitational torque
exerted by the Earth on the vehicle.

.. code-block:: python

  # Create torque settings dict
  torque_settings_vehicle = {}
  # Add aerodynamic torque exerted by the Earth on the vehicle
  torque_settings_vehicle["Earth"] = [propagation_setup.torque.second_degree_gravitational()]


    )";



    } else if(name == "custom_torque" && variant==0) {
        return R"(
        
Creates the settings for a custom torque.

Creates settings for a custom torque. This torque must be parameterized as a function of time
and expressed with an inertial orientation.


Parameters
----------
torque_function : callable[[float], list]
    Custom torque function with time as an independent variable.
scaling_function : callable[[float], float], default=None
    Scaling function with time as an independent variable to be multiplied by the custom torque function.
Returns
-------
TorqueSettings
    Torque settings object.





Examples
--------

In this example, we define a custom torque
exerted by the Earth on the vehicle.

.. code-block:: python

  # Create torque function
  def torque_function(time: float):
      # Compute torque
      torque = ...
      return torque

  # Create torque settings dict
  torque_settings_vehicle = {}
  # Add aerodynamic torque exerted by the Earth on the vehicle
  torque_settings_vehicle["Earth"] = [propagation_setup.torque.custom(torque_function)]


    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace thrust {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "ThrustFrames") {
         return R"(

        Enumeration of available thrust frame types.

        Enumeration of thrust frame types supported by tudat. The inertial frame has axes along the global orientation (as
        defined in the :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` object. The TNW is determined
        w.r.t. a central body as described in the function :func:`~tudatpy.astro.frame_conversion.tnw_to_inertial_rotation_matrix`
        (with the N-axis pointing away from the central body).





     )";


    } else if(name == "ThrustFrames.unspecified_thrust_frame_type") {
         return R"(
     )";


    } else if(name == "ThrustFrames.inertial_thrust_frame_type") {
         return R"(
     )";


    } else if(name == "ThrustFrames.tnw_thrust_frame_type") {
         return R"(
     )";



    } else if(name == "ThrustMagnitudeTypes") {
         return R"(

        Enumeration of available thrust magnitude types.





     )";


    } else if(name == "ThrustMagnitudeTypes.constant_thrust_magnitude") {
         return R"(
     )";


    } else if(name == "ThrustMagnitudeTypes.from_engine_properties_thrust_magnitude") {
         return R"(
     )";


    } else if(name == "ThrustMagnitudeTypes.thrust_magnitude_from_time_function") {
         return R"(
     )";


    } else if(name == "ThrustMagnitudeTypes.thrust_magnitude_from_dependent_variables") {
         return R"(
     )";


    } else if(name == "ThrustMagnitudeTypes.bang_bang_thrust_magnitude_from_mee_costates") {
         return R"(
     )";




    } else if(name == "ThrustMagnitudeSettings") {
         return R"(

        Functional base class to define settings for the thrust magnitude.


        Attributes
        ----------
        thrust_magnitude_type : ThrustMagnitudeType
            Thrust magnitude type object.
        thrust_origin_id : str
            Reference ID of the thrust origin that should be used (empty if N/A).




     )";





    } else if(name == "ConstantThrustMagnitudeSettings") {
         return R"(

        `ThrustMagnitudeSettings`-derived class to define settings for constant thrust magnitude.

        Derived class to provide settings for the thrust magnitude. This class should be used to define a constant thrust
        magnitude.


        Attributes
        ----------
        thrust_magnitude : float
            Value of the constant thrust magnitude.
        specific_impulse : float
            Value of the constant specific impulse.
        specific_impulse : numpy.ndarray
            Thrust direction vector expressed in the body-fixed reference frame.




     )";





    } else if(name == "CustomThrustMagnitudeSettings") {
         return R"(

        `ThrustMagnitudeSettings`-derived class to define settings for constant thrust magnitude.

        Derived class to provide settings for the thrust magnitude. This class should be used to define a thrust
        magnitude through a custom function.





     )";






    } else if(name == "get_propulsion_input_variables" && variant==0) {
        return R"(
        
Function to create a list of functions that compute and return independent variables for the thrust.

Function to create a list of functions that compute and return independent variables for thrust and/or specific
impulse. This parameterization is used to create a specific thrust magnitude type (see thrust magnitude from
dependent variables). This function retrieves all input functions from the environment and a list of user-defined
functions.


Parameters
----------
body_with_guidance : Body
    Body object whose thrust guidance should be defined.
independent_variables : list[ThrustIndependentVariables]
    Set of dependent variables that should be used to compute the thrust.
guidance_input_functions : list[callable[[], float], default=[]
    Set of functions to compute the thrust, each associated to a specific dependent variable.





    )";



    } else if(name == "constant_thrust_magnitude" && variant==0) {
        return R"(
        
Create thrust magnitude settings from a constant thrust magnitude and Isp.

Factory function that creates constant thrust magnitude settings. The specific impulse to use for the thrust is
also supplied when applying a mass rate model in the propagation of the vehicle dynamics, relating the thrust
to the mass decrease of the vehicle.


Parameters
----------
thrust_magnitude : float
    Value of the constant thrust magnitude.
specific_impulse : float
    Value of the constant specific impulse, used to link the thrust model to the mass propagation.
Returns
-------
ConstantThrustMagnitudeSettings
    Constant thrust magnitude settings object.





Examples
--------
In this example, we define constant thrust magnitude of 1.5 kN and a specific impulse of 315 s. 

.. code-block:: python 
  
  # Define constant thrust magnitude settings of 1.5kN, an Isp of 315s 
  thrust.constant_thrust_magnitude( 
      thrust_magnitude=1.5e3, 
      specific_impulse=315 
  ) 


    )";



    } else if(name == "custom_thrust_magnitude" && variant==0) {
        return R"(
        
Create thrust magnitude settings from a custom thrust force magnitude function.

Factory function that creates thrust magnitude from a custom thrust force magnitude function.
This model defines a thrust force and specific impulse that can vary with time. The thrust acceleration
is computed during the propagation by dividing the thrust force by the current vehicle mass.
The specific impulse can be used to apply a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass
decrease of the vehicle.


Parameters
----------
thrust_magnitude_function : callable[[float], float]
    Function of time returning the value of the thrust force magnitude.
specific_impulse_function : callable[[float], float]
    Function of time returning the value of the specific impulse, useful to link the mass propagation to the thrust model.
Returns
-------
FromFunctionThrustMagnitudeSettings
    From function thrust magnitude settings object.





Examples
--------
In this example, we define a thrust force magnitude based on a set of custom functions. 
The magnitude itself starts from 500N, and linearly increases with time. 
The specific impulse is constant, at 350s. Note that we use a `lambda` function to achieve this neatly. 
Finally, the engine is setup to work for 50s, and be turned off afterwards. 

.. code-block:: python 
  
  # Define the thrust magnitude function: thrust increases linearly with time 
  def thrust_magnitude_function(time): 
      return 500 + time/2 

  # Define a lambda specific impulse function: constant at 350s 
  specific_impulse_function = lambda time: 350 

  # Define the custom thrust magnitude settings based on the pre-defined functions 
  thrust.custom_thrust_magnitude(thrust_magnitude_function, specific_impulse_function ) 


    )";



    } else if(name == "custom_thrust_magnitude_fixed_isp" && variant==0) {
        return R"(
        
Same as :func:`~custom_thrust_magnitude`, but with a fixed value for the specific impulse.


Parameters
----------
thrust_magnitude_function : callable[[float], float]
    Function of time returning the value of the thrust force magnitude.
specific_impulse : float
    Constant value for specific impulse, useful to link the mass propagation to the thrust model.
Returns
-------
FromFunctionThrustMagnitudeSettings
    From function thrust magnitude settings object.






    )";



    } else if(name == "custom_thrust_acceleration_magnitude" && variant==0) {
        return R"(
        
Create thrust magnitude settings from a custom thrust acceleration magnitude function.

Factory function that creates thrust magnitude from a custom thrust acceleration magnitude function.
This model is similar to the :func:`~custom_thrust_magnitude`, with the difference being that this function
directly provides the thrust *acceleration*, not the thrust *force*.


Parameters
----------
thrust_acceleration_magnitude_function : callable[[float], float]
    Function of time returning the value of the thrust acceleration magnitude.
specific_impulse_function : callable[[float], float]
    Function of time returning the value of the specific impulse, useful to link the mass propagation to the thrust model.
Returns
-------
FromFunctionThrustMagnitudeSettings
    From function thrust magnitude settings object.






    )";



    } else if(name == "custom_thrust_acceleration_magnitude_fixed_isp" && variant==0) {
        return R"(
        
Same as :func:`~custom_thrust_acceleration_magnitude`, but with a fixed value for the specific impulse.


Parameters
----------
thrust_acceleration_magnitude_function : callable[[float], float]
    Function of time returning the value of the thrust acceleration magnitude.
specific_impulse : float
    Constant value for specific impulse, useful to link the mass propagation to the thrust model.
Returns
-------
FromFunctionThrustMagnitudeSettings
    From function thrust magnitude settings object.






    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace estimation {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "EstimatableParameterSet") {
         return R"(

        Class containing a consolidated set of estimatable parameters.

        Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
        The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.create_parameters_to_estimate` factory function.





     )";


    } else if(name == "EstimatableParameterSet.parameter_set_size") {
         return R"(

        **read-only**

        Size of the parameter set, i.e. amount of estimatable parameters contained in the set.

        :type: int
     )";


    } else if(name == "EstimatableParameterSet.initial_states_size") {
         return R"(

        **read-only**

        Amount of initial state parameters contained in the set.

        :type: int
     )";


    } else if(name == "EstimatableParameterSet.initial_single_arc_states_size") {
         return R"(

        **read-only**

        Amount of initial state parameters in the set, which are treated in a single-arc fashion.

        :type: int
     )";


    } else if(name == "EstimatableParameterSet.initial_multi_arc_states_size") {
         return R"(

        **read-only**

        Amount of initial state parameters in the set, which are treated in a multi-arc fashion.

        :type: int
     )";


    } else if(name == "EstimatableParameterSet.constraints_size") {
         return R"(

        **read-only**

        Total size of linear constraint that is to be applied during estimation.

        :type: int
     )";


    } else if(name == "EstimatableParameterSet.parameter_vector") {
         return R"(

        Vector containing the parameter values of all parameters in the set.

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )";




    } else if(name == "EstimatableParameterSet.indices_for_parameter_type" && variant==0) {
            return R"(

        Function to retrieve the indices of a given type of parameter.

        Function to retrieve the index of all parameters of a given type from the parameter set.
        This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.


        Parameters
        ----------
        parameter_type : Tuple[ :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterTypes`, Tuple[str, str] ]
            help
        Returns
        -------
        List[ Tuple[int, int] ]
            help





    )";




    } else if(name == "ObservationViabilityCalculator") {
         return R"(

        Template class for observation viability calculators.

        Template class for classes which conducts viability calculations on simulated observations.
        Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
        The user typically does not interact directly with this class.





     )";




    } else if(name == "ObservationViabilityCalculator.is_observation_viable" && variant==0) {
            return R"(

        Function to check whether an observation is viable.

        Function to check whether an observation is viable.
        The calculation is performed based on the given times and link end states.
        Note, that this function is called automatically during the simulation of observations.
        Direct calls to this function are generally not required.


        Parameters
        ----------
        link_end_states : List[ numpy.ndarray[numpy.float64[6, 1]] ]
            Vector of states of the link ends involved in the observation.
        link_end_times : List[float]
            Vector of times at the link ends involved in the observation.
        Returns
        -------
        bool
            True if observation is viable, false if not.





    )";




    } else if(name == "ObservationViabilityCalculator_1") {
         return R"(

        Class which conducts viability calculations on simulated observations of size 1.

        Class which conducts viability calculations on simulated observations of size 1, e.g. Doppler or range observations.
        Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
        The user typically does not interact directly with this class.





     )";





    } else if(name == "ObservationViabilityCalculator_2") {
         return R"(

        Class which conducts viability calculations on simulated observations.

        Class which conducts viability calculations on simulated observations of size 2, e.g. angular position observations.
        Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
        The user typically does not interact directly with this class.





     )";





    } else if(name == "ObservationViabilityCalculator_3") {
         return R"(

        Class which conducts viability calculations on simulated observations.

        Class which conducts viability calculations on simulated observations of size 3, e.g. Euler angle observations
        Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
        The user typically does not interact directly with this class.





     )";





    } else if(name == "ObservationViabilityCalculator_6") {
         return R"(

        Class which conducts viability calculations on simulated observations.

        Class which conducts viability calculations on simulated observations of size 6, e.g. (pseudo-) observations of the full cartesian state.
        Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
        The user typically does not interact directly with this class.





     )";





    } else if(name == "ObservationSimulator") {
         return R"(

        Class hosting the functionality for simulating observations.

        Class hosting the functionality for simulating a given observable over a defined link geometry.
        Instances of this class are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` objects upon instantiation of the :class:`~tudatpy.numerical_simulation.Estimator` class.





     )";





    } else if(name == "ObservationCollection") {
         return R"(

        Class collecting all observations and associated data for use in an estimation.

        Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
        this class is instantiated via a call to the :func:`~tudatpy.numerical_simulation.estimation.simulate_observations` function. More information is provided
        on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_





     )";


    } else if(name == "ObservationCollection.concatenated_times") {
         return R"(

        **read-only**

        Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )";


    } else if(name == "ObservationCollection.concatenated_observations") {
         return R"(

        **read-only**

        Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )";


    } else if(name == "ObservationCollection.concatenated_link_definition_ids") {
         return R"(

        **read-only**

        Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order of the present vector.

        :type: numpy.ndarray[ int ]
     )";


    } else if(name == "ObservationCollection.link_definition_ids") {
         return R"(

        **read-only**

        Dictionaty mapping a link end integer identifier to the specific link ends

        :type: dict[ int, dict[ LinkEndType, LinkEndId ] ]
     )";


    } else if(name == "ObservationCollection.observable_type_start_index_and_size") {
         return R"(

        **read-only**

        Dictionary defining per obervable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type

        :type: dict[ ObservableType, [ int, int ] ]
     )";


    } else if(name == "ObservationCollection.observation_set_start_index_and_size") {
         return R"(

        **read-only**

        The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the correspoding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per obervable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.

        :type: dict[ ObservableType, dict[ int, list[ int, int ] ] ]
     )";


    } else if(name == "ObservationCollection.sorted_observation_sets") {
         return R"(

        **read-only**

        The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_

        :type: dict[ ObservableType, dict[ int, list[ SingleObservationSet ] ] ]
     )";


    } else if(name == "ObservationCollection.observation_vector_size") {
         return R"(

        **read-only**

        Length of the total vector of observations

        :type: int
     )";




    } else if(name == "ObservationCollection.get_single_link_and_type_observations" && variant==0) {
            return R"(

        Function to get all observation sets for a given observable type and link definition.


        Parameters
        ----------
        observable_type : :class:`ObservableType`
            Observable type of which observations are to be simulated.
        link_ends : LinkDefinition
            Link ends for which observations are to be simulated.
        Returns
        -------
        list[ SingleObservationSet ]
            List of observation sets for given observable type and link definition.





    )";




    } else if(name == "SingleObservationSet") {
         return R"(

        Class collecting a single set of observations and associated data, of a given observable type, link ends, and ancilliary data.





     )";


    } else if(name == "SingleObservationSet.observable_type") {
         return R"(

        **read-only**

        Type of observable for which the object stores observations

        :type: ObservableType
     )";


    } else if(name == "SingleObservationSet.link_definition") {
         return R"(

        **read-only**

        Definition of the link ends for which the object stores observations

        :type: LinkDefinition
     )";


    } else if(name == "SingleObservationSet.reference_link_end") {
         return R"(

        **read-only**

        Reference link end for all stored observations

        :type: LinkEndType
     )";


    } else if(name == "SingleObservationSet.ancilliary_settings") {
         return R"(

        **read-only**

        Ancilliary settings all stored observations

        :type: ObservationAncilliarySimulationSettings
     )";


    } else if(name == "SingleObservationSet.concatenated_observations") {
         return R"(

        **read-only**

        Concatenated vector of all stored observations

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )";


    } else if(name == "SingleObservationSet.list_of_observations") {
         return R"(

        **read-only**

        List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1

        :type: list[ numpy.ndarray[numpy.float64[m, 1]] ]
     )";


    } else if(name == "SingleObservationSet.observation_times") {
         return R"(

        **read-only**

        Reference time for each of the observations in ``list_of_observations``

        :type: list[ float]
     )";


    } else if(name == "SingleObservationSet.observations_history") {
         return R"(

        **read-only**

        Dictionary of observations sorted by time. Created by making a dictionaty with ``observation_times`` as keys and  ``list_of_observations`` as values

        :type: dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
     )";





    } else if(name == "CombinedStateTransitionAndSensitivityMatrixInterface") {
         return R"(

        Class establishing an interface with the simulation's State Transition and Sensitivity Matrices.

        Class establishing an interface to the State Transition and Sensitivity Matrices.
        Instances of this class are instantiated automatically upon creation of :class:`~tudatpy.numerical_simulation.Estimator` objects,
        using the simulation information in the observation, propagation and integration settings that the :class:`~tudatpy.numerical_simulation.Estimator` instance is linked to.





     )";


    } else if(name == "CombinedStateTransitionAndSensitivityMatrixInterface.state_transition_size") {
         return R"(

        **read-only**

        Size of the (square) state transition matrix.

        :type: int
     )";


    } else if(name == "CombinedStateTransitionAndSensitivityMatrixInterface.sensitivity_size") {
         return R"(

        **read-only**

        Number of columns in the sensitivity matrix.

        :type: int
     )";


    } else if(name == "CombinedStateTransitionAndSensitivityMatrixInterface.full_parameter_size") {
         return R"(

        **read-only**

        Full amount of parameters w.r.t. which partials have been set up via State Transition and Sensitivity Matrices.

        :type: int
     )";




    } else if(name == "CombinedStateTransitionAndSensitivityMatrixInterface.state_transition_sensitivity_at_epoch" && variant==0) {
            return R"(

        Function to get the concatenated state transition and sensitivity matrix at a given time.

        Function to get the concatenated state transition and sensitivity matrix at a given time.
        Entries corresponding to parameters which are not active at the current arc are omitted.


        Parameters
        ----------
        time : float
            Time at which concatenated state transition and sensitivity matrix are to be retrieved.
        Returns
        -------
        numpy.ndarray[numpy.float64[m, n]]
            Concatenated state transition and sensitivity matrix at a given time.





    )";



    } else if(name == "CombinedStateTransitionAndSensitivityMatrixInterface.full_state_transition_sensitivity_at_epoch" && variant==0) {
            return R"(


        Parameters
        ----------
        time : float
            Time at which full concatenated state transition and sensitivity matrix are to be retrieved.
        Returns
        -------
        numpy.ndarray[numpy.float64[m, n]]
            Full concatenated state transition and sensitivity matrix at a given time.





    )";




    } else if(name == "EstimationConvergenceChecker") {
         return R"(

        Class defining the convergence criteria for an estimation.

        Class defining the convergence criteria for an estimation.
        The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker` factory function.





     )";





    } else if(name == "CovarianceAnalysisInput") {
         return R"(

        Class for defining all specific inputs to a covariance analysis.





     )";


    } else if(name == "CovarianceAnalysisInput.weight_matrix_diagonal") {
         return R"(

        **read-only**

        Complete diagonal of the weights matrix that is to be used

        :type: numpy.ndarray[numpy.float64[n, 1]]
     )";




    } else if(name == "CovarianceAnalysisInput.ctor" && variant==0) {
            return R"(

        Class constructor.

        Constructor through which the user can create instances of this class. Note that the weight are all initiated as 1.0, and the default settings of ``define_covariance_settings`` are used.


        Parameters
        ----------
        observations_and_times : ObservationCollection
            Total data structure of observations and associated times/link ends/type/etc.
        inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
            A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput`
            Instance of the :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput` class, defining the data and other settings to be used for the covariance analysis.





    )";



    } else if(name == "CovarianceAnalysisInput.set_constant_weight" && variant==0) {
            return R"(

        Function to set a constant weight matrix for all observables.

        Function to set a constant weight matrix for all observables.
        The weights are applied to all observations managed by the given PodInput object.


        Parameters
        ----------
        constant_weight : float
            Constant weight factor that is to be applied to all observations.
        Returns
        -------
        None
            Function modifies the object in-place.





    )";



    } else if(name == "CovarianceAnalysisInput.set_constant_weight_per_observable" && variant==0) {
            return R"(

        Function to set a constant weight matrix for a given type of observable.

        Function to set a constant weight matrix for a given type of observable.
        The weights are applied to all observations of the observable type specified by the `weight_per_observable` parameter.


        Parameters
        ----------
        constant_weight : Dict[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`, float ]
            Constant weight factor that is to be applied to all observations.
        Returns
        -------
        None
            Function modifies the object in-place.





    )";



    } else if(name == "CovarianceAnalysisInput.define_covariance_settings" && variant==0) {
            return R"(

        Function to define specific settings for covariance analysis process

        Function to define specific settings for covariance analysis process


        Parameters
        ----------
        reintegrate_equations : bool, default = True
            Boolean denoting whether the dynamics and variational equations are to be reintegrated
            or if existing values are to be used to perform first iteration.

        reintegrate_variational_equations : bool, default = True
            Boolean denoting whether the variational equations are to be reintegrated during estimation 
            (if this is set to False, and ``reintegrate_equations`` to true, only the dynamics are re-integrated)

        save_design_matrix : bool, default = True
            Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
            :math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.

        print_output_to_terminal : bool, default = True
            Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.

        Returns
        -------
        None
            Function modifies the object in-place.





    )";




    } else if(name == "EstimationInput") {
         return R"(

        Class for defining all inputs to the estimation.





     )";




    } else if(name == "EstimationInput.ctor" && variant==0) {
            return R"(

        Class constructor.

        Constructor through which the user can create instances of this class.


        Parameters
        ----------
        observations_and_times : ObservationCollection
            Total data structure of observations and associated times/link ends/type/etc.
        inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
            A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
        convergence_checker : :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker`, default = :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker`
            Object defining when the estimation is converged.
        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.EstimationInput`
            Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationInput` class, defining the data and other settings to be used for the estimation.





    )";



    } else if(name == "EstimationInput.define_estimation_settings" && variant==0) {
            return R"(

        Function to define specific settings for the estimation process

        Function to define specific settings for covariance analysis process


        Parameters
        ----------
        reintegrate_equations_on_first_iteration : bool, default = True
            Boolean denoting whether the dynamics and variational equations are to be reintegrated
            or if existing values are to be used to perform first iteration.

        reintegrate_variational_equations : bool, default = True
            Boolean denoting whether the variational equations are to be reintegrated during estimation 
            (if this is set to False, and ``reintegrate_equations_on_first_iteration`` to true, only the dynamics are re-integrated)

        save_design_matrix : bool, default = True
            Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
            :math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.

        print_output_to_terminal : bool, default = True
            Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.

        save_residuals_and_parameters_per_iteration : bool, default = True
            Boolean denoting whether the residuals and parameters from the each iteration are to be saved.

        save_state_history_per_iteration : bool, default = False
            Boolean denoting whether the state history and dependent variables are to be saved on each iteration.

        Returns
        -------
        None
            Function modifies the object in-place.





    )";




    } else if(name == "CovarianceAnalysisOutput") {
         return R"(

        Class collecting all outputs from the covariance analysis process.





     )";


    } else if(name == "CovarianceAnalysisOutput.inverse_covariance") {
         return R"(

        **read-only**

        (Unnormalized) inverse estimation covariance matrix :math:`\mathbf{P}^{-1}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )";


    } else if(name == "CovarianceAnalysisOutput.covariance") {
         return R"(

        **read-only**

        (Unnormalized) estimation covariance matrix :math:`\mathbf{P}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )";


    } else if(name == "CovarianceAnalysisOutput.inverse_normalized_covariance") {
         return R"(

        **read-only**

        Normalized inverse estimation covariance matrix :math:`\mathbf{\tilde{P}}^{-1}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )";


    } else if(name == "CovarianceAnalysisOutput.normalized_covariance") {
         return R"(

        **read-only**

        Normalized estimation covariance matrix :math:`\mathbf{\tilde{P}}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )";


    } else if(name == "CovarianceAnalysisOutput.formal_errors") {
         return R"(

        **read-only**

        Formal error vector :math:`\boldsymbol{\sigma}` of the estimation result (e.g. square root of diagonal entries of covariance)s

        :type: numpy.ndarray[numpy.float64[m, 1]]s
     )";


    } else if(name == "CovarianceAnalysisOutput.correlations") {
         return R"(

        **read-only**

        Correlation matrix of the estimation result. Entry :math:`i,j` is equal to :math:`P_{i,j}/(\sigma_{i}\sigma_{j})`

        :type: numpy.ndarray[numpy.float64[m, m]]
     )";


    } else if(name == "CovarianceAnalysisOutput.design_matrix") {
         return R"(

        **read-only**

        Matrix of unnormalized partial derivatives :math:`\mathbf{H}=\frac{\partial\mathbf{h}}{\partial\mathbf{p}}`.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )";


    } else if(name == "CovarianceAnalysisOutput.normalized_design_matrix") {
         return R"(

        **read-only**

        Matrix of normalized partial derivatives :math:`\tilde{\mathbf{H}}`.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )";


    } else if(name == "CovarianceAnalysisOutput.weighted_design_matrix") {
         return R"(

        **read-only**

        Matrix of weighted partial derivatives, equal to :math:`\mathbf{W}^{1/2}{\mathbf{H}}`

        :type: numpy.ndarray[numpy.float64[m, n]]
     )";


    } else if(name == "CovarianceAnalysisOutput.weighted_normalized_design_matrix") {
         return R"(

        **read-only**

        Matrix of weighted, normalized partial derivatives, equal to :math:`\mathbf{W}^{1/2}\tilde{\mathbf{H}}`

        :type: numpy.ndarray[numpy.float64[m, n]]
     )";


    } else if(name == "CovarianceAnalysisOutput.normalization_terms") {
         return R"(

        **read-only**

        Vector of normalization terms used for covariance and design matrix

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )";





    } else if(name == "EstimationOutput") {
         return R"(

        Class collecting all outputs from the iterative estimation process.





     )";


    } else if(name == "EstimationOutput.residual_history") {
         return R"(

        **read-only**

        Residual vectors, concatenated per iteration into a matrix; the :math:`i^{th}` column has the residuals from the :math:`i^{th}` iteration.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )";


    } else if(name == "EstimationOutput.parameter_history") {
         return R"(

        **read-only**

        Parameter vectors, concatenated per iteration into a matrix. Column 0 contains pre-estimation values. The :math:`(i+1)^{th}` column has the residuals from the :math:`i^{th}` iteration.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )";


    } else if(name == "EstimationOutput.simulation_results_per_iteration") {
         return R"(

        **read-only**

        List of complete numerical propagation results, with the :math:`i^{th}` entry of thee list thee results of the :math:`i^{th}` propagation

        :type: list[SimulationResults]
     )";


    } else if(name == "EstimationOutput.final_residuals") {
         return R"(

        **read-only**

        Vector of post-fit observation residuals, for the iteration with the lowest rms residuals.

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )";






    } else if(name == "simulate_observations" && variant==0) {
        return R"(
        
Function to simulate observations.

Function to simulate observations from set observation simulators and observation simulator settings.
Automatically iterates over all provided observation simulators, generating the full set of simulated observations.


Parameters
----------
observation_to_simulate : List[ :class:`ObservationSimulationSettings` ]
    List of settings objects, each object providing the observation time settings for simulating one type of observable and link end set.

observation_simulators : List[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
    List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one type of observable and link end set.

bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation.ObservationCollection`
    Object collecting all products of the observation simulation.






    )";



    } else if(name == "compute_target_angles_and_range" && variant==0) {
        return R"(
        
Function to compute the azimuth angle, elevation angle and range at a ground station.

Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
as observed at that station   


Parameters
----------
bodies : SystemOfBodies
    System of bodies that defines the full physical environment

station_id : tuple[ str, str]
    Identifier for the observing station, as a pair of strings: the body name and the station name.

target_body : str
    Name of body which is observed by ground station

observation_times : list[float]
    List of times at which the ground station observations are to be analyzed

is_station_transmitting : Bool
    Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station. 
    This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target

Returns
-------
dict[float,numpy.ndarray[numpy.float64[3, 1]]]
    Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range






    )";



    } else if(name == "propagate_covariance" && variant==0) {
        return R"(
        
Function to propagate system covariance through time.

Function to propagate the covariance of a given system through time.
The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


Parameters
----------
initial_covariance : numpy.ndarray[numpy.float64[m, n]]
    System covariance matrix (symmetric and positive semi-definite) at initial time.
    Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
    Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

output_times : List[ float ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the covariance propagation,
    which always adheres to the integrator settings that the `state_transition_interface` links to.
    Output times which do not coincide with integration time steps are calculated via interpolation.

Returns
-------
Dict[ float, numpy.ndarray[numpy.float64[m, n]] ]
    Dictionary reporting the propagated covariances at each output time.






    )";



    } else if(name == "propagate_formal_errors" && variant==0) {
        return R"(
        
Function to propagate system formal errors through time.

Function to propagate the formal errors of a given system through time.
Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


Parameters
----------
initial_covariance : numpy.ndarray[numpy.float64[m, n]]
    System covariance matrix (symmetric and positive semi-definite) at initial time.
    Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
    Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

output_times : List[ float ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the covariance propagation,
    which always adheres to the integrator settings that the `state_transition_interface` links to.
    Output times which do not coincide with integration time steps are calculated via interpolation.

Returns
-------
Dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
    Dictionary reporting the propagated formal errors at each output time.






    )";



    } else if(name == "estimation_convergence_checker" && variant==0) {
        return R"(
        
Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object.

Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object, which is required for defining the convergence criteria of an estimation.


Parameters
----------
maximum_iterations : int, default = 5
    Maximum number of allowed iterations for estimation.
minimum_residual_change : float, default = 0.0
    Minimum required change in residual between two iterations.
minimum_residual : float, default = 0.0
    Minimum value of observation residual below which estimation is converged.
number_of_iterations_without_improvement : int, default = 2
    Number of iterations without reduction of residual.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` class, defining the convergence criteria for an estimation.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace environment {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "AerodynamicsReferenceFrames") {
         return R"(

        Enumeration of reference frame identifiers typical for aerodynamic calculations.

        Enumeration of reference frame identifiers typical for aerodynamic calculations. Note that the frames are also defined
        in the absence of any aerodynamic forces and/or atmosphere. They define frames of a body w.r.t. a central body, with
        the details given by Mooij (1994). The chain of frames starts from the inertial frame, to the frame fixed to the
        central body (corotating), to the vertical frame (defined by the body's relative position), the trajectory and aerodynamic frames
        (defined by the body's relative velocity) and finally the body's own body-fixed frame.





     )";


    } else if(name == "AerodynamicsReferenceFrames.inertial_frame") {
         return R"(
The global orientation (which is by definition inertial).
     )";


    } else if(name == "AerodynamicsReferenceFrames.corotating_frame") {
         return R"(
The body-fixed frame of the central body.
     )";


    } else if(name == "AerodynamicsReferenceFrames.vertical_frame") {
         return R"(
Frame with z-axis pointing towards origin of central body, the x-axis lies in the meridian plane and points towards the central-body-fixed z-axis (the y-axis completes the frame).
     )";


    } else if(name == "AerodynamicsReferenceFrames.trajectory_frame") {
         return R"(
The (airspeed-based) trajectory frame has the x-axis in the direction of the velocity vector relative to the atmosphere (airspeed-based velocity vector), z-axis lies in the vertical plane and points downwards (the y-axis completes the frame).
     )";


    } else if(name == "AerodynamicsReferenceFrames.aerodynamic_frame") {
         return R"(
The (airspeed-based) aerodynamic frame has the x-axis in the direction of the velocity vector relative to the atmosphere (airspeed-based velocity vector), z-axis co-linear with the aerodynamic lift vector, pointing in the opposite direction (the y-axis completes the frame)..
     )";


    } else if(name == "AerodynamicsReferenceFrames.body_frame") {
         return R"(
The body-fixed frame of the body itself.
     )";



    } else if(name == "AerodynamicCoefficientFrames") {
         return R"(

        Enumeration of reference frames used for definition of aerodynamic coefficients.

        Enumeration of reference frames used for definition of aerodynamic coefficients. There is a partial overlap between this enum
        and the :class:`~tudatpy.numerical_simulation.environment.AerodynamicsReferenceFrames`. This enum combines a subset of those
        frames (which are typically used for aerodynamic coefficient definition), and a swap in sign. For instance, aerodynamic
        force coefficients are often defined positive along *negative* axes of the aerodynamic frame (drag, side force and lift coefficients)





     )";


    } else if(name == "AerodynamicCoefficientFrames.positive_body_fixed_frame_coefficients") {
         return R"(
The coefficients are defined in the body-fixed frame, with the directions the same as the body-fixed axes. For aerodynamic forces and moments, this results in the typical :math:`C_{x}, C_{y}, C_{y}` (force) and :math:`C_{l}, C_{m}, C_{n}` (moment) coefficients
     )";


    } else if(name == "AerodynamicCoefficientFrames.negative_body_fixed_frame_coefficients") {
         return R"(
Same as ``positive_body_fixed_frame_coefficients``, but opposite in direction (so axes along negative body-fixed frame axes)
     )";


    } else if(name == "AerodynamicCoefficientFrames.positive_aerodynamic_frame_coefficients") {
         return R"(
Same as ``negative_aerodynamic_frame_coefficients``, but opposite in direction (so axes along positive aerodynamic frame axes)
     )";


    } else if(name == "AerodynamicCoefficientFrames.negative_aerodynamic_frame_coefficients") {
         return R"(
The coefficients are defined in aerodynamic frame, with the directions the same as the negative axes. For aerodynamic forces, this results in the typical :math:`C_{D}, C_{S}, C_{D}` force coefficients
     )";



    } else if(name == "AerodynamicCoefficientsIndependentVariables") {
         return R"(

        Enumeration of the independent variables that can be used to compute aerodynamic coefficients.





     )";


    } else if(name == "AerodynamicCoefficientsIndependentVariables.mach_number_dependent") {
         return R"(
Mach number of the propagated vehicle.
     )";


    } else if(name == "AerodynamicCoefficientsIndependentVariables.angle_of_attack_dependent") {
         return R"(
Angle of attack of the propagated vehicle.
     )";


    } else if(name == "AerodynamicCoefficientsIndependentVariables.sideslip_angle_dependent") {
         return R"(
Sideslip angle of the propagated vehicle.
     )";


    } else if(name == "AerodynamicCoefficientsIndependentVariables.altitude_dependent") {
         return R"(
Altitude of the propagated vehicle.
     )";


    } else if(name == "AerodynamicCoefficientsIndependentVariables.time_dependent") {
         return R"(
Current simulation epoch.
     )";


    } else if(name == "AerodynamicCoefficientsIndependentVariables.control_surface_deflection_dependent") {
         return R"(
Angle of deflection of the control surface of the propagated vehicle.
     )";


    } else if(name == "AerodynamicCoefficientsIndependentVariables.undefined_independent_variable") {
         return R"(
Can be used for a custom coefficient interface with other variables, at the expense of being able to use the FlightConditions class to automatically updates the aerodynamic coefficients during propagation.
     )";




    } else if(name == "AerodynamicCoefficientInterface") {
         return R"(

        Base class for computing the current aerodynamic coefficients of the body


        Base class for computing the current aerodynamic coefficients of the body. The implementation of the computation
        depends on the choice of aerodynamic coefficient model (see :ref:`\`\`aerodynamic_coefficients\`\`` for available options).
        During the propagation, this object is automatically updated to the current state by the :class:`~AtmosphericFlightConditions` object.
        The user may override the current aerodynamic coefficients when using, for instance, a custom aerodynamic guidance model
        (see `here <https://docs.tudat.space/en/latest/_src_getting_started/_src_examples/notebooks/propagation/reentry_trajectory.html>`_ for an example).
        using the member functions of this class.





     )";


    } else if(name == "AerodynamicCoefficientInterface.reference_area") {
         return R"(

        **read-only**

        The aerodynamic reference area :math:`A` of the coefficients


        :type: float
     )";


    } else if(name == "AerodynamicCoefficientInterface.current_force_coefficients") {
         return R"(

        **read-only**

        The current aerodynamic force coefficients, in the frame defined by the :attr:`~force_coefficient_frame` attribute,
        as computed by the last call to the :meth:`~update_coefficients` function.


        :type: np.ndarray
     )";


    } else if(name == "AerodynamicCoefficientInterface.current_moment_coefficients") {
         return R"(

        **read-only**

        The current aerodynamic moment coefficients, in the frame defined by the :attr:`~moment_coefficient_frame` attribute,
        as computed by the last call to the :meth:`~update_coefficients` function.   


        :type: np.ndarray
     )";


    } else if(name == "AerodynamicCoefficientInterface.current_coefficients") {
         return R"(

        **read-only**

        Concatenation of :attr:`~current_force_coefficients` and :attr:`~current_moment_coefficients`


        :type: np.ndarray
     )";


    } else if(name == "AerodynamicCoefficientInterface.force_coefficient_frame") {
         return R"(

        **read-only**

        Reference frame in which the  :attr:`~current_force_coefficients` are defined


        :type: AerodynamicCoefficientFrames
     )";


    } else if(name == "AerodynamicCoefficientInterface.moment_coefficient_frame") {
         return R"(

        **read-only**

        Reference frame in which the  :attr:`~current_moment_coefficients` are defined


        :type: AerodynamicCoefficientFrames
     )";


    } else if(name == "AerodynamicCoefficientInterface.independent_variable_names") {
         return R"(

        **read-only**

        List of independent variables from which the aerodynamic coefficients are computed (e.g. required input to :meth:`~update_coefficients` function).


        :type: list[AerodynamicCoefficientsIndependentVariables]
     )";


    } else if(name == "AerodynamicCoefficientInterface.current_control_surface_free_force_coefficients") {
         return R"(

        **read-only**

        Same as :attr:`current_force_coefficients`, but without contribution (if any) from control surfaces


        :type: np.ndarray
     )";


    } else if(name == "AerodynamicCoefficientInterface.current_control_surface_free_moment_coefficients") {
         return R"(

        **read-only**

        Same as :attr:`current_moment_coefficients`, but without contribution (if any) from control surfaces


        :type: np.ndarray
     )";


    } else if(name == "AerodynamicCoefficientInterface.control_surface_independent_variable_names") {
         return R"(

        **read-only**

        List of independent variables from which the aerodynamic coefficients of each control surface are computed, with dictionary key being the control surface name (e.g. required input to :meth:`~update_full_coefficients` function).


        :type: dict[str,list[AerodynamicCoefficientsIndependentVariables]]
     )";




    } else if(name == "AerodynamicCoefficientInterface.current_control_surface_force_coefficient_increment" && variant==0) {
            return R"(

        Function to get the contribution from a single control surface to the aerodynamic force coefficient, as compute by last call to :meth:`~update_full_coefficients`



        Parameters
        ----------
        control_surface_name : str
            The name of the control surface for which the contribution is to be retrieved

        Returns
        -------
        np.ndarray
            Contribution from the requested control surface to the aerodynamic force coefficient





    )";



    } else if(name == "AerodynamicCoefficientInterface.current_control_surface_moment_coefficient_increment" && variant==0) {
            return R"(

        Function to get the contribution from a single control surface to the aerodynamic moment coefficients, as compute by last call to :meth:`~update_full_coefficients`



        Parameters
        ----------
        control_surface_name : str
            The name of the control surface for which the contribution is to be retrieved

        Returns
        -------
        np.ndarray
            Contribution from the requested control surface to the aerodynamic moment coefficients





    )";



    } else if(name == "AerodynamicCoefficientInterface.update_coefficients" && variant==0) {
            return R"(

        Function to update the aerodynamic coefficients of the body only


        Function to update the aerodynamic coefficients of the body only (without the control surface contribution),
        based on the current state. This function may be called by the user, but will set *only* the
        :attr:`~current_force_coefficients` and :attr:`~current_moment_coefficients` (while leaving the
        :attr:`~current_control_surface_free_force_coefficients` and :attr:`~current_control_surface_free_moment_coefficients` unchanged)


        Parameters
        ----------
        independent_variables : list[float]
            List of inputs from which the aerodynamic coefficients are to be computed, with each entry corresponding to the
            value of the physical variable defined by the :attr:`independent_variable_names` attribute.

        time : float
            Current time (in seconds since J2000)

        Returns
        -------
        np.ndarray
            Contribution from the requested control surface to the aerodynamic moment coefficients





    )";



    } else if(name == "AerodynamicCoefficientInterface.update_full_coefficients" && variant==0) {
            return R"(

        Function to update the aerodynamic coefficients, from both the body and its control surfaces


        Function to update the aerodynamic coefficients of both the body and its control surfaces,
        based on the current state. This function will call the :meth:`~update_coefficients` function to update the body coefficients.
        This function may be called by the user, and will set the following attributes:
        :attr:`~current_force_coefficients`, :attr:`~current_moment_coefficients` , 
        :attr:`~current_control_surface_free_force_coefficients` and :attr:`~current_control_surface_free_moment_coefficients`.
        In addition, it will modify the coefficients returned by the :meth:`~current_control_surface_force_coefficient_increment` and
        :meth:`~current_control_surface_moment_coefficient_increment` functions


        Parameters
        ----------
        independent_variables : list[float]
            List of inputs from which the aerodynamic coefficients of the body are to be computed, with each entry corresponding to the
            value of the physical variable defined by the :attr:`independent_variable_names` attribute.

        control_surface_independent_variables : dict[str,list[float]]
            List of inputs from which the control surface aerodynamic coefficients are to be computed (with dictionary key the control surface name), 
            with each entry corresponding to the
            value of the physical variable defined by the :attr:`control_surface_independent_variable_names` attribute.

        time : float
            Current time (in seconds since J2000)

        check_force_contribution : bool, default = True
            Boolean that determines if the force contribution to the aerodynamic moments should be added. Note that this input is
            only used if the :attr:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings.add_force_contribution_to_moments` attribute is set to True.





    )";




    } else if(name == "HypersonicLocalInclinationAnalysis") {
         return R"(

        Class to calculate the hypersonic aerodynamic coefficients using local inclination methods.

        Class to calculate the hypersonic aerodynamic coefficients using local inclination methods. This class
        implements several local inclination methods, and calculates the aerodynamic force and moment coefficients
        as a function of angle of attack, sideslip angle and Mach number





     )";




    } else if(name == "HypersonicLocalInclinationAnalysis.ctor" && variant==0) {
            return R"(

        Class constructor, taking the shape of the vehicle, and various analysis options as input.


        Parameters
        ----------
        independent_variable_points : list[list[float]]
            List containing three lists, with each sublist containing the data points of each of the 
            independent variables for the coefficient generation. The physical meaning of each of the 
            three independent variables is: 0 = mach number, 1 = angle of attack, 2 = angle of sideslip. 
            Each of the subvectors must be sorted in ascending order.

        body_shape : SurfaceGeometry
            Class that defines the shape of the vehicle as a continuous surface. The local inclination analysis
            discretizes the surface of the vehicle into quadrilateral panels, defined by the other inputs to
            this constructor. In case the :class:`tudat.geometry.SurfaceGeometry` object is made up of multiple
            sub-shapes, different settings may be used for each

        number_of_lines : List[ float ]
            Number of discretization points in the first independent surface variable of each of the subparts of body_shape. 	
            The size of this list should match the number of parts of which the body_shape is composed. The first independent
            variable of a subpart typically runs along the longitudinal vehicle direction

        number_of_points : List[ float ]
            Number of discretization points in the second independent surface variable of each of the subparts of body_shape. 	
            The size of this list should match the number of parts of which the body_shape is composed. The first independent
            variable of a subpart typically runs along the lateral vehicle direction

        invert_orders : List[ bool ]
            Booleans to denote whether the surface normals of the panels of each discretized body_shape subpart are to be inverted
            (i.e. inward-facing->outward facing or vice versa). The size of this list should match the number of parts of which the body_shape is composed. 

        selected_methods : List[ List[ int ] ]
            Double list of selected local inclination methods, the first index (outer list) represents compression or expansion (0 and 1), 
            the second index (inner list) denotes the vehicle part index. The size of this inner list should match the number of parts of which the body_shape is composed. 
            The int defining the method type is interpreted as follows. 
            For the compression methods, the following are available:
            *  0: Newtonian Method.
            *  1: Modified Newtonian.
            *  2 and 3: not available at this moment.
            *  4: Tangent-wedge method.
            *  5: Tangent-cone method.
            *  6: Modified Dahlem-Buck method.
            *  7: VanDyke unified pressure method.
            *  8: Smyth Delta Wing method.
            *  9: Hankey flat surface method
            The expansion method has the following options:
            *  0: Vacuum Pressure coefficient method.
            *  1: Zero Pressure function.
            *  4: High Mach base pressure method.
            *  3 or 5: Prandtl-Meyer method.
            *  6: ACM empirical pressure coefficient.

        reference_area : float
            Reference area used to non-dimensionalize aerodynamic forces and moments.

        moment_reference_point : numpy.ndarray
            Reference point wrt which aerodynamic moments are calculated.

        save_pressure_coefficients : Bool
            Boolean denoting whether to save the pressure coefficients that are computed to files





    )";




    } else if(name == "FlightConditions") {
         return R"(

        Object that calculates various state-derived quantities typically
        relevant for flight dynamics.


        Object that calculates various state-derived quantities typically
        relevant for flight dynamics, such as latitude, longitude,
        altitude, etc. It also contains an
        :py:class:`~AerodynamicAngleCalculator` that computes derived
        angles (flight path, heading angle, etc.). This object is limited
        to non-atmospheric flight. For flight through Body objects
        endowed with an atmosphere model, the derived class
        :py:class:`~AtmosphericFlightConditions` is used. This object is
        stored inside a Body object, and represents the flight conditions
        of a single body w.r.t. a single central body.





     )";


    } else if(name == "FlightConditions.aerodynamic_angle_calculator") {
         return R"(

        **read-only**

        The object that is responsible for computing various relevant
        flight dynamics angles and frame rotations.


        :type: AerodynamicAngleCalculator
     )";


    } else if(name == "FlightConditions.longitude") {
         return R"(

        **read-only**

        The body-fixed longitude of the body w.r.t. its central body.


        :type: float
     )";


    } else if(name == "FlightConditions.latitude") {
         return R"(

        **read-only**

        The body-fixed geographic latitude of the body w.r.t. its
        central body.


        :type: float
     )";


    } else if(name == "FlightConditions.geodetic_latitude") {
         return R"(

        **read-only**

        The body-fixed geodetic latitude of the body w.r.t. its central
        body, using an :py:class:`~OblateSpheroidBodyShapeModel` of the
        central body. If no such model is defined in the central body,
        this attribute equals the latitude.


        :type: float
     )";


    } else if(name == "FlightConditions.altitude") {
         return R"(

        **read-only**

        The altitude of this body w.r.t. the central body (using the
        :py:class:`~ShapeModel` of the central body).


        :type: float
     )";


    } else if(name == "FlightConditions.time") {
         return R"(

        **read-only**

        The current time, at which this object was last updated


        :type: float
     )";


    } else if(name == "FlightConditions.body_centered_body_fixed_state") {
         return R"(

        **read-only**

        Cartesian translational state, expressed in a frame centered
        at, and fixed to, the central body. Note that, due to the
        rotation of the central body, the norm of the body-fixed,
        body-centered, velocity differs from the norm of the inertial
        body-centered velocity. 


        :type: numpy.ndarray
     )";





    } else if(name == "AtmosphericFlightConditions") {
         return R"(

        Object that calculates various state-derived quantities typically
        relevant for flight dynamics, for flight in an atmosphere.


        Object that calculates various state-derived quantities typically
        relevant for flight dynamics, for flight in an atmosphere, such
        as latitude,  longitude, altitude, density, Mach number etc. It
        also contains an ``AerodynamicAngleCalculator`` that computes
        derived angles (flight path, heading angle, etc.). This object is
        derived from ``FlightConditions``, which performs computations for
        non-atmospheric flight only. This object is stored inside a Body
        object, and represents the flight conditions of a single body
        w.r.t. a single central body.





     )";


    } else if(name == "AtmosphericFlightConditions.density") {
         return R"(

        **read-only**

        The freestream atmospheric density at the body's current
        location.


        :type: float
     )";


    } else if(name == "AtmosphericFlightConditions.temperature") {
         return R"(

        **read-only**

        The freestream atmospheric temperature at the body's current
        location.


        :type: float
     )";


    } else if(name == "AtmosphericFlightConditions.dynamic_pressure") {
         return R"(

        **read-only**

        The freestream atmospheric dynamic pressure at the body's
        current location.


        :type: float
     )";


    } else if(name == "AtmosphericFlightConditions.pressure") {
         return R"(

        **read-only**

        The freestream atmospheric static pressure at the body's
        current location.


        :type: float
     )";


    } else if(name == "AtmosphericFlightConditions.speed_of_sound") {
         return R"(

        **read-only**

        The freestream atmospheric speed of sound at the body's current
        location.


        :type: float
     )";


    } else if(name == "AtmosphericFlightConditions.airspeed") {
         return R"(

        **read-only**

        The airspeed of the body w.r.t. the atmosphere.


        :type: float
     )";


    } else if(name == "AtmosphericFlightConditions.airspeed_velocity") {
         return R"(

        **read-only**

        The velocity vector of the body w.r.t. the freestream
        atmosphere (e.g. vectorial counterpart of airspeed).


        :type: numpy.ndarray
     )";


    } else if(name == "AtmosphericFlightConditions.mach_number") {
         return R"(

        **read-only**

        The freestream Mach number of the body.


        :type: float
     )";


    } else if(name == "AtmosphericFlightConditions.aero_coefficient_independent_variables") {
         return R"(

        **read-only**

        List of current values of independent variables of aerodynamic
        coefficients. This list is only defined if the body has an
        :py:class:`~AerodynamicCoefficientInterface` that has
        dependencies on environmental variables (e.g. Mach number,
        angle of attack, etc.).


        :type: numpy.ndarray
     )";


    } else if(name == "AtmosphericFlightConditions.control_surface_aero_coefficient_independent_variables") {
         return R"(

        **read-only**

        List of lists current values of independent variables of
        aerodynamic coefficients for control surfaces. The outer list
        defines the control surface, the inner list the values of the
        independent variables. This list is only defined if the body
        has an :py:class:`~AerodynamicCoefficientInterface` with
        control surfaces that have dependencies on environmental
        variables (e.g. Mach number, angle of attack, etc.).


        :type: numpy.ndarray
     )";


    } else if(name == "AtmosphericFlightConditions.aerodynamic_coefficient_interface") {
         return R"(

        **read-only**

        Object extracted from the same Body object as this
        :py:class:`~AtmosphericFlightConditions` object, which defines
        the aerodynamic coefficients.


        :type: AerodynamicCoefficientInterface
     )";





    } else if(name == "AerodynamicAngleCalculator") {
         return R"(

        Object to calculate (aerodynamic) orientation angles, and frame transformations,
        from current vehicle state.


        Object to calculate (aerodynamic) orientation angles (list given by the :class:`~AerodynamicsReferenceFrameAngles` enum)
        and transformations between frames (list given by the :class:`~AerodynamicsReferenceFrames` enum) from current vehicle state.





     )";




    } else if(name == "AerodynamicAngleCalculator.get_rotation_matrix_between_frames" && variant==0) {
            return R"(

        Function to get the rotation matrix between two frames. 


        Function to get the rotation matrix between two frames. This function
        is meant to be used only *during* a numerical propagation, in particular
        for the definition of a custom (e.g. guidance) model.


        Parameters
        ----------
        original_frame : AerodynamicsReferenceFrames
            The frame :math:`A` from which the rotation matrix is to be calculated

        target_frame : AerodynamicsReferenceFrames
            The frame :math:`B` to which the rotation matrix is to be calculated

        Returns
        -------
        np.ndarray
            Rotation matrix :math:`\mathbf{R}^{B/A}` from frame :math:`A` to frame `B`





    )";



    } else if(name == "AerodynamicAngleCalculator.get_angle" && variant==0) {
            return R"(

        Function to get a single orientation angle


        Function to get a single orientation angle. This function
        is meant to be used only *during* a numerical propagation, in particular
        for the definition of a custom (e.g. guidance) model.


        Parameters
        ----------
        original_frame : AerodynamicsReferenceFrameAngles
            The identifier for the angle that is to be returnd

        Returns
        -------
        double
            Value of requested angle





    )";




    } else if(name == "RotationalEphemeris") {
         return R"(

        Object that stores the rotational state of the bodies.


        Object that stores the rotational state of the bodies. This object can be used to calculate rotation matrices,
        which are used to transform coordinates between reference frames.





     )";


    } else if(name == "RotationalEphemeris.body_fixed_frame_name") {
         return R"(

        **read-only**

        The identifier of the body-fixed frame, used in other parts of the simulation to identify it.


        :type: str
     )";


    } else if(name == "RotationalEphemeris.inertial_frame_name") {
         return R"(

        **read-only**

        The identifier of the inertial frame, used in other parts of the simulation to identify it.


        :type: str
     )";




    } else if(name == "RotationalEphemeris.body_fixed_to_inertial_rotation" && variant==0) {
            return R"(

        Function to get rotation matrix from body-fixed frame to inertial frame over time.


        Function to get rotation matrix from body-fixed (target) frame to inertial (base) frame over time. 
        The calculation of this rotation matrix depends on the specific rotation model that has been defined,
        either from an a priori definition (see :ref:`\`\`rotation_model\`\`` submodule) or from processing
        the results of propagation of the rotational equations of motion.


        Parameters
        ----------
        current_time : float
            The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix is evaluated

        Returns
        -------
        np.ndarray
            Rotation matrix :math:`\mathbf{R}^{I/B}` from body-fixed frame :math:`B` to inertial frame `I`





    )";



    } else if(name == "RotationalEphemeris.inertial_to_body_fixed_rotation" && variant==0) {
            return R"(

        Function to get rotation matrix from inertial frame to body-fixed frame over time.


        Function computes the inverse (equal to transpose) rotation of the ``body_fixed_to_inertial_rotation`` function.


        Parameters
        ----------
        current_time : float
            The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix is evaluated





    )";



    } else if(name == "RotationalEphemeris.time_derivative_body_fixed_to_inertial_rotation" && variant==0) {
            return R"(

        Function to get time derivative of rotation matrix from body-fixed frame to inertial frame over time.


        Function to get time derivative of rotation matrix from body-fixed frame to inertial frame over time (see ``body_fixed_to_inertial_rotation``), 
        denoted :math:`\dot{\mathbf{R}}^{(I/B)}`,


        Parameters
        ----------
        current_time : float
            The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix derivative is evaluated

        Returns
        -------
        np.ndarray
            Rotation matrix :math:`\dot{\mathbf{R}}^{I/B}` from body-fixed frame :math:`B` to inertial frame `I`





    )";



    } else if(name == "RotationalEphemeris.time_derivative_inertial_to_body_fixed_rotation" && variant==0) {
            return R"(

        Function to get time derivative of rotation matrix from inertial frame to body-fixed frame over time.


        Function to get time derivative of rotation matrix from inertial frame to body-fixed frame over time (see ``inertial_to_body_fixed_rotation``), 
        denoted :math:`\dot{\mathbf{R}}^{(B/I)}`,


        Parameters
        ----------
        current_time : float
            The time (in seconds since epoch J2000, TDB time scale) at which the rotation matrix derivative is evaluated

        Returns
        -------
        np.ndarray
            Rotation matrix :math:`\dot{\mathbf{R}}^{B/I}` from inertial frame `I` to body-fixed frame :math:`B`





    )";



    } else if(name == "RotationalEphemeris.angular_velocity_in_body_fixed_frame" && variant==0) {
            return R"(

        Function to get the body's angular velocity vector, expressed in the body-fixed frame.


        Function to get the body's angular velocity vector :math:`\boldsymbol{\omega}^{(B)}`, expressed in the body-fixed frame :math:`B`.
        The calculation of the angular velocity depends on the specific rotation model that has been defined,
        either from an a priori definition (see :ref:`\`\`rotation_model\`\`` submodule) or from processing
        the results of propagation of the rotational equations of motion. 
        Note that when numerically propagating rotational dynamics, this angular velocity vector is typically directly defined
        in the last three entries of the state vector.


        Parameters
        ----------
        current_time : float
            The time (in seconds since epoch J2000, TDB time scale) at which the angular velocity vector is evaluated

        Returns
        -------
        np.ndarray
            Angular velocity vector of the body  :math:`\boldsymbol{\omega}^{(B)}` expressed in the body-fixed frame :math:`B`





    )";



    } else if(name == "RotationalEphemeris.angular_velocity_in_inertial_frame" && variant==0) {
            return R"(

        Function to get the body's angular velocity vector, expressed in the inertial frame.


        Function to get the body's angular velocity vector :math:`\boldsymbol{\omega}^{(I)}`, expressed in the body-fixed frame :math:`I`.
        This quantity is computed from :math:`\mathbf{R}^{I/B}\boldsymbol{\omega}^{(B)}`, see the ``angular_velocity_in_body_fixed_frame`` and
        ``body_fixed_to_inertial_rotation`` functions.


        Parameters
        ----------
        current_time : float
            The time (in seconds since epoch J2000, TDB time scale) at which the angular velocity vector is evaluated

        Returns
        -------
        np.ndarray
            Angular velocity vector of the body  :math:`\boldsymbol{\omega}^{(B)}` expressed in the body-fixed frame :math:`B`





    )";




    } else if(name == "VehicleSystems") {
         return R"(

        Object used to store physical (hardware) properties of a vehicle.       






     )";




    } else if(name == "VehicleSystems.set_control_surface_deflection" && variant==0) {
            return R"(

        Function to set the current deflection of an aerodynamic control surface.


        Function to set the current deflection of an aerodynamic control surface,
        identified by its name. To set the control surface deflection, the control 
        surface has to exist. A control surface is created whenever control surfaces are
        defined in a body's aerodynamic coefficient interface.


        Parameters
        ----------
        control_surface_id : str
            The identified (name) of the given control surface

        deflection_angle : float
            The deflection (in radians) that the control surface is to be set to. This will 
            typically influence the aerodynamic coefficients of the vehicle





    )";


    } else if(name == "VehicleSystems.set_control_surface_deflection" && variant==1) {
            return R"(

        Function to retrieve the current deflection of an aerodynamic control surface.


        Function to retrieve the current deflection of an aerodynamic control surface,
        identified by its name. To extract the control surface deflection, the control 
        surface has to exist. A control surface is created whenever control surfaces are
        defined in a body's aerodynamic coefficient interface.


        Parameters
        ----------
        control_surface_id : str
            The identified (name) of the given control surface           

        Returns
        -------
        float
            Current deflection (in radians) that the control surface





    )";



    } else if(name == "VehicleSystems.set_control_surface_deflection" && variant==0) {
            return R"(

        Function to set the current deflection of an aerodynamic control surface.


        Function to set the current deflection of an aerodynamic control surface,
        identified by its name. To set the control surface deflection, the control 
        surface has to exist. A control surface is created whenever control surfaces are
        defined in a body's aerodynamic coefficient interface.


        Parameters
        ----------
        control_surface_id : str
            The identified (name) of the given control surface

        deflection_angle : float
            The deflection (in radians) that the control surface is to be set to. This will 
            typically influence the aerodynamic coefficients of the vehicle





    )";


    } else if(name == "VehicleSystems.set_control_surface_deflection" && variant==1) {
            return R"(

        Function to retrieve the current deflection of an aerodynamic control surface.


        Function to retrieve the current deflection of an aerodynamic control surface,
        identified by its name. To extract the control surface deflection, the control 
        surface has to exist. A control surface is created whenever control surfaces are
        defined in a body's aerodynamic coefficient interface.


        Parameters
        ----------
        control_surface_id : str
            The identified (name) of the given control surface           

        Returns
        -------
        float
            Current deflection (in radians) that the control surface





    )";



    } else if(name == "VehicleSystems.get_engine_model" && variant==0) {
            return R"(

        Function to retrieve an engine model from the vehicle



        Parameters
        ----------
        engine_name : str
            The identifier for the engine model that is to be retrieved

        Returns
        -------
        EngineModel
            Model for the engine that is requested





    )";




    } else if(name == "Ephemeris") {
         return R"(

        Object that computes the state of a body as a function of time


        Object (typically stored inside a :class:`~Body` object) that computes the state of a body as a function of time,
        both outside of a propagation, and during a propagation if the given body's translational state is not propagated.
        Note that this object computes the state w.r.t. its own origin (defined by ``frame_origin``), which need not be the same as the global frame origin
        of the environment.





     )";


    } else if(name == "Ephemeris.frame_origin") {
         return R"(

        **read-only**

        Name of the reference body/point w.r.t. which this object provides its states


        :type: str
     )";


    } else if(name == "Ephemeris.frame_orientation") {
         return R"(

        **read-only**

        Name of the frame orientation w.r.t which this object provides its states
              


        :type: str
     )";




    } else if(name == "Ephemeris.cartesian_state" && variant==0) {
            return R"(

        This function returns the Cartesian state (position and velocity) at the given time, w.r.t. the ``frame_origin``.


        Parameters
        ----------
        current_time : float
            Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.

        Returns
        -------
        np.ndarray
            Requested Cartesian state





    )";



    } else if(name == "Ephemeris.cartesian_position" && variant==0) {
            return R"(

        As ``cartesian_state``, but only the three position components


        Parameters
        ----------
        current_time : float
            Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.

        Returns
        -------
        np.ndarray
            Requested Cartesian position





    )";



    } else if(name == "Ephemeris.cartesian_velocity" && variant==0) {
            return R"(

        As ``cartesian_state``, but only the three velocity components


        Parameters
        ----------
        current_time : float
            Time (in seconds since J2000 in TDB time scale) at which the state is to be computed.

        Returns
        -------
        np.ndarray
            Requested Cartesian velocity





    )";




    } else if(name == "EngineModel") {
         return R"(

        Object used to store properties of an engine, and compute the thrust magnitude and body-fixed thrust direction     






     )";


    } else if(name == "EngineModel.thrust_magnitude_calculator") {
         return R"(

        **read-only**

        Object used to manage and compute the magnitude of the thrust provided by the engine


        :type: ThrustMagnitudeWrapper
     )";





    } else if(name == "Body") {
         return R"(

        Object that stores the environment properties and current state of
        a single body.


        Object that stores the environment properties and current state
        of a single celestial body (natural or artificial). Each separate
        environment model (gravity field, ephemeris, etc.) is stored as a
        member object in this class. During each time step, the Body gets
        updated to the current time/propagated state, and the current
        properties, in as much as they are time-dependent, can be
        extracted from this object





     )";


    } else if(name == "Body.state") {
         return R"(

        **read-only**

        The translational state of the Body, as set during the current
        step of the numerical propagation. The translational state
        stored here is always in Cartesian elements, w.r.t. the global
        frame origin, with axes along the global frame orientation. If
        the body's translational state is numerically propagated, this
        property gets extracted from the propagated state vector. If it
        is not propagated, the state is extracted from this body's
        ephemeris. In both cases, any required state transformations
        are automatically applied. Note that this function  is *only*
        valid during the numerical propagation if any aspects of the
        dynamics or dependent variables require the body's state.


        :type: numpy.ndarray
     )";


    } else if(name == "Body.position") {
         return R"(

        **read-only**

        The translational position of the Body, as set during the
        current step of the numerical propagation
        (see :py:attr:`~state`).


        :type: numpy.ndarray
     )";


    } else if(name == "Body.velocity") {
         return R"(

        **read-only**

        The translational velocity of the Body, as set during the
        current step of the numerical propagation
        (see :py:attr:`~state`).


        :type: numpy.ndarray
     )";


    } else if(name == "Body.inertial_to_body_fixed_frame") {
         return R"(

        **read-only**

        The rotation from inertial frame (with global frame
        orientation) to this Body's body-fixed frame. The rotation is
        always returned here as a rotation matrix.  If the body's
        rotational state is numerically propagated, this property gets
        extracted from the propagated state vector. If it is not
        propagated, the state is extracted from this body's rotational
        ephemeris.

        .. note:: This function is **only** valid during the
                  numerical propagation if any aspects of the dynamics
                  or dependent variables require the body's rotational
                  state.


        :type: numpy.ndarray
     )";


    } else if(name == "Body.body_fixed_to_inertial_frame") {
         return R"(

        **read-only**

        The rotation from this Body's body-fixed frame to inertial
        frame (see :py:attr:`~inertial_to_body_fixed_frame`).


        :type: numpy.ndarray
     )";


    } else if(name == "Body.inertial_to_body_fixed_frame_derivative") {
         return R"(

        **read-only**

        Time derivative of rotation matrix from inertial frame to this
        Body's body-fixed frame
        (see :py:attr:`~inertial_to_body_fixed_frame`).


        :type: numpy.ndarray
     )";


    } else if(name == "Body.body_fixed_to_inertial_frame_derivative") {
         return R"(

        **read-only**

        Time derivative of rotation matrix from this Body's body-fixed
        frame to inertial frame
        (see :py:attr:`~inertial_to_body_fixed_frame`).


        :type: numpy.ndarray
     )";


    } else if(name == "Body.inertial_angular_velocity") {
         return R"(

        **read-only**

        Angular velocity vector of the body, expressed in inertial
        frame (see :py:attr:`~inertial_to_body_fixed_frame`).


        :type: numpy.ndarray
     )";


    } else if(name == "Body.body_fixed_angular_velocity") {
         return R"(

        **read-only**

        Angular velocity vector of the body, expressed in body-fixed
        frame (see :py:attr:`~inertial_to_body_fixed_frame`).


        :type: numpy.ndarray
     )";


    } else if(name == "Body.mass") {
         return R"(

        The current mass :math:`m` of the vehicle, as used in the calculation of
        non-conservative acceleration. This attribute is a shorthand for accessing the
        mass as computed/stored in the ``rigid_body_properties`` attribute. For certain
        types of rigid-body properties, this attribute cannot be used to (re)set the current
        mass. If the body has no ``rigid_body_properties``, and this function is used to
        set a mass, a new object is automatically created, with settings analogous to the
        the :func:`~tudatpy.numerical_simulation.environment_setup.rigid_body.constant_rigid_body_properties` setting.

        Unlike the attributes containing the state, orientation, angular velocity
        of the Body, this attribute may be used to retrieve the state during the
        propagation *and* to define the mass of a vehicle.


        :type: float
     )";


    } else if(name == "Body.inertia_tensor") {
         return R"(

        The current inertia tensor :math:`\mathbf{I}` of the vehicle, as used in the calculation of
        (for instance) the reponse to torques. This attribute is a shorthand for accessing the
        inertia tensor as computed/stored in the ``rigid_body_properties`` attribute. For certain
        types of rigid-body properties, this attribute cannot be used to (re)set the current
        mass. 

        Unlike the attributes containing the state, orientation, angular velocity
        of the Body, this attribute may be used to retrieve the state during the
        propagation *and* to define the mass of a vehicle.


        :type: np.ndarray
     )";


    } else if(name == "Body.ephemeris") {
         return R"(

        Ephemeris model of this body, used to calculate its current state as a function of time.
        Depending on the selected type of model, the type of this attribute
        is of type Ephemeris, or a derived class thereof.


        :type: Ephemeris
     )";


    } else if(name == "Body.flight_conditions") {
         return R"(

        Object used to calculated and store the current flight conditions of a vehicle (altitude, latitude, longitude,
        flight-path angle, etc.) w.r.t. a central body. In case the central body contains an atmosphere, this object
        also stores current local density, Mach number, etc. This object is typically used for aerodynamic accelerations,
        guidance models or other central-body-related custom models. 


        :type: FlightConditions
     )";


    } else if(name == "Body.atmosphere_model") {
         return R"(

        Atmosphere model of this body, used to calculate density, temperature, etc. at a given
        state/time. Depending on the selected type of model, the type of this attribute
        is of type AtmosphereModel, or a derived class thereof.


        :type: AtmosphereModel
     )";


    } else if(name == "Body.shape_model") {
         return R"(

        Shape model of this body, used to define the exterior shape of the body, for instance for
        the calculation of vehicle's altitude. Depending on the selected type of model, the type of this attribute
        is of type BodyShapeModel, or a derived class thereof.


        :type: BodyShapeModel
     )";


    } else if(name == "Body.gravity_field_model") {
         return R"(

        Gravity field model of this body, used to define the exterior gravitational potential, and
        its gradient(s). Depending on the selected type of model, the type of this attribute
        is of type GravityFieldModel, or a derived class thereof.


        :type: GravityFieldModel
     )";


    } else if(name == "Body.aerodynamic_coefficient_interface") {
         return R"(

        Object defining the aerodynamic coefficients of a vehicle (force-only, or force and moment)
        as a function of any number of independent variables. Depending on the selected type of model, the type of this attribute
        is of type AerodynamicCoefficientInterface, or a derived class thereof.


        :type: AerodynamicCoefficientInterface
     )";


    } else if(name == "Body.rotation_model") {
         return R"(

        Object defining the orientation of a body, used to calculate the rotation to/from a body-fixed
        frame (and its derivate). Depending on the selected type of model, the type of this attribute
        is of type RotationalEphemeris, or a derived class thereof.


        :type: RotationalEphemeris
     )";


    } else if(name == "Body.flight_conditions") {
         return R"(

        Object used to calculated and store the current flight conditions of a vehicle (altitude, latitude, longitude,
        flight-path angle, etc.) w.r.t. a central body. In case the central body contains an atmosphere, this object
        also stores current local density, Mach number, etc. This object is typically used for aerodynamic accelerations,
        guidance models or other central-body-related custom models. 


        :type: FlightConditions
     )";


    } else if(name == "Body.system_models") {
         return R"(

        Object used to store physical (hardware) properties of a vehicle, such as engines, control surfaces, etc. This
        object is typically created automatically whenever such a hardware model needs to be assigned to a vehicle.


        :type: VehicleSystems
     )";


    } else if(name == "Body.ground_station_list") {
         return R"(

        Dictionary of all ground stations that exist in the body, with dictionary key being the name of the station,
        and the ground station object the key of the dictionary.


        :type: dict[str,GroundStation]
     )";


    } else if(name == "Body.gravitational_parameter") {
         return R"(

        **read-only**

        Attribute of convenience, equivalent to ``.gravity_field_model.gravitational_parameter``


        :type: float
     )";




    } else if(name == "Body.get_ground_station" && variant==0) {
            return R"(

        This function extracts a ground station object from the body.

        This function extracts a ground station object, for a station of a given name, from the body.
        If no station of this name exists, an exception is thrown


        Parameters
        ----------
        station_name : str
            Name of the ground station that is to be retrieved.

        Returns
        -------
        GroundStation
            Ground station object of the station of requested name





    )";




    } else if(name == "SystemOfBodies") {
         return R"(

        Object that contains a set of Body objects and associated frame
        information.


        Object that contains a set of Body objects and associated frame
        information. This object stored the entire environment for a
        typical Tudat numerical simulation, and is fundamental for the
        overall Tudat architecture.





     )";




    } else if(name == "SystemOfBodies.get" && variant==0) {
            return R"(

        This function extracts a single Body object from the SystemOfBodies.


        Parameters
        ----------
        body_name : str
            Name of the Body that is to be retrieved.

        Returns
        -------
        Body
            Body object of the requested name





    )";



    } else if(name == "SystemOfBodies.get_body" && variant==0) {
            return R"(

        Deprecated version of :py:func:`~get`





    )";



    } else if(name == "SystemOfBodies.create_empty_body" && variant==0) {
            return R"(

        This function creates a new empty body.

        This function creates a new empty body, and adds it to the
        :py:class:`~SystemOfBodies`. Since the body is empty, it will
        not have any environment models defined. These must all be
        added manually by a user.


        Parameters
        ----------
        body_name : string
            Name of the Body that is to be added

        process_body : bool, default=True
            Variable that defines whether this new Body will have its
            global frame origin/orientation set to conform to rest of
            the environment.

            .. warning:: Only in very rare cases should
                         this variable be anything other than ``True``.
                         Users are recommended to keep this default value
                         intact.





        Examples
        --------

        This function is often used early on in the environment
        creation segment of a simulation, following the creation of
        a :py:class:`~SystemOfBodies` from the default settings
        for celestial bodies.

        .. code-block:: python
           :emphasize-lines: 18

           # Define string names for bodies to be created from default.
           bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

           # Use "Earth"/"J2000" as global frame origin and orientation.
           global_frame_origin = "Earth"
           global_frame_orientation = "J2000"

           # Create default body settings, usually from `spice`.
           body_settings = environment_setup.get_default_body_settings(
               bodies_to_create,
               global_frame_origin,
               global_frame_orientation)

           # Create system of selected celestial bodies
           bodies = environment_setup.create_system_of_bodies(body_settings)

           # Create vehicle objects.
           bodies.create_empty_body("Delfi-C3")

    )";



    } else if(name == "SystemOfBodies.add_body" && variant==0) {
            return R"(

        This function adds an existing body, which the user has
        separately created, to the :py:class:`~SystemOfBodies`.



        Parameters
        ----------
        body_to_add : Body
            Body object that is to be added.

        body_name : numpy.ndarray
            Name of the Body that is to be added.

        process_body : bool, default=True
            Variable that defines whether this new Body will have its
            global frame origin/orientation set to conform to rest of
            the environment.

            .. warning:: Only in very rare cases should this variable be
                         anything other than ``True``. Users are
                         recommended to keep this default value intact.





    )";



    } else if(name == "SystemOfBodies.remove_body" && variant==0) {
            return R"(

        This function removes an existing body from the
        :py:class:`~SystemOfBodies`.



        .. warning:: This function does *not* necessarily delete the
                     Body object, it only removes it from this object.
                     If any existing models in the simulation refer to
                     this Body, it will persist in memory.


        Parameters
        ----------
        body_name : numpy.ndarray
            Name of the Body that is to be removed.





    )";





    } else if(name == "save_vehicle_mesh_to_file" && variant==0) {
        return R"(
        
Function to save the mesh used for a hypersonic local inclination analysis to a file.

Function to save the mesh used for a hypersonic local inclination analysis to a file. This function saves 
two files to the specified directory, with filenames: "ShapeFile.dat" and "SurfaceNormalFile.dat", where these
files names may be prefixed by an optional string (see below). The first of these files contains four columns defining
the surface points that define mesh, with Column 0: point index; Column 1: x-position of point; Column 1: y-position of point; 
Column 2: z-position of point. The second file contains four columns with Column 0: point index; Column 1: x-component of surface normal; 
Column 1: y-position of surface normal; Column 2: z-position of surface normal.


Parameters
----------
local_inclination_analysis_object : HypersonicLocalInclinationAnalysis
    Object used to calculate the aerodynamics of the vehicle

output_directory : str
    Directory to which the files are to be saved

output_file_prefix : str, default=''
    Optional prefix of output file names






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace propagation {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "PropagationTerminationReason") {
         return R"(

        Enumeration of types of termination of propagation.





     )";


    } else if(name == "PropagationTerminationReason.propagation_never_run") {
         return R"(
     )";


    } else if(name == "PropagationTerminationReason.unknown_reason") {
         return R"(
     )";


    } else if(name == "PropagationTerminationReason.termination_condition_reached") {
         return R"(
     )";


    } else if(name == "PropagationTerminationReason.runtime_error_caught_in_propagation") {
         return R"(
     )";


    } else if(name == "PropagationTerminationReason.nan_or_inf_detected_in_state") {
         return R"(
     )";




    } else if(name == "SimulationResults") {
         return R"(

        Base class for objects that store all results of a numerical propagation.

        Base class for objects that store all results of a numerical propagation. Derived class are implemented for single-, multi- and hybrid-arc propagation of botj dynamics and variational equations





     )";





    } else if(name == "SingleArcSimulationResults") {
         return R"(

        Class that stores all the results (including logging data) of a single-arc propagation






     )";


    } else if(name == "SingleArcSimulationResults.state_history") {
         return R"(

        **read-only**

        Numerical solution of the equations of motion as key-value pairs. The key denotes the epoch. The value contains the
        numerically calculated state at this epoch. For this function, the states are always converted to so-called
        'processed' formulations (e.g. Cartesian states for translational dynamics), see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/processed_propagated_elements.html>`_
        for details. For the history of the states that were actually propagated, use the ``unprocessed_state_history``.

        .. note:: The propagated state at each epoch contains the state types in the following order: Translational ( **T** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
                  When propagating two bodies, an example of what the output state would look like is for instance:
                  [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ] The specifics can be retrieved using the :attr:`state_ids` attribute of this class

        .. note:: For propagation of translational dynamics using cowell
                  propagator, the conventional and propagated
                  coordinates are identical.


        :type: dict[float, numpy.ndarray]
     )";


    } else if(name == "SingleArcSimulationResults.unprocessed_state_history") {
         return R"(

        **read-only**

        Numerical solution of the equations of motion as key-value pairs, without any processing applied. The key denotes the epoch. The value contains the
        numerically calculated state at this epoch. This attribute contains the states of the propagated bodies expressed in the
        "raw" form in which the propagation took place. For instance, when using a Gauss-Kepler propagation scheme, this
        attribute will contain the numerically propagated Keplerian elements at each time epoch


        :type: dict[float, numpy.ndarray]
     )";


    } else if(name == "SingleArcSimulationResults.dependent_variable_history") {
         return R"(

        **read-only**

        Dependent variables computed during the propagation as key-value pairs.
        The vector of all dependent variables concatenated into a single vector as value, with the epoch as key.
        They order of the concatenated dependent variables in a single value is provided by the ``dependent_variable_ids`` attribute of this object. 


        :type: dict[float, numpy.ndarray]
     )";


    } else if(name == "SingleArcSimulationResults.cumulative_computation_time_history") {
         return R"(

        **read-only**

        History of cumulative computation time in seconds needed during the propagation as key-value
        pairs. At each epoch (key) the computation time (value) in seconds is the total computation time
        used up to and including that time step. This includes the total time up to and including the current time step,
        since the beginning of the (single-arc) propagation.


        :type: dict[float, float]
     )";


    } else if(name == "SingleArcSimulationResults.cumulative_number_of_function_evaluations") {
         return R"(

        **read-only**

        This function returns the history of cumulative number of function evaluations taken during the propagation as key-value
        pairs. At each epoch (key), the total number of computed function evaluations (value) are given up to and including that time step.
        This includes all function evaluations up to and including the current time step, since the beginning of the (single-arc) propagation.


        :type: dict[float, int]
     )";


    } else if(name == "SingleArcSimulationResults.termination_details") {
         return R"(

        **read-only**

        Object describing the details of the event that triggered the termination of the last propagation.


        :type: PropagationTerminationDetails
     )";


    } else if(name == "SingleArcSimulationResults.dependent_variable_ids") {
         return R"(

        **read-only**

        Key-value container with the starting entry of the dependent variables saved (key), along with associated ID (value).


        :type: dict[[int,int], str]
     )";


    } else if(name == "SingleArcSimulationResults.state_ids") {
         return R"(

        **read-only**

        Key-value container with the starting entry of the states (key), along with associated ID (value).


        :type: dict[[int,int] str]
     )";


    } else if(name == "SingleArcSimulationResults.integration_completed_successfully") {
         return R"(

        **read-only**

        Boolean defining whether the last propagation was finished
        successfully, as defined by the termination conditions, or if
        it was terminated prematurely (for instance due to an
        exception, or an Inf/NaN state entry being detected).


        :type: bool
     )";





    } else if(name == "PropagationTerminationDetails") {
         return R"(

        Class that provides information on the reason for the
        termination of the propagation.






     )";


    } else if(name == "PropagationTerminationDetails.termination_reason") {
         return R"(

        Enum defining the reason the propagation was terminated


        :type: PropagationTerminationReason
     )";


    } else if(name == "PropagationTerminationDetails.terminated_on_exact_condition") {
         return R"(

        Boolean defining whether the propagation was terminated on an *exact* final condition,
        or once the propagation went *past* the determined final condition. The choice of behaviour is
        defined by the termination settings provided as input to the Simulator object. This variable only
        has a meaningful definition if the ``termination_reason`` has value ``termination_condition_reached``


        :type: bool
     )";





    } else if(name == "PropagationTerminationDetailsFromHybridCondition") {
         return R"(

        Class that provides information on the reason for the termination of the propagation, for hybrid termination conditions


        Derived class from :class:`PropagationTerminationDetails` that provides information on the reason for the termination of the propagation,
        for the case of hybrid termination conditions (defined using the :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.hybrid_termination`)
        function





     )";


    } else if(name == "PropagationTerminationDetailsFromHybridCondition.was_condition_met_when_stopping") {
         return R"(

        List of booleans defining, per entry in ``termination_settings`` when calling :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.hybrid_termination`,
        whether the corresponding entry of the hybrid termination settings was met or not


        :type: list[bool]
     )";





    } else if(name == "RotationalProperModeDampingResults") {
         return R"(

        Object that stores the results of the algorithm to damp the proper mode of rotational dynamics for an initial state,
        as computed by the :func:`~get_damped_proper_mode_initial_rotational_state` function






     )";


    } else if(name == "RotationalProperModeDampingResults.damped_initial_state") {
         return R"(

        Initital state produced by the damping algorithm, for which the signature of the proper mode should be
        removed (or at least, substantially reduced). Note that this initial state corresponds to the *full* state vector
        that is provided to the ``get_damped_proper_mode_initial_rotational_state`` function (e.g. is size 7
        for rotational dynamics of a single body, size 13 for coupled orbital-rotational dynamics of a single body, etc.)


        :type: np.ndarray
     )";


    } else if(name == "RotationalProperModeDampingResults.forward_backward_states") {
         return R"(

        Data structure that contains the full state histories used by the damping algorithm. The contents are are as follows:

        * The :math:`i^{th}` entry of the list corresponds to the :math:`i^{th}` iteration of the forward-backward propagation
        * Each tuple in the list contains two dictionaries, the first one corresponding to the forward propagation results, the seconds one to the backward propagation results


        :type: list[tuple[dict[float,np.ndarray],dict[float,np.ndarray]]]
     )";


    } else if(name == "RotationalProperModeDampingResults.forward_backward_dependent_variables") {
         return R"(

        As ``forward_backward_states``, but for the dependent variables.


        :type: list[tuple[dict[float,np.ndarray],dict[float,np.ndarray]]]
     )";






    } else if(name == "get_state_of_bodies" && variant==0) {
        return R"(
        
Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time.

Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time. This function
is typically used to extract an initial state for a propagation of a set of bodies, for which the initial state is extracted from the
existing ephemerides of the bodies.


Parameters
----------
bodies_to_propagate : list[str]
    List of names of bodies for which the state is to be extracted
central_bodies : list[str]
    List of central bodies, w.r.t. which the states are to be computed (in the same order as ``bodies_to_propagate``)
bodies_to_propagate : SystemOfBodies
    System of bodies that define the environment
initial_time : float
    Time at which the states are to be extracted from the environment
Returns
-------
numpy.ndarray
    Vector of size :math:`6\times N`, with the translational states of each entry of body from 
    ``bodies_to_propagate`` w.r.t. the corresponding central body in ``central_bodies``. 







    )";



    } else if(name == "get_damped_proper_mode_initial_rotational_state" && variant==0) {
        return R"(
        
Function to compute an initial rotational state for which the proper mode of rotation is damped.

Function to compute an initial rotational state for which the proper mode of rotation is damped, using the algorithm 
used by Rambaux et al. (2010) to compute an initial rotational state for Phobos. This algorithm propagates the
dynamics of the system a number of times, with the settings specified by the user and a specific modification to
damp the proper mode. Since a number of propagations are performed by this function, it may take some time to run.
Specifically, the algorithm works as follows:

* Introduce a damping torque (see below) to damp the proper mode, with damping time :math:`\tau_{d}`
* Propagate the dynamics forward in time for a period of :math:`10\tau_{d}`
* Remove the virtual torque, and propagate the dynamics back to the initial time :math:`t_{0}`
* Repeat the above for the list of damping times provided by the user

The state after the final backwards propagation to :math:`t_{0}` is provided as output by this function, to be
used as damped initial state. The output from this function also provides the user access to the full state history
and dependent variable history of the forward and backward propagations, to allow a user to track and validate
the pgress of the algorithm.

The damping torque :math:`\Gamma` is defined as follows:

.. math::
   \boldsymbol{\Gamma}= -\frac{1}{\tau_{d}}\mathbf{I}\begin{pmatrix}\omega_{x}\\ \omega_{y}\\ \omega_{x}-\omega_{p} \end{pmatrix}

where :math:\mathbf{I}` is the body's inertia tensor (in its body-fixed frame), :math:`\tau_{d}` the damping time of the
current propagation, and :math:`\omega_{x}, \omega_{y}, \omega_{z}` the body's current rotation about its
body-fixed, x-, y- and z-axes, respectively. The damping torque is implemented to damp out all rotations along 
the body-fixed x- and y-axes, and any deviations from constant rotation with frequency :\omega_{p}: about the body-fixed z-axis. 

.. note:: The mean rotation rate of the body :math:`\omega_{p}` is a user-defined input, and must be tuned to the dynamics of the system.


Parameters
----------
bodies : SystemOfBodies
    Set of body objects that defines the environment
propagator_settings : SingleArcPropagatorSettings
    Propagator settings for the dynamics of which the initial rotational state is to be damped. These propagator
    settings must be for rotational dynamics only, or for multi-type rotational dynamics that contains rotational
    dynamics for a single body (e.g. translational-rotational dynamics for a single body)

body_mean_rotational_rate : float
    Mean rotational rate :math:`\omega_{p}` to which the damping algorithm will force the body-fixed rotation about its z-axis.
dissipation_times : list[ float ]
    List of damping times :math:`\tau_{d}` for which the algorithm is to be run. Note that this list should be organized in ascending order for the algorithm to perform properly
propagate_undamped : bool, default = True
    Boolean defining whether the first forward/backward propagation performed by the damping algorithm has damping turned off (damping turned off if True, damping turned on if False).
    Propagating without any damping before starting the damping algorithm is useful for verification purposes, but not required for the algorithm itself.

Returns
-------
DampedInitialRotationalStateResults
    Object that contains the results of the damping algorithm (final damped rotational state, and forward/backward propagation results).






    )";



    } else if(name == "combine_initial_states" && variant==0) {
        return R"(
        
Function to retrieve the initial state for a list of propagator settings.

Function to retrieve the initial state for a list of propagator settings. This way, the initial state for
different quantities to be propagated (e.g., translational state, rotational state, mass) are retrieved and
organized in a single container.


Parameters
----------
propagator_settings_per_type : dict
    Propagator settings where the type of propagation is reported as key and the respective list of propagator settings as value.
Returns
-------
numpy.ndarray
    Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the vector of SingleArcPropagatorSettings of given type.






    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace trajectory_design {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else {
        return "No documentation found.";
    }

}


    
namespace transfer_trajectory {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "TransferLegTypes") {
         return R"(

        Enumeration of available leg types.





     )";


    } else if(name == "TransferLegTypes.unpowered_unperturbed_leg_type") {
         return R"(
     )";


    } else if(name == "TransferLegTypes.dsm_position_based_leg_type") {
         return R"(
     )";


    } else if(name == "TransferLegTypes.dsm_velocity_based_leg_type") {
         return R"(
     )";


    } else if(name == "TransferLegTypes.spherical_shaping_low_thrust_leg") {
         return R"(
     )";


    } else if(name == "TransferLegTypes.hodographic_low_thrust_leg") {
         return R"(
     )";




    } else if(name == "TransferLeg") {
         return R"(

        Base class for defining a transfer leg.

        Functional (base) class for transfer legs, requiring the leg type, departure body ephemeris and arrival body ephemeris.
        Transfer node classes requiring additional information must be created using an object derived from this class.





     )";





    } else if(name == "SphericalShapingLeg") {
         return R"(

        Class for defining low-thrust spherical-shaping leg.

        `TransferLeg` derived class for defining a low-thrust leg described using spherical shaping [3]_.





     )";





    } else if(name == "HodographicShapingLeg") {
         return R"(

        Class for defining low-thrust hodographic-shaping leg.

        `TransferLeg` derived class for defining a low-thrust leg described using hodographic shaping [4]_.





     )";





    } else if(name == "TransferNodeSettings") {
         return R"(

        Base class for providing settings for transfer nodes.

        Functional (base) class for settings of transfer nodes that require no information in addition to their type.
        Transfer node classes requiring additional information must be created using an object derived from this class.





     )";





    } else if(name == "SwingbyNodeSettings") {
         return R"(

        Class for defining settings of swingby node.

        `TransferNodeSettings` derived class for providing settings for swingby nodes, which consist of the minimum periapsis
        radius.





     )";





    } else if(name == "EscapeAndDepartureNodeSettings") {
         return R"(

        Class for defining settings of escape and departure node.

        `TransferNodeSettings` derived class for providing settings for escape and departure nodes, which consist of the
        departure semi-major axis and eccentricity.





     )";





    } else if(name == "CaptureAndInsertionNodeSettings") {
         return R"(

        Class for defining settings of capture and insertion node.

        `TransferNodeSettings` derived class for providing settings for capture and insertion nodes, which consist of the
        capture semi-major axis and eccentricity.





     )";





    } else if(name == "TransferLegSettings") {
         return R"(

        Base class for providing settings for transfer legs.

        Functional (base) class for settings of transfer legs that require no information in addition to their type.





     )";





    } else if(name == "TransferTrajectory") {
         return R"(

        Class defining a transfer trajectory constituted by transfer legs and nodes.

        Class defining a transfer trajectory constituted by transfer legs and nodes. The object is tipically created using the `create_transfer_trajectory` function.




     )";


    } else if(name == "TransferTrajectory.delta_v") {
         return R"(

        **read-only**

        Total Delta V used in the transfer trajectory.

        :type: float
     )";


    } else if(name == "TransferTrajectory.time_of_flight") {
         return R"(

        **read-only**

        Total time of flight of the transfer trajectory.

        :type: float
     )";


    } else if(name == "TransferTrajectory.delta_v_per_node") {
         return R"(

        **read-only**

        List of the Delta V applied in each node.

        :type: list[float]
     )";


    } else if(name == "TransferTrajectory.delta_v_per_leg") {
         return R"(

        **read-only**

        List of the Delta V applied in each leg.

        :type: list[float]
     )";


    } else if(name == "TransferTrajectory.number_of_nodes") {
         return R"(

        **read-only**

        Number of nodes in the transfer trajectory.

        :type: float
     )";


    } else if(name == "TransferTrajectory.number_of_legs") {
         return R"(

        **read-only**

        Number of legs in the transfer trajectory.

        :type: float
     )";




    } else if(name == "TransferTrajectory.evaluate" && variant==0) {
            return R"(

        Evaluate transfer trajectory.

        Function to evaluate the transfer trajectory, which consists of computing the transfer Delta V and time of
        flight, for the specified set of parameters.


        Parameters
        ----------
        node_times : list[float]
            List of the time at each node.
        leg_parameters : list[list[float]]
            List of lists with the parameters characterizing each leg. Each inner list corresponds to the
            parameters of one leg; if a leg does not require any parameter, its list can contain any value(s),
            therefore it is recommended to leave it empty.

        node_parameters : list[list[float]]
            List of lists with the parameters characterizing each node. Each inner list corresponds to the
            parameters of one node; if a node does not require any parameter, its list can contain any value(s),
            therefore it is recommended to leave it empty.

        Returns
        -------
        None
            None





    )";



    } else if(name == "TransferTrajectory.single_node_delta_v" && variant==0) {
            return R"(

        Retrieves the Delta V applied in the specified node.


        Parameters
        ----------
        node_index : int
            Index of the node for which the Delta V to be is retrieved.
        Returns
        -------
        float
            Delta V for the specified node.





    )";



    } else if(name == "TransferTrajectory.single_leg_delta_v" && variant==0) {
            return R"(

        Retrieves the Delta V applied in the specified leg.


        Parameters
        ----------
        leg_index : int
            Index of the leg for which the Delta V is to be retrieved.
        Returns
        -------
        float
            Delta V for the specified leg.





    )";



    } else if(name == "TransferTrajectory.states_along_trajectory" && variant==0) {
            return R"(

        Returns the state history throughout the trajectory.

        Function that returns the state history throughout the trajectory, using the same number of data points in
        each leg. For each leg, the retrieved states are equally spaced in time.


        Parameters
        ----------
        number_of_data_points_per_leg : int
            Number of data points used to describe each leg.
        Returns
        -------
        tuple[numpy.ndarray,numpy.ndarray]
            Tuple of (state history, time history).





    )";



    } else if(name == "TransferTrajectory.inertial_thrust_accelerations_along_trajectory" && variant==0) {
            return R"(

        Returns the inertial thrust acceleration history throughout the trajectory.

        Function that returns the inertial thrust acceleration history throughout the trajectory, using the same number of data points in
        each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
        For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.


        Parameters
        ----------
        number_of_data_points_per_leg : int
            Number of data points used to describe each leg.
        Returns
        -------
        tuple[numpy.ndarray,numpy.ndarray]
            Tuple of (state history, time history).





    )";



    } else if(name == "TransferTrajectory.rsw_thrust_accelerations_along_trajectory" && variant==0) {
            return R"(

        Returns the thrust acceleration history in the RSW frame throughout the trajectory.

        Function that returns the thrust acceleration history in the RSW frame throughout the trajectory, using the same number of data points in
        each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
        For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.


        Parameters
        ----------
        number_of_data_points_per_leg : int
            Number of data points used to describe each leg.
        Returns
        -------
        tuple[numpy.ndarray,numpy.ndarray]
            Tuple of (state history, time history).





    )";



    } else if(name == "TransferTrajectory.tnw_thrust_accelerations_along_trajectory" && variant==0) {
            return R"(

        Returns the thrust acceleration history in the TNW frame throughout the trajectory.

        Function that returns the thrust acceleration history in the TNW frame throughout the trajectory, using the same number of data points in
        each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
        For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.


        Parameters
        ----------
        number_of_data_points_per_leg : int
            Number of data points used to describe each leg.
        Returns
        -------
        tuple[numpy.ndarray,numpy.ndarray]
            Tuple of (state history, time history).





    )";





    } else if(name == "mga_settings_unpowered_unperturbed_legs" && variant==0) {
        return R"(
        
Function to get the legs and nodes settings of a transfer with just upowered legs.


Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
one initial node (departure or swingby), unpowered transfer leg(s) connected by swingby nodes, one final
(capture or swingby) node.
If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
node.


Parameters
----------
body_order : list[str]
    List of bodies to visit, including departure body, swingby bodies and arrival body.
departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
    node as a swingby node (instead of a departure node).

arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
    node as a swingby node (instead of a capture node).

minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
    Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
    value. Default values from Izzo [1]_.

Returns
-------
tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
    Tuple specifying the settings of each transfer leg and node.






    )";



    } else if(name == "mga_settings_dsm_velocity_based_legs" && variant==0) {
        return R"(
        
Function to get the legs and nodes settings of a transfer constituted by legs with 1 impulsive deep space maneuver (DSM)
described using the velocity formulation.


Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
one initial node (departure or swingby), velocity-based DSM transfer leg(s) connected by swingby nodes, one final
(capture or swingby) node.
If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
node.


Parameters
----------
body_order : list[str]
    List of bodies to visit, including departure body, swingby bodies and arrival body.
departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
    node as a swingby node (instead of a departure node).

arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
    node as a swingby node (instead of a capture node).

minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
    Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
    value. Default values from Izzo [1]_.

Returns
-------
tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
    Tuple specifying the settings of each transfer leg and node.






    )";



    } else if(name == "mga_settings_dsm_position_based_legs" && variant==0) {
        return R"(
        
Function to get the legs and nodes settings of a transfer constituted by legs with 1 impulsive deep space maneuver (DSM)
described using the position formulation.


Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
one initial node (departure or swingby), position-based DSM transfer leg(s) connected by swingby nodes, one final
(capture or swingby) node.
If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
node.


Parameters
----------
body_order : list[str]
    List of bodies to visit, including departure body, swingby bodies and arrival body.
departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
    node as a swingby node (instead of a departure node).

arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
    node as a swingby node (instead of a capture node).

minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
    Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
    value. Default values from Izzo [1]_.

Returns
-------
tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
    Tuple specifying the settings of each transfer leg and node.






    )";



    } else if(name == "mga_settings_spherical_shaping_legs" && variant==0) {
        return R"(
        
Function to get the legs and nodes settings of a transfer constituted by low-thrust spherical shaping legs.


Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
one initial node (departure or swingby), spherical shaping leg(s) connected by swingby nodes, one final
(capture or swingby) node.
If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
node.


Parameters
----------
body_order : list[str]
    List of bodies to visit, including departure body, swingby bodies and arrival body.
root_finder_settings : RootFinderSettings
    Settings of the root finder used by the spherical shaping leg when computing the value of the free coefficient
    that allows meeting the desired time of flight.

departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
    node as a swingby node (instead of a departure node).

arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
    node as a swingby node (instead of a capture node).

lower_bound_free_coefficient : float, default=TUDAT_NAN
    Lower bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
    it must be specified if the selected root finder requires the definition of a lower bound.

upper_bound_free_coefficient : float, default=TUDAT_NAN
    Upper bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
    it must be specified if the selected root finder requires the definition of an upper bound.

initial_value_free_coefficient : float, default=TUDAT_NAN
    Initial guess for the free coeficient. Parameter is potentially used by the root finder:
    it must be specified if the selected root finder requires the definition of an initial guess.

minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
    Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
    value. Default values from Izzo [1]_.

Returns
-------
tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
    Tuple specifying the settings of each transfer leg and node.






    )";



    } else if(name == "mga_settings_hodographic_shaping_legs" && variant==0) {
        return R"(
        
Function to get the legs and nodes settings of a transfer constituted by low-thrust hodographic shaping legs,
with user-provided velocity shaping functions.


Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
one initial node (departure or swingby), hodographic shaping leg(s) connected by swingby nodes, one final
(capture or swingby) node.
If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
node.


Parameters
----------
body_order : list[str]
    List of bodies to visit, including departure body, swingby bodies and arrival body.
radial_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
    List with the lists of radial velocity function components used in each leg.

normal_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
    List with the lists of normal velocity function components used in each leg.

axial_velocity_function_components_per_leg : list[list[BaseFunctionHodographicShaping]]
    List with the lists of axial velocity function components used in each leg.

departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
    node as a swingby node (instead of a departure node).

arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
    node as a swingby node (instead of a capture node).

minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
    Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
    value. Default values from Izzo [1]_.

Returns
-------
tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
    Tuple specifying the settings of each transfer leg and node.






    )";



    } else if(name == "mga_settings_hodographic_shaping_legs_with_recommended_functions" && variant==0) {
        return R"(
        
Function to get the legs and nodes settings of a transfer constituted by low-thrust hodographic shaping legs,
using the recommended velocity shaping functions.


Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
one initial node (departure or swingby), hodographic shaping leg(s) connected by swingby nodes, one final
(capture or swingby) node. Hodographic shaping legs use the velocity-shaping functions recommended by Musegaas,
2012 [2]_.
If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
node.


Parameters
----------
body_order : list[str]
    List of bodies to visit, including departure body, swingby bodies and arrival body.
time_of_flight_per_leg : list[float]
    List with the time of flight of each leg.

time_of_flight_per_leg : list[float]
    List with the number of revolutions of each leg.

departure_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
    node as a swingby node (instead of a departure node).

arrival_orbit : tuple[float, float], default=(TUDAT_NAN, TUDAT_NAN)
    Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
    node as a swingby node (instead of a capture node).

minimum_pericenters : dict[str, float], default=DEFAULT_MINIMUM_PERICENTERS
    Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
    value. Default values from Izzo [1]_.

Returns
-------
tuple[ list[TransferLegSettings], list[TransferNodeSettings] ]
    Tuple specifying the settings of each transfer leg and node.






    )";



    } else if(name == "unpowered_leg" && variant==0) {
        return R"(
        
Factory function for creating the settings of an unpowered leg.


Factory function for creating the settings of an unpowered leg; the settings consist of just the leg type.
Given the departure and arrival position, this leg computes the departure and arrival velocity using a Lambert
targeter.
The calculations performed in this leg do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Musegaas, 2012 [2]_.

Returns
-------
TransferLegSettings
    Transfer leg settings.






    )";



    } else if(name == "dsm_position_based_leg" && variant==0) {
        return R"(
        
Factory function for creating the settings of a transfer leg with 1 impulsive deep space maneuver (DSM) described using
the position formulation.


Factory function for creating the settings of a transfer leg with 1 position-based DSM; the settings consist of just the leg type.
Given the departure position and the DSM position this leg uses a Lambert targeter to compute the departure
velocity and the velocity before the DSM. Given the DSM position and the arrival position, the leg also uses
a Lambert targeter to compute the velocity after the DSM and the arrival velocity. The Delta V applied in the
DSM is computed using the velocity before and after the DSM.
The calculations performed in this leg do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Musegaas, 2012 [2]_.

Returns
-------
TransferLegSettings
    Transfer leg settings.






    )";



    } else if(name == "dsm_velocity_based_leg" && variant==0) {
        return R"(
        
Factory function for creating the settings of a transfer leg with 1 impulsive deep space maneuver (DSM) described using
the velocity formulation.


Factory function for creating the settings of a transfer leg with 1 velocity-based DSM; the settings consist of just the leg type.
Given the departure position and velocity this leg "propagates" the Kepler elements until the instant of application
of the DSM (giving the position at the DSM and the velocity before the DSM). Given the position of the DSM
and the arrival position, it computes the velocity after the DSM
(which is used to compute the Delta V applied in the DSM) and the arrival velocity using a Lambert targeter.
The calculations performed in this leg do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Musegaas, 2012 [2]_.

Returns
-------
TransferLegSettings
    Transfer leg settings.






    )";



    } else if(name == "spherical_shaping_leg" && variant==0) {
        return R"(
        
Factory function for creating the settings of a low-thrust spherical shaping leg.


Factory function for creating the settings of a low-thrust spherical shaping leg; the settings consist of
variables necessary for setting up the root finder and variable to set up an interpolator.
The trajectory is determined via spherical shaping, which shapes the position and time history throughout the
transfer. The trajectory depends on a single parameter, which is selected using the root finder in order to meet
a user-specified time of flight.
The calculations performed in this leg do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Roegiers, 2014 [3]_.


Parameters
----------
root_finder_settings : RootFinderSettings
    Settings of the root finder used by the spherical shaping leg when computing the value of the free coefficient
    that allows meeting the desired time of flight.

lower_bound_free_coefficient : float, default=TUDAT_NAN
    Lower bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
    it must be specified if the selected root finder requires the definition of a lower bound.

upper_bound_free_coefficient : float, default=TUDAT_NAN
    Upper bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
    it must be specified if the selected root finder requires the definition of an upper bound.

initial_value_free_coefficient : float, default=TUDAT_NAN
    Initial guess for the free coeficient. Parameter is potentially used by the root finder:
    it must be specified if the selected root finder requires the definition of an initial guess.

time_to_azimuth_interpolator_step_size : float, default=physical_constants::JULIAN_DAY
    Time step size used as reference to define the azimuth values at which the epoch is computed, when defining an
    interpolator to convert between epoch and azimuth.

Returns
-------
TransferLegSettings
    Transfer leg settings.






    )";



    } else if(name == "hodographic_shaping_leg" && variant==0) {
        return R"(
        
Factory function for creating the settings of a low-thrust hodographic shaping leg.


Factory function for creating the settings of a low-thrust hodographic shaping leg; the settings consist of
the functions used to shape the velocity throughout the transfer.
Note that shape functions with at least 3 terms (i.e. 3 degrees of freedom) must be provided for each velocity
component (radial, normal and axial); this is required to obtain a trajectory that satisfies the boundary
conditions.
The calculations performed in this leg do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Gondelach, 2012 [4]_.


Parameters
----------
radial_velocity_function_components : list[BaseFunctionHodographicShaping]
    List with components of the radial velocity shaping function, which determine the radial velocity profile
    throughout the transfer.

normal_velocity_function_components : list[BaseFunctionHodographicShaping]
    List with components of the normal velocity shaping function, which determine the normal velocity profile
    throughout the transfer.

axial_velocity_function_components : list[BaseFunctionHodographicShaping]
    List with components of the axial velocity shaping function, which determine the axial velocity profile
    throughout the transfer.

Returns
-------
TransferLegSettings
    Transfer leg settings.






    )";



    } else if(name == "swingby_node" && variant==0) {
        return R"(
        
Factory function for creating the settings of a swingby node.

Factory function for creating the settings of a swingby node. The settings consist consist of the minimum
allowed periapsis radius.
The minimum periapsis radius can be set to infinity. In that case, the swingby does not affect the velocity of the
spacecraft (e.g. might be relevant for swingbys of small bodies).
The exact behavior of this node depends on the types of legs that precede and follow it. Given a known incoming
and unknown outgoing velocity, the node forward propagates the gravity assist, possibly with a Delta V at
the a periapsis. Given an unknown incoming and known outgoing velocity, the node backward propagates the
gravity assist, possibly with a Delta V at the a periapsis. Given known incoming and outgoing velocities,
the node computes the Delta V required to meet those.
The calculations performed in this node do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Musegaas, 2012 [2]_.


Parameters
----------
minimum_periapsis : float, default=TUDAT_NAN
    Minimum periapsis radius. The minimum periapsis only needs to be specified if the types of swignby nodes that 
    requires it is used. If that is the case and no minimum periapsis was selected an error is thrown.

Returns
-------
SwingbyNodeSettings
    Swingby node settings.






    )";



    } else if(name == "departure_node" && variant==0) {
        return R"(
        
Factory function for creating the settings of an escape or departure node.

Factory function for creating the settings of an escape or departure node. The settings consist of the
departure orbit eccentricity and semi-major axis.
Given the initial orbit and the departure velocity, the node computes the Delta V that needs to be applied
at the periapsis of the initial orbit to enter the escape trajectory.
The calculations performed in this node do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Musegaas, 2012 [2]_.


Parameters
----------
departure_semi_major_axis : float
    Departure orbit semi-major axis.
departure_eccentricity : float
    Departure orbit eccentricity.
Returns
-------
EscapeAndDepartureNodeSettings
    Escape or departure node settings.






    )";



    } else if(name == "capture_node" && variant==0) {
        return R"(
        
Factory function for creating the settings of a capture or insertion node.

Factory function for creating the settings of a capture or insertion node. The settings consist of the
capture orbit eccentricity and semi-major axis.
Given the the arrival velocity and the final orbit, the node computes the Delta V that needs to be applied
at the periapsis of the final orbit to exit the capture trajectory.
The calculations performed in this node do not involve any numerical integration, and are solved by
(semi-)analytical models. For details see Musegaas, 2012 [2]_.


Parameters
----------
capture_semi_major_axis : float
    Capture orbit semi-major axis.
capture_eccentricity : float
    Capture orbit eccentricity.
Returns
-------
CaptureAndInsertionNodeSettings
    Capture or insertion node settings.






    )";



    } else if(name == "print_parameter_definitions" && variant==0) {
        return R"(
        
Prints the list of parameters required to define the transfer trajectory, according to the
specified node and leg settings.



Parameters
----------
leg_settings : list[TransferLegSettings]
    List of transfer leg settings.
node_settings : list [TransferNodeSettings]
    List of transfer node settings.
Returns
-------
None
    None






    )";



    } else if(name == "create_transfer_trajectory" && variant==0) {
        return R"(
        
Factory function for creating a transfer trajectory consisting of the specified sequence of transfer nodes and
transfer legs.


The function creates a transfer trajectory based on the provided transfer nodes settings and transfer legs
settings. The number of nodes should be equal to the number of legs plus 1.
This function creates an instance of the `TransferTrajectory` class.


Parameters
----------
bodies : SystemOfBodies
    System of bodies to be used in the transfer trajectory.
leg_settings : list[TransferLegSettings]
    List of transfer leg settings.
node_settings : list [TransferNodeSettings]
    List of transfer node settings.
node_names : list [str]
    Sequence of bodies used as transfer nodes.
central_body : str
    Central body with respect to which the two-body trajectory of the spacecraft
    is calculated.

Returns
-------
TransferTrajectory
    Transfer trajectory object.






    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace shape_based_thrust {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "BaseFunctionHodographicShaping") {
         return R"(

        Base class for defining settings of the shape functions for hodographic shaping method.

        Base class for defining settings of the shape functions for Hodograph shaping method. Objects derived
        from this class are created by calling the dedicated factory functions in this module





     )";






    } else if(name == "recommended_radial_hodograph_functions" && variant==0) {
        return R"(
        
Factory function for creating the default radial hodographic trajectory shaping functions.

Factory function for creating the default radial hodographic trajectory shaping functions. This function 
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
    List of default settings object for radial hodographic shaping






    )";



    } else if(name == "recommended_normal_hodograph_functions" && variant==0) {
        return R"(
        
Factory function for creating the default normal hodographic trajectory shaping functions.

Factory function for creating the default normal hodographic trajectory shaping functions. This function 
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
    List of default settings object for axial hodographic shaping






    )";



    } else if(name == "recommended_axial_hodograph_functions" && variant==0) {
        return R"(
        
Factory function for creating the default axial hodograph	ic trajectory shaping functions.

Factory function for creating the default axial hodographic trajectory shaping functions. This function 
(and its counterparts radial and normal components) provided three shaping functions that have been found in
literature to work well for this method. For a given time-of-flight :math:`T` and number of revolutions :math:`N`, this function returns a list of
three shaping functions:

* Cosine term, see :func:`hodograph_cosine` with frequency = :math:`\frac{2\pi(N+1/2)}{T}`
* Power cosine function term, see :func:`hodograph_power_cosine` with  exponent = 3.0, frequency = :math:`\frac{2\pi(N+1/2)}{T}`, scale_factor = :math:`1/T`
* Power sine function term, see :func:`hodograph_power_sine` with  exponent = 3.0, frequency = :math:`\frac{2\pi(N+1/2)}{T}`, scale_factor = :math:`1/T`


Parameters
----------
time_of_flight : float
    Total time of flight (in seconds) of the trajectory that is to be generated.
number_of_revolutions : int
    Number of full revolutions around the central body that are to be used.
Returns
-------
list[BaseFunctionHodographicShaping]
    List of default settings object for axial hodographic shaping






    )";



    } else if(name == "hodograph_constant" && variant==0) {
        return R"(
        
Factory function for creating a constant contribution to hodographic trajectory shaping.

Factory function for creating a constant contribution to hodographic trajectory shaping. This adds a contribution 
:math:`K` to the selected velocity component, with :math:`K` a free parameter.

Returns
-------
BaseFunctionHodographicShaping
    Settings object for a constant contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_sine" && variant==0) {
        return R"(
        
Factory function for creating a sine contribution to hodographic trajectory shaping.

Factory function for creating a sine contribution to hodographic trajectory shaping. For a 
provided frequency :math:`f`, this adds a contribution :math:`K\sin(f\cdot t)` to the selected
velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


Parameters
----------
frequency : float
    Frequency of the sine contribution to the shape function.
Returns
-------
BaseFunctionHodographicShaping
    Settings object for a cosine contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_cosine" && variant==0) {
        return R"(
        
Factory function for creating a cosine contribution to hodographic trajectory shaping.

Factory function for creating a cosine contribution to hodographic trajectory shaping. For a 
provided frequency :math:`f`, this adds a contribution :math:`K\cos(f\cdot T)` to the selected
velocity component, with :math:`T` the time since departure, and :math:`K` a free parameter.


Parameters
----------
frequency : float
    Frequency of the cosine contribution to the shape function.
Returns
-------
BaseFunctionHodographicShaping
    Settings object for a cosine contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_exponential" && variant==0) {
        return R"(
        
Factory function for creating a exponential contribution to hodographic trajectory shaping.

Factory function for creating a exponential contribution to hodographic trajectory shaping. For a 
provided exponent :math:`r` and (optional) scale factor :math:`c`, this adds a contribution :math:`K\exp(cr\cdot t)` to the selected
velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


Parameters
----------
exponent : float
    Exponent of the exponential contribution to the shape function.
scale_factor : float, default = 1.0
    Optional scale factor, which can be used to scale the physical meaning of the exponent :math:`r`.
Returns
-------
BaseFunctionHodographicShaping
    Settings object for a exponential contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_exponential_sine" && variant==0) {
        return R"(
        
Factory function for creating a exponential sine contribution to hodographic trajectory shaping.

Factory function for creating a exponential sine contribution to hodographic trajectory shaping. For a 
provided exponent :math:`r`, scale factor :math:`c` and frequency :math:`f`, this adds a contribution :math:`K\sin(f\cdot t)\exp(cr\cdot t)` to the selected
velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


Parameters
----------
frequency : float
    Frequency of the sine contribution to the shape function.
exponent : float
    Exponent of the exponential contribution to the shape function.
scale_factor : float, default = 1.0
    Optional scale factor, which can be used to scale the physical meaning of the exponent :math:`r`.
Returns
-------
BaseFunctionHodographicShaping
    Settings object for a exponential sine contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_exponential_cosine" && variant==0) {
        return R"(
        
Factory function for creating a exponential cosine contribution to hodographic trajectory shaping.

Factory function for creating a exponential cosine contribution to hodographic trajectory shaping. For a 
provided exponent :math:`r`, scale factor :math:`c`  and frequency :math:`f`, this adds a contribution :math:`K\cos(f\cdot t)\exp(cr\cdot t)` to the selected
velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


Parameters
----------
frequency : float
    Frequency of the cosine contribution to the shape function.
exponent : float
    Exponent of the exponential contribution to the shape function.
scale_factor : float, default = 1.0
    Optional scale factor, which can be used to scale the physical meaning of the exponent :math:`r`.
Returns
-------
BaseFunctionHodographicShaping
    Settings object for a exponential cosine contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_power" && variant==0) {
        return R"(
        
Factory function for creating a power function contribution to hodographic trajectory shaping.

Factory function for creating a power function contribution to hodographic trajectory shaping. For a 
provided exponent :math:`r` and (optional) scale factor :math:`c` this adds a contribution :math:`K\cdot c\cdot t^{r}` to the selected
velocity component, with :math:`t` the time since departure, and :math:`K` a free parameter.


Parameters
----------
exponent : float
    Exponent of the power function contribution to the shape function.
scale_factor : float, default = 1.0
    Optional scale factor, which can be used to scale the physical meaning of the free parameter :math:`K`.
Returns
-------
BaseFunctionHodographicShaping
    Settings object for a power function contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_power_sine" && variant==0) {
        return R"(
        
Factory function for creating a power sine function contribution to hodographic trajectory shaping.

Factory function for creating a power sine function contribution to hodographic trajectory shaping. For a 
provided exponent :math:`r`, (optional) scale factor :math:`c` and frequency :math:`f`, this adds a contribution :math:`K\cdot c\sin(f\cdot t)\cdot t^{r}` to the selected
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
    Settings object for a power sine function contribution to hodographic shaping.






    )";



    } else if(name == "hodograph_power_cosine" && variant==0) {
        return R"(
        
Factory function for creating a power cosine function contribution to hodographic trajectory shaping.

Factory function for creating a power cosine function contribution to hodographic trajectory shaping. For a 
provided exponent :math:`r`, (optional) scale factor :math:`c` and frequency :math:`f`, this adds a contribution :math:`K\cdot c\cos(f\cdot t)\cdot t^{r}` to the selected
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
    Settings object for a power cosine function contribution to hodographic shaping.






    )";



    } else {
        return "No documentation found.";
    }

}


}




}




    
namespace plotting {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";





    } else if(name == "plot_blue_marble_ground_track" && variant==0) {
        return R"(
        





    )";



    } else if(name == "plot_miller_ground_track" && variant==0) {
        return R"(
        





    )";



    } else if(name == "dual_y_axis" && variant==0) {
        return R"(
        





    )";



    } else if(name == "trajectory_3d" && variant==0) {
        return R"(
        





    )";



    } else if(name == "pareto_front" && variant==0) {
        return R"(
        





    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace util {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";




    } else if(name == "redirect_std") {
         return R"(





     )";






    } else if(name == "result2array" && variant==0) {
        return R"(
        





    )";



    } else if(name == "compare_results" && variant==0) {
        return R"(
        





    )";



    } else if(name == "pareto_optimums" && variant==0) {
        return R"(
        





    )";



    } else if(name == "split_history" && variant==0) {
        return R"(
        





    )";



    } else if(name == "vector2matrix" && variant==0) {
        return R"(
        





    )";



    } else {
        return "No documentation found.";
    }

}


}




    
namespace io {

static inline std::string get_docstring(std::string name, int variant=0) {

    if (name == "test") {
        return "test";



    } else if(name == "StaticCoefficientNames") {
         return R"(

        Enumeration of Missile DATCOM static aerodynamic coefficient types.





     )";


    } else if(name == "StaticCoefficientNames.cn") {
         return R"(
Normal force coefficient.
     )";


    } else if(name == "StaticCoefficientNames.cm") {
         return R"(
Pitching moment coefficient.
     )";


    } else if(name == "StaticCoefficientNames.ca") {
         return R"(
Axial force coefficient.
     )";


    } else if(name == "StaticCoefficientNames.cy") {
         return R"(
Side force coefficient.
     )";


    } else if(name == "StaticCoefficientNames.cln") {
         return R"(
Yawing moment coefficient.
     )";


    } else if(name == "StaticCoefficientNames.cll") {
         return R"(
Rolling moment coefficient.
     )";


    } else if(name == "StaticCoefficientNames.cna") {
         return R"(
Normal force coefficient change w.r.t. angle of attack.
     )";


    } else if(name == "StaticCoefficientNames.cma") {
         return R"(
Pitching moment coefficient change w.r.t. angle of attack.
     )";


    } else if(name == "StaticCoefficientNames.cyb") {
         return R"(
Side force coefficient change w.r.t. sideslip angle.
     )";


    } else if(name == "StaticCoefficientNames.cnb") {
         return R"(
Yawing moment coefficient change w.r.t. sideslip angle.
     )";


    } else if(name == "StaticCoefficientNames.clb") {
         return R"(
Rolling moment coefficient change w.r.t. sideslip angle.
     )";



    } else if(name == "DynamicCoefficientNames") {
         return R"(

        Enumeration of Missile DATCOM dynamic aerodynamic coefficient types.





     )";


    } else if(name == "DynamicCoefficientNames.cnq") {
         return R"(
Normal force coefficient change w.r.t. pitch rate.
     )";


    } else if(name == "DynamicCoefficientNames.cmq") {
         return R"(
Pitching moment coefficient change w.r.t. pitch rate.
     )";


    } else if(name == "DynamicCoefficientNames.caq") {
         return R"(
Pitching moment coefficient change w.r.t. pitch rate.
     )";


    } else if(name == "DynamicCoefficientNames.cyq") {
         return R"(
Side force coefficient change w.r.t. pitch rate.
     )";


    } else if(name == "DynamicCoefficientNames.clnq") {
         return R"(
Yawing moment coefficient change w.r.t. pitch rate.
     )";


    } else if(name == "DynamicCoefficientNames.cllq") {
         return R"(
Rolling moment coefficient change w.r.t. pitch rate.
     )";


    } else if(name == "DynamicCoefficientNames.cnr") {
         return R"(
Normal force coefficient change w.r.t. yaw rate.
     )";


    } else if(name == "DynamicCoefficientNames.cmr") {
         return R"(
Pitching moment coefficient change w.r.t. yaw rate.
     )";


    } else if(name == "DynamicCoefficientNames.car") {
         return R"(
Pitching moment coefficient change w.r.t. yaw rate.
     )";


    } else if(name == "DynamicCoefficientNames.cyr") {
         return R"(
Side force coefficient change w.r.t. yaw rate.
     )";


    } else if(name == "DynamicCoefficientNames.clnr") {
         return R"(
Yawing moment coefficient change w.r.t. yaw rate.
     )";


    } else if(name == "DynamicCoefficientNames.cllr") {
         return R"(
Rolling moment coefficient change w.r.t. yaw rate.
     )";


    } else if(name == "DynamicCoefficientNames.cnp") {
         return R"(
Normal force coefficient change w.r.t. roll rate.
     )";


    } else if(name == "DynamicCoefficientNames.cmp") {
         return R"(
Pitching moment coefficient change w.r.t. roll rate.
     )";


    } else if(name == "DynamicCoefficientNames.cap") {
         return R"(
Pitching moment coefficient change w.r.t. roll rate.
     )";


    } else if(name == "DynamicCoefficientNames.cyp") {
         return R"(
Side force coefficient change w.r.t. roll rate.
     )";


    } else if(name == "DynamicCoefficientNames.clnp") {
         return R"(
Yawing moment coefficient change w.r.t. roll rate.
     )";


    } else if(name == "DynamicCoefficientNames.cllp") {
         return R"(
Rolling moment coefficient change w.r.t. roll rate.
     )";


    } else if(name == "DynamicCoefficientNames.cnad") {
         return R"(
Normal force coefficient change w.r.t. angle of attack change rate.
     )";


    } else if(name == "DynamicCoefficientNames.cmad") {
         return R"(
Pitching moment coefficient change w.r.t. angle of attack change rate.
     )";




    } else if(name == "missile_DATCOM_data") {
         return R"(

        Class containing data and methods interfacing the Missile DATCOM software.

        This class is the main method that can be used to interface tudat with the Missile DATCOM software.
        It can be initialised with the output file from Missile DATCOM, and provides methods to convert these results
        into tudat-compatible data.

        .. note:: The Missile DATCOM software from which outputs can be interfaced to TUDAT is an entirely separate software from Tudat(Py).
                  Please refer to Missile DATCOM user manuals for information on how to use it. These can be accessed on the US Defence Technical
                  Information Center at accession numbers `ADA267447 <https://apps.dtic.mil/sti/citations/ADA267447>`_ and
                  `ADA503576 <https://apps.dtic.mil/sti/citations/ADA503576>`_.

        .. note:: The interfacing of Missile DATCOM to tudat assumes that aerodynamic coefficients are computed as a function of both
                  Mach number and angle of attack.





     )";




    } else if(name == "missile_DATCOM_data.ctor" && variant==0) {
            return R"(

        Class constructor.

        Function used to construct and initialise the class. In essence, it can be used to read and extract the aerodynamic coefficients
        computed by Missile DATCOM, and save them in different formats.


        Parameters
        ----------
        file_name_and_path : str
            Full path and file name of the `for004.dat` Missile DATCOM results output file.




    )";



    } else if(name == "missile_DATCOM_data.get_static_coefficient" && variant==0) {
            return R"(

        Get a specific static coefficient from the result database.


        Parameters
        ----------
        mach_index : int
            Index of the Mach number for which to get the static coefficient.
        angle_of_attack_index : int
            Index of the angle of attack for which to get the static coefficient.
        coefficient_index : tudatpy.io.StaticCoefficientNames
            Type of the static aerodynamic coefficient.
        Returns
        -------
        float
            Static aerodynamic coefficient.





    )";



    } else if(name == "missile_DATCOM_data.get_dynamic_coefficient" && variant==0) {
            return R"(

        Get a specific dynamic coefficient from the result database.


        Parameters
        ----------
        mach_index : int
            Index of the Mach number for which to get the static coefficient.
        angle_of_attack_index : int
            Index of the angle of attack for which to get the static coefficient.
        coefficient_index : tudatpy.io.DynamicCoefficientNames
            Type of the dynamic aerodynamic coefficient.
        Returns
        -------
        float
            Dynamic aerodynamic coefficient.





    )";



    } else if(name == "missile_DATCOM_data.get_angle_of_attacks" && variant==0) {
            return R"(

        Get the list of angle of attacks at which Missile DATCOM has been run.

        Returns
        -------
        numpy.ndarray
            List of angle of attacks.





    )";



    } else if(name == "missile_DATCOM_data.get_mach_numbers" && variant==0) {
            return R"(

        Get the list of Mach numbers at which Missile DATCOM has been run.

        Returns
        -------
        numpy.ndarray
            List of Mach numbers.





    )";



    } else if(name == "missile_DATCOM_data.get_Reynolds_numbers" && variant==0) {
            return R"(

        Get the list of Reynolds numbers at which Missile DATCOM has been run.

        Returns
        -------
        numpy.ndarray
            List of Reynolds numbers.





    )";



    } else if(name == "missile_DATCOM_data.write_all_coefficients_to_files" && variant==0) {
            return R"(

        Write all the aerodynamic coefficients to CSV files.


        Parameters
        ----------
        file_name_base : str
            Full base path and name of the file that will be saved. The name of each aerodynamic coefficient will be included at the end of the file name.
        base_precision : int, optional, default=15
            Number of digits to represent the base of the floating-point number.
        exponent_width : int, optional, default=2
            Number of digits to represent the exponent of the floating-point number.




    )";



    } else if(name == "missile_DATCOM_data.write_force_and_moment_coefficients_to_files" && variant==0) {
            return R"(

        Write the force and moment coefficients to a file in the format taken by the :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_from_files` function.


        Parameters
        ----------
        file_name_base : str
            Full base path and name of the file that will be saved. The name of each aerodynamic coefficient will be included at the end of the file name.
        base_precision : int, optional, default=15
            Number of digits to represent the base of the floating-point number.
        exponent_width : int, optional, default=2
            Number of digits to represent the exponent of the floating-point number.




    )";





    } else if(name == "save2txt" && variant==0) {
        return R"(
        





    )";



    } else if(name == "save_time_history_to_file" && variant==0) {
        return R"(
        





    )";



    } else if(name == "get_resource_path" && variant==0) {
        return R"(
        
Get the path at which tudat resources are located.

Returns
-------
str
    Local path at which tudat resources are located.






    )";



    } else if(name == "get_ephemeris_path" && variant==0) {
        return R"(
        
Get the path at which the ephemeris used by tudat are located.

Returns
-------
str
    Local path at which the tudat ephemeris resources are located.






    )";



    } else if(name == "get_earth_orientation_path" && variant==0) {
        return R"(
        
Get the path at which the Earth orientation resources used by tudat are located.

Returns
-------
str
    Local path at which tudat Earth orientation resources are located.






    )";



    } else if(name == "get_quadrature_path" && variant==0) {
        return R"(
        
Get the path at which the Gaussian quadrature resources are located.

Returns
-------
str
    Local path at which tudat Gaussian quadrature resources are located.






    )";



    } else if(name == "get_spice_kernel_path" && variant==0) {
        return R"(
        
Get the path at which the SPICE kernel used by tudat is located.

Returns
-------
str
    Local path at which the SPICE kernel is located.






    )";



    } else if(name == "get_atmosphere_tables_path" && variant==0) {
        return R"(
        
Get the path at which tudat atmosphere tables are located.

Returns
-------
str
    Local path at which tudat atmosphere tables are located.






    )";



    } else if(name == "get_gravity_models_path" && variant==0) {
        return R"(
        
Get the path at which tudat gravity models are located.

Returns
-------
str
    Local path at which tudat gravity models are located.






    )";



    } else if(name == "get_space_weather_path" && variant==0) {
        return R"(
        
Get the path at which tudat space weather is located.

Returns
-------
str
    Local path at which tudat space weather is located.






    )";



    } else if(name == "read_vector_history_from_file" && variant==0) {
        return R"(
        
Read a vector history from a file.


Parameters
----------
vector_size : int
    Size of the vector at each epoch.
file_name : str
    Name of the file containing the vector history.
Returns
-------
Dict[float, numpy.ndarray]
    Dictionary mapping epochs to the vector at the given epoch.






    )";



    } else if(name == "read_matrix_history_from_file" && variant==0) {
        return R"(
        
Read a matrix history from a file.


Parameters
----------
matrix_rows : int
    Number of rows in the matrix at each epoch.
matrix_columns : int
    Number of columns in the matrix at each epoch.
file_name : str
    Name of the file containing the matrix history.
Returns
-------
Dict[float, numpy.ndarray]
    Dictionary mapping epochs to the matrix at the given epoch.






    )";



    } else {
        return "No documentation found.";
    }

}


}




}

