/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_rotation_model_setup.h"

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup.h>

// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{
namespace rotation_model
{

void expose_rotation_model_setup( py::module &m )
{
    /////////////////////////////////////////////////////////////////////////////
    // createRotationalModel.h
    /////////////////////////////////////////////////////////////////////////////
    py::enum_<tss::RotationModelType>( m,
                                       "RotationModelType",
                                       R"doc(

        Enumeration of rotation model types.

        Enumeration of rotation model types supported by tudat.





     )doc" )
        .value( "simple_rotational_model", tss::RotationModelType::simple_rotation_model,
                R"doc(No documentation found.)doc" )
        .value( "spice_rotation_model",
                tss::RotationModelType::spice_rotation_model,
                R"doc(
     )doc" )
        .value( "gcrs_to_itrs_rotation_model",
                tss::RotationModelType::gcrs_to_itrs_rotation_model,
                R"doc(
     )doc" )
        .value( "synchronous_rotation_model",
                tss::RotationModelType::synchronous_rotation_model,
                R"doc(
     )doc" )
        .value( "planetary_rotation_model",
                tss::RotationModelType::planetary_rotation_model,
                R"doc(
     )doc" )
        .export_values( );

    py::enum_<tba::IAUConventions>( m,
                                    "IAUConventions",
                                    R"doc(

        Enumeration of IAU conventions for Earth rotation.

        Enumeration of IAU conventions for Earth rotation supported by tudat.





     )doc" )
        .value( "iau_2000_a",
                tba::IAUConventions::iau_2000_a,
                R"doc(
     )doc" )
        .value( "iau_2000_b",
                tba::IAUConventions::iau_2000_b,
                R"doc(
     )doc" )
        .value( "iau_2006",
                tba::IAUConventions::iau_2006,
                R"doc(
     )doc" )
        .export_values( );

    py::class_<tss::RotationModelSettings, std::shared_ptr<tss::RotationModelSettings> >( m, "RotationModelSettings", R"doc(

        Base class for providing settings for automatic rotation model creation.

        This class is a functional base class for settings of rotation models that require no information in addition to their type.
        Basic rotation model has constant orientation of the rotation axis (body-fixed z-axis) and constant rotation rate about this axis.
        Rotation models requiring additional information must be created using the functions which create the specific object derived from this base class.





     )doc" )
        //            .def(py::init<const
        //            tss::RotationModelType, const std::string
        //            &,
        //                 const std::string &>(),
        //                 py::arg("rotation_type"),
        //                 py::arg("base_frame"),
        //                 py::arg("target_frame"))
        .def_property_readonly( "rotation_type",
                                &tss::RotationModelSettings::getRotationType,
                                R"doc(

        **read-only**

        Type of rotation model that is to be created.

        :type: RotationModelType
     )doc" )
        .def_property( "base_frame",
                       &tss::RotationModelSettings::getOriginalFrame,
                       &tss::RotationModelSettings::resetOriginalFrame,
                       R"doc(

        Name of the base frame of rotation model.

        :type: str
     )doc" )
        .def_property_readonly( "target_frame",
                                &tss::RotationModelSettings::getTargetFrame,
                                R"doc(

        **read-only**

        Name of the target frame of rotation model.

        :type: str
     )doc" );

    py::class_<tss::SimpleRotationModelSettings, std::shared_ptr<tss::SimpleRotationModelSettings>, tss::RotationModelSettings>(
        m, "SimpleRotationModelSettings", R"doc(No documentation found.)doc" );

    py::class_<tss::PlanetaryRotationModelSettings, std::shared_ptr<tss::PlanetaryRotationModelSettings>, tss::RotationModelSettings>(
        m, "PlanetaryRotationModelSettings", R"doc(No documentation found.)doc" );

    py::class_<tss::IauRotationModelSettings, std::shared_ptr<tss::IauRotationModelSettings>, tss::RotationModelSettings>(
        m, "IAURotationModelSettings", R"doc(No documentation found.)doc" );


    m.def( "simple",
           py::overload_cast<const std::string &, const std::string &, const Eigen::Matrix3d &, const double, const double>(
               &tss::simpleRotationModelSettings ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "initial_orientation" ),
           py::arg( "initial_time" ),
           py::arg( "rotation_rate" ),
           R"doc(

Function for creating simple rotation model settings.

Function for settings object, defining a basic rotation model with constant orientation of the rotation axis and constant rotation rate about this axis.
Rotation from original (inertial) to target (body-fixed) frame at some reference time ``initial_time`` (:math:`t_{0}`) is defined by the ``initial_orientation`` (:math:`\mathbf{R}^{(B/I)}(t_{0})`) rotation matrix.
Rotation about the body-fixed z-axis is defined by the ``rotation_rate`` (:math:`\omega`) float variable (in rad/s). The rotation matrix is computed from:

.. math::
   \mathbf{R}^{(B/I)}(t)=\mathbf{R}_{z}(\omega(t-t_{0}))(t_{0})\mathbf{R}^{(B/I)}(t_{0})

where :math:`\mathbf{R}^{(B/I)}` denotes the rotation matrix from inertial to body-fixed frame, and :math:`\mathbf{R}_{z}` denotes a rotation matrix about the z-axis.

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


    )doc" );

    m.def( "simple_from_spice",
           &tss::simpleRotationModelFromSpiceSettings,
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "target_frame_spice" ),
           py::arg( "initial_time" ),
           R"doc(

Function for creating simple rotation model settings using initial orientation and rotation rates from Spice.

Function for settings object, defining a :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` rotation model with the added functionality that the initial orientation and rotation rate are extracted from Spice, as opposed to provided manually.
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


    )doc" );

    m.def( "synchronous",
           &tss::synchronousRotationModelSettings,
           py::arg( "central_body_name" ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           R"doc(

Function for creating synchronous rotational ephemeris settings.

Function for settings object, defining a synchronous rotation model where rotation of a body is defined from its relative orbit w.r.t. some central body. Specifically
- the body-fixed x-axis is *always* pointing towards the central body
- the body-fixed z-axis is *always* perpendicular to the orbital plane (along the direction of :math:`\mathbf{x}\times\mathbf{v}` )
- the body-fixed y-axis completes the right-handed reference frame

Such a model can be useful for, for instance, approximate rotation of tidally locked natural satellites or nadir-pointing spacecraft.


Parameters
----------
central_body_name : str
    Name of the central body of synchronous rotation.
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


    )doc" );

    m.def( "spice",
           &tss::spiceRotationModelSettings,
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "spice_frame_name" ) = "",
           R"doc(

Function for creating rotation model settings from the Spice interface.

Function for settings object, defining a rotation model directly (and entirely) from Spice interface.


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


    )doc" );

    m.def( "gcrs_to_itrs",
           &tss::gcrsToItrsRotationModelSettings,
           py::arg( "precession_nutation_theory" ) = tba::iau_2006,
           py::arg( "base_frame" ) = "GCRS",
           py::arg( "cio_interpolation_settings" ) = nullptr,
           py::arg( "tdb_to_tt_interpolation_settings" ) = nullptr,
           py::arg( "short_term_eop_interpolation_settings" ) = nullptr,
           R"doc(

Function for creating high-accuracy Earth rotation model settings.

Function for settings object, defining high-accuracy Earth rotation model according to the IERS Conventions 2010.  The model
computes the rotation from ITRS to GCRS (with rotation matrix :math:`\mathbf{R}^{(\text{GCRS}/\text{ITRS})}`) and its inverse from:

.. math::
   \mathbf{R}^{(\text{GCRS}/\text{ITRS})} = \mathbf{R}^{(\text{GCRS}/\text{CIRS})}(X,Y,s)\mathbf{R}^{(\text{CIRS}/\text{TIRS})}(\theta_{E})\mathbf{R}^{(\text{TIRS}/\text{ITRS})}(x_{p}, y_{p}, s')

using the intermediate frames TIRS (Terrestial Intermediate Reference System) and CIRS (Celestial Intermediate Reference System) where (with equations referring to the IERS 2010 Conventions) :math:`\mathbf{R}^{(\text{GCRS}/\text{CIRS})}` implements Eq. (5.10), :math:`\mathbf{R}^{(\text{CIRS}/\text{TIRS})}` implements Eq. (5.5), and
:math:`\mathbf{R}^{(\text{TIRS}/\text{ITRS})}` implements Eq. (5.3). The inputs to these rotation matrices are :

* :math:`X`, :math:`Y`: Celestial pole position elements
* :math:`s`: The CIO (Celestial Intermediate Origin)
* :math:`\theta_{E}`: Earth rotation angle (denoted as :math:`ERA` in IERS Conventions)
* :math:`x_{p}`, :math:`y_{p}`: Polar motion components
* :math:`s'`: The TIO (Terrestial Intermediate Origin)

Depending on the selected ``precession_nutation_theory`` input, the SOFA function ``iauXys00a``, ``iauXys00b`` or ``iauXys06a`` is used to compute :math:`X,Y,s`, when selecting
:class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.IAUConventions` ``iau_2000a``, ``iau_2000b`` or ``iau_2006``, respectively. Corrections to the nominal values of :math:`X,Y`
are applied using linear interpolation of daily corrections for :math:`X,Y` from the eopc04_14_IAU2000.62-now.txt file. The quantity :math:`s'` is computed from Eq. (5.13) (implemented in SOFA's ``iauSp00`` function).

The value of :math:`\theta_{E}` is computed directly from UTC-UT1, which is computed using settings given in :func:`~tudatpy.astro.time_conversions.default_time_scale_converter`, the computation of
:math:`theta_{E}` from this quantity follows from Eq. (5.15), implemented by SOFA's ``iauEra00`` function.

The polar motion components :math:`x_{p}`, :math:`y_{p}` are computed from:

* Corrections for semi-diurnal variations due to libration for a non-rigid Earth as per Table 5.1a (with :math:`n=2`) of IERS Conventions 2010
* Corrections diurnal and semidiurnal variations due to ocean tides as per Tables 8.1a and 8.1b of the IERS Conventions 2010
* Linear interpolation (correcting for discontinuities during days with leap seconds) of daily corrections for :math:`x_{p}, y_{p}`: from the eopc04_14_IAU2000.62-now.txt file in the tudat-resources directory

Note that for this model the original frame must be J2000 or GCRS (in the case of the former, the frame bias between GCRS and J2000 is automatically corrected for). The target frame (e.g. body-fixed frame) name is ITRS.
The target frame (e.g. body-fixed frame) name is ITRS.

Alternative options to modify the input (not exposed here) include the EOP correction file, input time scale, short period UT1 and polar motion variations.

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
using a high-accuracy Earth rotation model as defined by IERS Conventions 2010:


.. code-block:: python

   # define parameters describing the rotation between frames
   precession_nutation_theory = environment_setup.rotation_model.IAUConventions.iau_2006
   original_frame = "J2000"
   # create rotation model settings and assign to body settings of "Earth"
   body_settings.get( "Earth" ).rotation_model_settings = environment_setup.rotation_model.gcrs_to_itrs(
   precession_nutation_theory, original_frame)


    )doc" );

    m.def( "aerodynamic_angle_based",
           &tss::aerodynamicAngleRotationSettings,
           py::arg( "central_body" ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "angle_funcion" ) = nullptr,
           R"doc(

Function for creating rotation model settings based on custom aerodynamic angles (attack, sideslip, bank).

Function for creating rotation model settings based on custom aerodynamic angles:
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






    )doc" );

    m.def( "zero_pitch_moment_aerodynamic_angle_based",
           &tss::pitchTrimRotationSettings,
           py::arg( "central_body" ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "angle_funcion" ) = nullptr,
           R"doc(

Function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank.

Function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank. This function is
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






    )doc" );

    m.def( "custom_inertial_direction_based",
           &tss::bodyFixedDirectionBasedRotationSettings,
           py::arg( "inertial_body_axis_direction" ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "free_rotation_angle_function" ) = nullptr,
           R"doc(

Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction

Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction :math:`\hat{\mathbf{T}}_{I}`. Specifically, it ensures
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






    )doc" );

    m.def( "orbital_state_direction_based",
           &tss::orbitalStateBasedRotationSettings,
           py::arg( "central_body" ),
           py::arg( "is_colinear_with_velocity" ),
           py::arg( "direction_is_opposite_to_vector" ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ) = "",
           py::arg( "free_rotation_angle_function" ) = nullptr,
           R"doc(

Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector.

Function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector. This function is
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






    )doc" );

    m.def( "constant_rotation_model",
           py::overload_cast<const std::string &, const std::string &, const Eigen::Matrix3d &>(
               &tss::constantRotationModelSettings ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "initial_orientation" ),
           R"doc(

Function for creating simple rotation model settings for target-frames with constant orientation.

Function for settings object, defining simple rotation model setting objects with constant rotation matrix.
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


    )doc" );

    m.def( "custom_rotation_model",
           &tss::customRotationModelSettings,
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "custom_rotation_matrix_function" ),
           py::arg( "finite_difference_time_step" ),
           R"doc(

Function for creating rotation model settings based on custom definition of rotation matrix.

Function for creating rotation model settings based on custom definition of rotation matrix. The user provides a custom function that computes the rotation matrix
from body-fixed to inertial frame as a function of time. This function can
depend on any quantities of the user's choosing, for details on how to link the properties of the environment to this function, see `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/custom_models.html>`_.
Since this function only computes the rotation matrix directly, the rotation matrix time derivative (and consequently, the angular velocity) are computed numerically, using a second
order finite-difference method. Note that this computation of time derivative will only take into account the explicit time-dependence of thh custom rotation matrix.

Parameters
----------
base_frame : str
    Name of the base frame of rotation model.
target_frame : str
    Name of the target frame of rotation model.
custom_rotation_matrix_function: Callable[[float], numpy.ndarray[numpy.float64[3, 3]]]
    Function computing the body-fixed to inertial rotation matrix as a function of time
finite_difference_time_step: float
    Step size to use when computing the rotation matrix derivative numerically
Returns
-------
:class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings`
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model.



    )doc" );

    m.def( "mars_high_accuracy",
           &tss::getHighAccuracyMarsRotationModel,
           py::arg( "base_frame" ) = "ECLIPJ2000",
           py::arg( "target_frame" ) = "Mars_Fixed",
           R"doc(

Function for creating a high-accuracy Mars rotation model.

Function for creating a high-accuracy Mars rotation model, using the default parameters of `Konopliv et al. (2016) <https://www.sciencedirect.com/science/article/abs/pii/S0019103516001305>`_
and the mathematical model of `Konopliv et al. (2006) <https://www.sciencedirect.com/science/article/pii/S0019103506000297>`_ . The rotation matrix formulation is given in Eq. (13)-(19) of that paper.
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






    )doc" );

    m.def( "iau_rotation_model",
           &tss::iauRotationModelSettings,
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "nominal_meridian" ),
           py::arg( "nominal_pole" ),
           py::arg( "rotation_rate" ),
           py::arg( "pole_precession" ),
           py::arg( "merdian_periodic_terms" ),
           py::arg( "pole_periodic_terms" )
    );
}
}  // namespace rotation_model
}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
