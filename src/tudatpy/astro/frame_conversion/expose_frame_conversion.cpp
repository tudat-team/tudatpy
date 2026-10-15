/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_frame_conversion.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/astro/reference_frames/terrestrialFrameTransformation.h>

namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;

namespace py = pybind11;

namespace tudatpy
{

namespace astro
{
namespace frame_conversion
{

// Register Cartesian frame conversions and terrestrial frame transformations.
void expose_frame_conversion( py::module& m )
{
    m.def( "inertial_to_rsw_rotation_matrix",
           &trf::getInertialToRswSatelliteCenteredFrameRotationMatrix,
           py::arg( "inertial_cartesian_state" ),
           R"doc(

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






     )doc" );

    m.def( "rsw_to_inertial_rotation_matrix",
           &trf::getRswSatelliteCenteredToInertialFrameRotationMatrix,
           py::arg( "inertial_cartesian_state" ),
           R"doc(

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






     )doc" );

    m.def( "tnw_to_inertial_rotation_matrix",
           py::overload_cast< const Eigen::Vector6d&, const bool >( &trf::getTnwToInertialRotation ),
           py::arg( "inertial_cartesian_state" ),
           py::arg( "n_axis_points_away_from_central_body" ) = true,
           R"doc(

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






     )doc" );

    m.def( "inertial_to_tnw_rotation_matrix",
           py::overload_cast< const Eigen::Vector6d&, const bool >( &trf::getInertialToTnwRotation ),
           py::arg( "inertial_cartesian_state" ),
           py::arg( "n_axis_points_away_from_central_body" ) = true,
           R"doc(

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






     )doc" );

    m.def( "inertial_to_body_fixed_rotation_matrix",
           py::overload_cast< const double, const double, const double >( &trf::getInertialToPlanetocentricFrameTransformationMatrix ),
           py::arg( "pole_declination" ),
           py::arg( "pole_right_ascension" ),
           py::arg( "prime_meridian_longitude" ),
           R"doc(

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
 Euler angle rotation (see :cite:t:`Archinal2018`).


     )doc" );

    m.def( "body_fixed_to_inertial_rotation_matrix",
           py::overload_cast< const double, const double, const double >(
                   &trf::getRotatingPlanetocentricToInertialFrameTransformationMatrix ),
           py::arg( "pole_declination" ),
           py::arg( "pole_right_ascension" ),
           py::arg( "pole_meridian" ),
           R"doc(

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
 Euler angle rotation (see :cite:t:`Archinal2018`).




     )doc" );

    m.def( "rotate_state_to_frame",
           py::overload_cast< const Eigen::Vector6d&, const Eigen::Matrix3d&, const Eigen::Matrix3d& >(
                   &te::transformStateToFrameFromRotations< double > ),
           py::arg( "original_state" ),
           py::arg( "rotation_matrix" ),
           py::arg( "rotation_matrix_time_derivative" ) = Eigen::Matrix3d::Zero( ),
           R"doc(

 Rotates a Cartesian state (position and velocity) from one frame :math:`B` to another frame :math:`A`, using the rotation matrix :math:`\mathbf{R}^{(A/B)}` from frame :math:`B` to :math:`A`, and its time derivative
 :math:`\dot{\mathbf{R}}^{(A/B)}`.
 This function computes:

 .. math::
    \mathbf{r}^{(A)}=\mathbf{R}^{(A/B)}\mathbf{r}^{(B)}+\dot{\mathbf{R}}^{(A/B)}\mathbf{v}^{(B)}\\
    \mathbf{v}^{(A)}=\mathbf{R}^{(A/B)}\mathbf{v}^{(B)}

 Parameters
 ----------
 original_state : ndarray[numpy.float64[6, 1]]
     Cartesian state vector :math:`\mathbf{x}^{(B)}=[\mathbf{r}^{(B)};\mathbf{v}^{(B)}]` in frame :math:`B`
 rotation_matrix: ndarray[numpy.float64[3, 3]]
     Rotation matrix :math:`\mathbf{R}^{(A/B)}` from frame :math:`B` to :math:`A`
 rotation_matrix_time_derivative: ndarray[numpy.float64[3, 3]], default = numpy.zeros((3, 3))
     Time derivative of rotation matrix (:math:`\dot{\mathbf{R}}^{(A/B)})` from frame :math:`B` to :math:`A`; default zero indicates that frames :math:`A` and :math:`B` have a constant orientation w.r.t. one another.
 Returns
 -------
 numpy.ndarray
     Input state in frame :math:`B`, rotated to frame :math:`A`


     )doc" );

    // Terrestrial catalogue operations are numerical frame conversions, independent of environment settings.
    py::class_< trf::TerrestrialFrameAffineTransformation >( m, "TerrestrialFrameAffineTransformation", R"doc(

Transformation between terrestrial frame realizations described by a translation, a matrix and their time derivatives.

A position :math:`\mathbf{x}` is transformed as :math:`\mathbf{y}=\mathbf{T}+\mathbf{A}\mathbf{x}`, and a state as
:math:`\dot{\mathbf{y}}=\dot{\mathbf{T}}+\dot{\mathbf{A}}\mathbf{x}+\mathbf{A}\dot{\mathbf{x}}`. For a Helmert
transformation, the matrix :math:`\mathbf{A}` contains the scale and the linearized rotation, and is therefore not
orthogonal. Objects of this class are created by :func:`evaluate_linearized_helmert` and
:meth:`HelmertTransformationModel.get_transformation`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "translation", &trf::TerrestrialFrameAffineTransformation::translation, R"doc(

         Translation :math:`\mathbf{T}` in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "matrix", &trf::TerrestrialFrameAffineTransformation::matrix, R"doc(

         Matrix :math:`\mathbf{A}` (dimensionless).

         :type: numpy.ndarray[numpy.float64[3, 3]]
      )doc" )
            .def_readwrite( "translation_rate", &trf::TerrestrialFrameAffineTransformation::translationRate, R"doc(

         Time derivative of the translation, in metres per second (of the time used for the input velocity).

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "matrix_rate", &trf::TerrestrialFrameAffineTransformation::matrixRate, R"doc(

         Time derivative of the matrix, per second (of the time used for the input velocity).

         :type: numpy.ndarray[numpy.float64[3, 3]]
      )doc" )
            .def( "transform_position",
                  &trf::TerrestrialFrameAffineTransformation::transformPosition,
                  py::arg( "position" ),
                  R"doc(

         Function to transform a Cartesian position.

         Parameters
         ----------
         position : numpy.ndarray[numpy.float64[3, 1]]
             Cartesian position in the source frame, in metres.

         Returns
         -------
         numpy.ndarray[numpy.float64[3, 1]]
             Cartesian position in the target frame, in metres.

      )doc" )
            .def( "transform_state",
                  &trf::TerrestrialFrameAffineTransformation::transformState,
                  py::arg( "state" ),
                  R"doc(

         Function to transform a Cartesian state, including the time derivative of the transformation.

         Parameters
         ----------
         state : numpy.ndarray[numpy.float64[6, 1]]
             Cartesian state in the source frame, in metres and metres per second.

         Returns
         -------
         numpy.ndarray[numpy.float64[6, 1]]
             Cartesian state in the target frame, in metres and metres per second.

      )doc" )
            .def( "inverse", &trf::TerrestrialFrameAffineTransformation::inverse, R"doc(

         Function to compute the inverse transformation (using the exact inverse of the matrix, and its derivative).

         Returns
         -------
         TerrestrialFrameAffineTransformation
             Transformation from the target frame to the source frame.

      )doc" )
            .def( "compose_after",
                  &trf::TerrestrialFrameAffineTransformation::composeAfter,
                  py::arg( "before" ),
                  R"doc(

         Function to compute the transformation that first applies ``before``, and then this transformation.

         Parameters
         ----------
         before : TerrestrialFrameAffineTransformation
             Transformation that is applied first.

         Returns
         -------
         TerrestrialFrameAffineTransformation
             Combined transformation, including the derivatives of both transformations.

      )doc" );

    m.def( "evaluate_linearized_helmert",
           &trf::evaluateLinearizedHelmert,
           py::arg( "parameters" ),
           py::arg( "decimal_year" ),
           py::arg( "decimal_year_rate_with_respect_to_input_seconds" ),
           R"doc(

Function to evaluate a linearized Helmert transformation at a given epoch.

Function to evaluate the 14-parameter (linearized) Helmert transformation defined by a set of published parameters at
a given epoch. The parameters are propagated linearly from their reference epoch to the requested decimal year. No
time scale or calendar conversion is performed: the epoch is given directly as a decimal year.


Parameters
----------
parameters : TerrestrialFrameTransformationParameters
    Published transformation parameters (see :ref:`terrestrial_frame_data`).
decimal_year : float
    Epoch at which the transformation is evaluated, as a decimal year.
decimal_year_rate_with_respect_to_input_seconds : float
    Derivative of the decimal year with respect to the time (in seconds) used for the velocities that are to be
    transformed (for instance :math:`1/31557600` for a Julian year). Used only for the derivative of the transformation.

Returns
-------
TerrestrialFrameAffineTransformation
    Transformation from the source to the target frame of the parameters.

      )doc" );

    py::class_< trf::HelmertTransformationModel >( m, "HelmertTransformationModel", R"doc(

Model for linearized Helmert transformations with parameter rates between terrestrial frame realizations.

The model combines the transformations of the table through its reference realization, so that any two realizations
in the table can be transformed into each other. A transformation between two identical realization names is the
identity. No time scale or calendar conversion is performed: epochs are given directly as decimal years.

      )doc" )
            .def( py::init< const tudat::input_output::TerrestrialFrameTransformationData& >( ),
                  py::arg( "data" ),
                  R"doc(

         Constructor.

         Parameters
         ----------
         data : TerrestrialFrameTransformationData
             Table of transformation parameters (see :ref:`terrestrial_frame_data`).

      )doc" )
            .def( "get_transformation",
                  &trf::HelmertTransformationModel::getTransformation,
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "decimal_year" ),
                  py::arg( "decimal_year_rate_with_respect_to_input_seconds" ),
                  R"doc(

         Function to compute the transformation between two realizations at a given epoch.

         Parameters
         ----------
         source_realization : str
             Name of the realization from which to transform.
         target_realization : str
             Name of the realization to which to transform.
         decimal_year : float
             Epoch of the transformation, as a decimal year.
         decimal_year_rate_with_respect_to_input_seconds : float
             Derivative of the decimal year with respect to the time (in seconds) used for the input velocities.

         Returns
         -------
         TerrestrialFrameAffineTransformation
             Transformation from the source to the target realization.

      )doc" )
            .def( "transform_position",
                  &trf::HelmertTransformationModel::transformPosition,
                  py::arg( "position" ),
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "decimal_year" ),
                  R"doc(

         Function to transform a Cartesian position between two realizations.

         Parameters
         ----------
         position : numpy.ndarray[numpy.float64[3, 1]]
             Cartesian position in the source realization, in metres.
         source_realization : str
             Name of the realization from which to transform.
         target_realization : str
             Name of the realization to which to transform.
         decimal_year : float
             Epoch of the transformation, as a decimal year.

         Returns
         -------
         numpy.ndarray[numpy.float64[3, 1]]
             Cartesian position in the target realization, in metres.

      )doc" )
            .def( "transform_state",
                  &trf::HelmertTransformationModel::transformState,
                  py::arg( "state" ),
                  py::arg( "source_realization" ),
                  py::arg( "target_realization" ),
                  py::arg( "decimal_year" ),
                  py::arg( "decimal_year_rate_with_respect_to_input_seconds" ),
                  R"doc(

         Function to transform a Cartesian state between two realizations.

         Parameters
         ----------
         state : numpy.ndarray[numpy.float64[6, 1]]
             Cartesian state in the source realization, in metres and metres per second.
         source_realization : str
             Name of the realization from which to transform.
         target_realization : str
             Name of the realization to which to transform.
         decimal_year : float
             Epoch of the transformation, as a decimal year.
         decimal_year_rate_with_respect_to_input_seconds : float
             Derivative of the decimal year with respect to the time (in seconds) used for the input velocity.

         Returns
         -------
         numpy.ndarray[numpy.float64[6, 1]]
             Cartesian state in the target realization, in metres and metres per second.


         Examples
         --------
         .. code-block:: python

            from tudatpy.data_input.environment_data import terrestrial_frames
            from tudatpy.astro import frame_conversion

            model = frame_conversion.HelmertTransformationModel(
                terrestrial_frames.official_itrf2020_parameters( ) )
            target_state = model.transform_state(
                source_state, "ITRF2020", "ITRF2014", 2025.0, 1.0 / 31557600.0 )

      )doc" )
            .def_property_readonly(
                    "data",
                    []( const trf::HelmertTransformationModel& model ) { return model.getData( ); },
                    R"doc(

         **read-only**

         Table of transformation parameters used by the model.

         :type: TerrestrialFrameTransformationData
      )doc" );

    m.attr( "transform_cartesian_state_to_frame" ) = m.attr( "rotate_state_to_frame" );
}

}  // namespace frame_conversion
}  // namespace astro
}  // namespace tudatpy
