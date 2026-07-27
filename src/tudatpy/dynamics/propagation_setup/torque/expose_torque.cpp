/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "expose_torque.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/deprecationWarnings.h>

#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/propagation_setup/createMassRateModels.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/environmentUpdater.h"
#include "tudat/simulation/propagation_setup/propagationOutput.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
#include "tudat/simulation/propagation_setup/torqueSettings.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudat
{
namespace simulation_setup
{
inline std::shared_ptr< TorqueSettings > customTorqueSettingsDeprecated(
        const std::function< Eigen::Vector3d( const double ) > torqueFunction,
        const std::function< double( const double ) > scalingFunction = nullptr )
{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning(
                "tudatpy.dynamics.propagation_setup."
                "acceleration.custom",
                "tudatpy.dynamics.propagation_setup."
                "acceleration.custom_torque" );
        isWarningPrinted = true;
    }

    return customTorqueSettings( torqueFunction, scalingFunction );
}
}  // namespace simulation_setup
}  // namespace tudat
namespace tudatpy
{
namespace dynamics
{
namespace propagation_setup
{
namespace torque
{

void expose_torque_setup( py::module& m )
{
    py::enum_< tba::AvailableTorque >( m,
                                       "AvailableTorque",
                                       R"doc(

         Enumeration of available torque types.

         Enumeration of torque types supported by tudat.





      )doc" )
            .value( "torque_free_type",
                    tba::AvailableTorque::torque_free,
                    R"doc(
      )doc" )
            .value( "underfined_type", tba::AvailableTorque::underfined_torque, R"doc(No documentation found.)doc" )
            .value( "second_order_gravitational_type",
                    tba::AvailableTorque::second_order_gravitational_torque,
                    R"doc(No documentation found.)doc" )
            .value( "aerodynamic_type", tba::AvailableTorque::aerodynamic_torque, R"doc(No documentation found.)doc" )
            .value( "radiation_pressure_torque_type", tba::AvailableTorque::radiation_pressure_torque, R"doc(No documentation found.)doc" )
            .value( "spherical_harmonic_gravitational_type",
                    tba::AvailableTorque::spherical_harmonic_gravitational_torque,
                    R"doc(No documentation found.)doc" )
            .value( "inertial_type", tba::AvailableTorque::inertial_torque, R"doc(No documentation found.)doc" )
            .value( "dissipative_type", tba::AvailableTorque::dissipative_torque, R"doc(No documentation found.)doc" )
            .value( "full_two_body_spherical_harmonic_gravitational_type",
                    tba::AvailableTorque::full_two_body_spherical_harmonic_gravitational_torque,
                    R"doc(No documentation found.)doc" )
            .value( "fourth_degree_full_two_body_gravitational_type",
                    tba::AvailableTorque::fourth_degree_full_two_body_gravitational_torque,
                    R"doc(No documentation found.)doc" )
            .export_values( );

    py::class_< tss::TorqueSettings, std::shared_ptr< tss::TorqueSettings > >( m,
                                                                               "TorqueSettings",
                                                                               R"doc(

         Functional base class to define settings for torques.

         This is a functional base class to define settings for torques that require no information in addition to their type.
         Classes defining settings for torque models requiring additional information must be
         derived from this class.
         Bodies exerting and undergoing torque are set outside of this class.
         This class can be used for the easy setup of torque models
         (see createTorqueModels.h), but users may also chose to do so manually.
         (Derived) Class members are all public, for ease of access and modification.





      )doc" );

    py::class_< tss::SphericalHarmonicTorqueSettings, std::shared_ptr< tss::SphericalHarmonicTorqueSettings >, tss::TorqueSettings >(
            m,
            "SphericalHarmonicTorqueSettings",
            R"doc(

         `TorqueSettings`-derived class to define settings for torques caused by spherical harmonic gravity.

         `TorqueSettings`-derived class to define settings for torques caused by spherical harmonic gravity.





      )doc" );

    py::class_< tss::FullTwoBodySphericalHarmonicTorqueSettings,
                std::shared_ptr< tss::FullTwoBodySphericalHarmonicTorqueSettings >,
                tss::TorqueSettings >( m,
                                       "FullTwoBodySphericalHarmonicTorqueSettings",
                                       R"doc(

         `TorqueSettings`-derived class to define settings for full two-body spherical harmonic gravitational torques.

         `TorqueSettings`-derived class for torques resulting from the full two-body spherical harmonic interaction,
         including figure-figure couplings between the gravity-field coefficients of both bodies, following
         :cite:t:`dirkx2019`.





      )doc" );

    py::class_< tss::FourthDegreeFullTwoBodyGravitationalTorqueSettings,
                std::shared_ptr< tss::FourthDegreeFullTwoBodyGravitationalTorqueSettings >,
                tss::TorqueSettings >( m,
                                       "FourthDegreeFullTwoBodyGravitationalTorqueSettings",
                                       R"doc(

         `TorqueSettings`-derived class to define settings for the closed-form degree-two by degree-two gravitational torque.

         `TorqueSettings`-derived class for the fourth-degree figure-figure gravitational torque generated by the
         interaction of the degree-two gravity fields of two extended bodies.





      )doc" );

    m.def( "aerodynamic",
           &tss::aerodynamicTorque,
           R"doc(

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


     )doc" );

    m.def( "radiation_pressure_torque", &tss::radiationPressureTorque, R"doc(No documentation found.)doc" );

    m.def( "second_degree_gravitational",
           &tss::secondDegreeGravitationalTorque,
           R"doc(

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


     )doc" );

    m.def( "spherical_harmonic_gravitational",
           &tss::sphericalHarmonicGravitationalTorque,
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           R"doc(

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


     )doc" );

    m.def( "full_two_body_spherical_harmonic_gravitational",
           py::overload_cast< const int, const int, const int, const int >( &tss::fullTwoBodySphericalHarmonicGravitationalTorque ),
           py::arg( "maximum_degree_body_undergoing_torque" ),
           py::arg( "maximum_order_body_undergoing_torque" ),
           py::arg( "maximum_degree_body_exerting_torque" ),
           py::arg( "maximum_order_body_exerting_torque" ),
           R"doc(

 Creates the settings for the full two-body spherical harmonic gravitational torque.

 Creates settings for the torque on one extended body due to the full two-body spherical harmonic interaction with
 another extended body. The model includes figure-figure couplings between spherical harmonic coefficients of both
 bodies, following :cite:t:`dirkx2019`. The corresponding acceleration settings use the same degree/order limits.

 The torque is obtained from the same mutual potential and effective-coefficient construction used by
 :func:`~tudatpy.dynamics.propagation_setup.acceleration.full_two_body_spherical_harmonic_gravity`.
 In the notation of :cite:t:`dirkx2019`, for body 2 expressed in the frame :math:`F_1` of body 1, the torque
 follows schematically from the angular-momentum operator applied to the potential:

 .. math::

    \mathbf{M}^{F_1}_{2} =
    -\hat{\mathcal{J}}\left(V_{1-2}\right)

 .. math::

    \mathbf{M}^{F_1}_{2} =
    -G M_1 M_2
    \sum_{l_1=0}^{\infty}\sum_{m_1=-l_1}^{l_1}
    \sum_{l_2=0}^{\infty}\sum_{m_2=-l_2}^{l_2}
    \hat{\mathcal{J}}\left(\bar{\mathcal{M}}^{2,F_1}_{l_2,m_2}\right)
    \hat{u}^{l_1,m_1}_{l_2,m_2}
    \frac{Y_{l_1+l_2,m_1+m_2}(\vartheta,\varphi)}
         {r^{l_1+l_2+1}}

 where :math:`\mathbf{M}^{F_1}_{2}` is the torque associated with the rotation of body 2 coefficients expressed
 in :math:`F_1`, :math:`\hat{\mathcal{J}}` is the angular-momentum operator, :math:`\bar{\mathcal{M}}^{2,F_1}_{l_2,m_2}`
 is the complex normalized spherical harmonic coefficient of body 2 after transformation to :math:`F_1`,
 :math:`\hat{u}^{l_1,m_1}_{l_2,m_2}` contains the remaining mass, radius, normalization, and body-1 coefficient
 factors, :math:`Y_{lm}` is the spherical harmonic basis function, and :math:`r`, :math:`\vartheta`, and
 :math:`\varphi` define the relative position. The transformed coefficient derivative is computed from the Wigner-D
 representation:

 .. math::

    \hat{\mathcal{J}}\left(\bar{\mathcal{M}}^{2,F_1}_{l_2,m_2}\right)
    =
    \sum_{k_2=-l_2}^{l_2}
    \bar{\nu}_{lmk}\hat{\mathcal{J}}\left(D^{l_2}_{m_2,k_2}\right)
    \bar{\mathcal{M}}^{2,F_2}_{l_2,k_2}

 where :math:`D^{l_2}_{m_2,k_2}` is a Wigner D-matrix entry, :math:`\bar{\nu}_{lmk}` is the associated
 normalization/sign factor, and :math:`\bar{\mathcal{M}}^{2,F_2}_{l_2,k_2}` is the original body-2 coefficient.
 See :cite:t:`dirkx2019` for the complete coefficient definitions and sign conventions.

 Algorithmically, the model:

 * uses the same selected coefficient-pair interactions as the full two-body acceleration;
 * transforms body-2 coefficients into the frame of body 1;
 * evaluates angular-momentum-operator derivatives of the transformed coefficients through the Wigner-D cache;
 * sums the selected torque contributions consistently with the effective potential terms;
 * returns the requested physical torque in the propagation setup convention.


 Parameters
 ----------
 maximum_degree_body_undergoing_torque : int
     Maximum spherical harmonic degree of the body undergoing the torque.
 maximum_order_body_undergoing_torque : int
     Maximum spherical harmonic order of the body undergoing the torque.
 maximum_degree_body_exerting_torque : int
     Maximum spherical harmonic degree of the body exerting the torque.
 maximum_order_body_exerting_torque : int
     Maximum spherical harmonic order of the body exerting the torque.
 Returns
 -------
 FullTwoBodySphericalHarmonicTorqueSettings
     Full two-body spherical harmonic torque settings object.





     )doc" );

    m.def( "full_two_body_spherical_harmonic_gravitational_from_coefficient_combinations",
           py::overload_cast< const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& >(
                   &tss::fullTwoBodySphericalHarmonicGravitationalTorque ),
           py::arg( "coefficient_combinations" ),
           R"doc(

 Creates the settings for the full two-body spherical harmonic gravitational torque from explicit coefficient combinations.

 Creates settings for the full two-body spherical harmonic gravitational torque using an explicit list of coefficient
 combinations. Each entry in ``coefficient_combinations`` is a tuple
 ``(degree_body_undergoing_torque, order_body_undergoing_torque, degree_body_exerting_torque, order_body_exerting_torque)``.
 Only the requested coefficient-pair interactions are included in the model.
 This is the explicit-combination counterpart of
 :func:`~tudatpy.dynamics.propagation_setup.torque.full_two_body_spherical_harmonic_gravitational`;
 see that function and
 :func:`~tudatpy.dynamics.propagation_setup.acceleration.full_two_body_spherical_harmonic_gravity`
 for the governing potential and torque formulation.


 Parameters
 ----------
 coefficient_combinations : list[tuple[int, int, int, int]]
     Coefficient combinations retained in the full two-body interaction.
 Returns
 -------
 FullTwoBodySphericalHarmonicTorqueSettings
     Full two-body spherical harmonic torque settings object.





     )doc" );

    m.def( "fourth_degree_full_two_body_gravitational",
           &tss::fourthDegreeFullTwoBodyGravitationalTorque,
           R"doc(

 Creates the settings for the closed-form degree-two by degree-two figure-figure gravitational torque.

 Creates settings for the fourth-degree full two-body gravitational torque generated by the interaction of the
 degree-two gravity fields of two extended bodies. This is a specialized closed-form model for the degree-two by
 degree-two figure-figure torque based on the fourth-order mutual potential and torque formulation of
 :cite:t:`schutz1981`.

 In compact tensor notation equivalent to the degree-two by degree-two part of the model, define the symmetric,
 trace-free gravity-field tensor :math:`\mathbf{K}_i` of body :math:`i` from its unnormalized degree-two
 coefficients as

 .. math::

    \mathbf{K}_i =
    \begin{bmatrix}
    C_{20}/3 - 2C_{22} & -2S_{22} & -C_{21}\\
    -2S_{22} & C_{20}/3 + 2C_{22} & -S_{21}\\
    -C_{21} & -S_{21} & -2C_{20}/3
    \end{bmatrix}_i

 With :math:`\mathbf{n}=\mathbf{r}/r`, :math:`\mathbf{r}` the position of body 2 relative to body 1 in the
 body-1-fixed frame, and body-2 quantities rotated into that same frame, define
 :math:`a=\mathbf{n}^T\mathbf{K}_1\mathbf{n}`,
 :math:`b=\mathbf{n}^T\mathbf{K}_2\mathbf{n}`,
 :math:`c=\mathbf{n}^T\mathbf{K}_1\mathbf{K}_2\mathbf{n}`, and
 :math:`d=\mathrm{tr}(\mathbf{K}_1\mathbf{K}_2)`. The degree-two by degree-two figure-figure part of the
 mutual potential can then be written as

 .. math::

    V_{22} =
    \frac{G M_1 M_2 R_1^2 R_2^2}{4r^5}
    \left(105ab-60c+6d\right)

 and the torque on body 1 can be written as

 .. math::

    \mathbf{T}_1 =
    -\frac{G M_1 M_2 R_1^2 R_2^2}{2r^5}
    \mathrm{vex}\left(\mathbf{K}_1\mathbf{H}_1-\mathbf{H}_1\mathbf{K}_1\right)

 with

 .. math::

    \mathbf{H}_1 =
    105b\mathbf{n}\mathbf{n}^T
    -30\left(\mathbf{K}_2\mathbf{n}\mathbf{n}^T+
    \mathbf{n}\mathbf{n}^T\mathbf{K}_2\right)
    +6\mathbf{K}_2

 Here :math:`\mathrm{vex}(\mathbf{S})=[S_{32},S_{13},S_{21}]^T` maps a skew-symmetric matrix to its axial vector.
 The implemented model evaluates the corresponding closed-form Schutz tensor-component expression in the
 body-1-fixed frame.


 Returns
 -------
 FourthDegreeFullTwoBodyGravitationalTorqueSettings
     Fourth-degree full two-body gravitational torque settings object.





     )doc" );

    m.def( "custom_torque",
           &tss::customTorqueSettings,
           py::arg( "torque_function" ),
           py::arg_v( "scaling_function", std::function< double( const double ) >( ), "None" ),
           R"doc(No documentation found.)doc" );

    m.def( "custom",
           &tss::customTorqueSettingsDeprecated,
           py::arg( "torque_function" ),
           py::arg_v( "scaling_function", std::function< double( const double ) >( ), "None" ) );

    // NOTE: the only unexposed torque model is
    // dissipativeTorque, but it is probably obsolete
}

}  // namespace torque
}  // namespace propagation_setup
}  // namespace dynamics
}  // namespace tudatpy
