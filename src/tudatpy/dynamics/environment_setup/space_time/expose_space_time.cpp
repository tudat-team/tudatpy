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
#include "expose_space_time.h"

#include <pybind11/stl.h>
#include <tudat/astro/relativity/metric.h>
#include <tudat/simulation/environment_setup/createBodiesSettings.h>
#include <tudat/simulation/environment_setup/createMetric.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tr = tudat::relativity;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace space_time
{

void expose_space_time_setup( py::module& m )
{
    py::class_< tr::PPNParameterSet, std::shared_ptr< tr::PPNParameterSet > >( m, "PPNParameterSet", R"doc(

        Container class for PPN parameters used in relativistic models.
        These parameters are assigned through
        :attr:`~tudatpy.dynamics.environment.SystemOfBodies.space_time_properties`.

        The first-order PPN parameters :math:`\gamma` and :math:`\beta` follow
        the conventions summarized by :cite:t:`willlrr2014`.
        The parameter :math:`\epsilon` is used here as a second-order
        post-Newtonian parameter.

     )doc" )
            .def( py::init( []( const double parameter_gamma, const double parameter_beta, const double parameter_epsilon ) {
                      return std::make_shared< tr::PPNParameterSet >( parameter_gamma, parameter_beta, 0.0, parameter_epsilon );
                  } ),
                  py::arg( "parameter_gamma" ),
                  py::arg( "parameter_beta" ),
                  py::arg( "parameter_epsilon" ) = 0.0 )
            .def_property( "parameter_gamma",
                           &tr::PPNParameterSet::getParameterGamma,
                           &tr::PPNParameterSet::setParameterGamma,
                           R"doc(PPN parameter :math:`\gamma`.)doc" )
            .def_property( "parameter_beta",
                           &tr::PPNParameterSet::getParameterBeta,
                           &tr::PPNParameterSet::setParameterBeta,
                           R"doc(PPN parameter :math:`\beta`.)doc" )
            .def_property( "parameter_epsilon",
                           &tr::PPNParameterSet::getParameterEpsilon,
                           &tr::PPNParameterSet::setParameterEpsilon,
                           R"doc(Second-order post-Newtonian parameter :math:`\epsilon`.)doc" );

    py::enum_< tss::SpaceTimeMetricTypes >( m, "SpaceTimeMetricType" )
            .value( "schwarzschild_metric", tss::schwarzschild_metric )
            .value( "solar_system_metric", tss::solar_system_metric );

    py::class_< tss::SpaceTimeMetricSettings, std::shared_ptr< tss::SpaceTimeMetricSettings > >( m, "SpaceTimeMetricSettings", R"doc(

        Base class for space-time metric settings.

     )doc" );

    py::class_< tss::SchwarzschildSpaceTimeMetricSettings,
                std::shared_ptr< tss::SchwarzschildSpaceTimeMetricSettings >,
                tss::SpaceTimeMetricSettings >( m, "SchwarzschildSpaceTimeMetricSettings", R"doc(

        Settings for a harmonic Schwarzschild metric.

     )doc" );

    py::class_< tss::SolarSystemSpaceTimeMetricSettings,
                std::shared_ptr< tss::SolarSystemSpaceTimeMetricSettings >,
                tss::SpaceTimeMetricSettings >( m, "SolarSystemSpaceTimeMetricSettings", R"doc(

        Settings for a solar-system post-Newtonian metric.

     )doc" );

    py::class_< tss::SpaceTimePropertiesSettings, std::shared_ptr< tss::SpaceTimePropertiesSettings > >(
            m, "SpaceTimePropertiesSettings", R"doc(

        Settings assigned to
        :attr:`~tudatpy.dynamics.environment_setup.BodyListSettings.space_time_settings`,
        for initializing
        :attr:`~tudatpy.dynamics.environment.SystemOfBodies.space_time_properties`.

     )doc" )
            .def( py::init< const std::shared_ptr< tss::SpaceTimeMetricSettings >&,
                            const std::shared_ptr< tr::PPNParameterSet >&,
                            const double >( ),
                  py::arg_v( "metric_settings", std::shared_ptr< tss::SpaceTimeMetricSettings >( ), "None" ).none( true ),
                  py::arg_v( "ppn_parameter_set", std::shared_ptr< tr::PPNParameterSet >( ), "None" ).none( true ),
                  py::arg( "equivalence_principle_lpi_violation_parameter" ) = 0.0 )
            .def_property( "metric_settings",
                           &tss::SpaceTimePropertiesSettings::getMetricSettings,
                           &tss::SpaceTimePropertiesSettings::setMetricSettings )
            .def_property( "ppn_parameter_set",
                           &tss::SpaceTimePropertiesSettings::getPpnParameterSet,
                           &tss::SpaceTimePropertiesSettings::setPpnParameterSet )
            .def_property( "equivalence_principle_lpi_violation_parameter",
                           &tss::SpaceTimePropertiesSettings::getEquivalencePrincipleLpiViolationParameter,
                           &tss::SpaceTimePropertiesSettings::setEquivalencePrincipleLpiViolationParameter );

    m.def( "ppn_parameter_set",
           &tss::ppnParameterSet,
           py::arg( "parameter_gamma" ) = 1.0,
           py::arg( "parameter_beta" ) = 1.0,
           py::arg( "parameter_epsilon" ) = 0.0,
           R"doc(

 Create settings for PPN parameters used by system-level space-time properties.

 The returned object is typically assigned to
 :attr:`~tudatpy.dynamics.environment_setup.space_time.SpaceTimePropertiesSettings.ppn_parameter_set`, which is then used to build
 :attr:`~tudatpy.dynamics.environment.SystemOfBodies.space_time_properties`.

 Parameters
 ----------
 parameter_gamma : float, optional
     First-order PPN parameter :math:`\gamma`.
 parameter_beta : float, optional
     First-order PPN parameter :math:`\beta`.
 parameter_epsilon : float, optional
     Second-order post-Newtonian parameter :math:`\epsilon`.

 Returns
 -------
 PPNParameterSet
     Settings object for PPN parameters.

        )doc" );

    m.def( "space_time_properties_settings",
           &tss::spaceTimePropertiesSettings,
           py::arg_v( "metric_settings", std::shared_ptr< tss::SpaceTimeMetricSettings >( ), "None" ).none( true ),
           py::arg_v( "ppn_parameter_set", std::shared_ptr< tr::PPNParameterSet >( ), "None" ).none( true ),
           py::arg( "equivalence_principle_lpi_violation_parameter" ) = 0.0,
           R"doc(

 Create settings for system-level space-time properties.

 This function creates a :class:`~SpaceTimePropertiesSettings` object that is assigned to
 :attr:`~tudatpy.dynamics.environment_setup.BodyListSettings.space_time_settings`. During
 :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`, these settings define
 :attr:`~tudatpy.dynamics.environment.SystemOfBodies.space_time_properties`.

 Parameters
 ----------
 metric_settings : SpaceTimeMetricSettings, optional
     Settings that define the base metric model to be used.
 ppn_parameter_set : PPNParameterSet, optional
     PPN parameter settings used by relativistic models.
 equivalence_principle_lpi_violation_parameter : float, optional
     Local-position-invariance violation parameter for equivalence-principle models.

 Returns
 -------
 SpaceTimePropertiesSettings
     Settings object used to initialize ``SystemOfBodies.space_time_properties``.

        )doc" );

    m.def( "schwarzschild_metric_settings",
           &tss::schwarzschildSpaceTimeMetricSettings,
           py::arg( "body" ),
           py::arg( "include_second_post_newtonian_order" ) = false,
           R"doc(

Create settings for a harmonic Schwarzschild metric.

The implemented covariant metric perturbation is

.. math::

   h_{00} = 2\frac{U}{c^2} - 2\beta\frac{U^2}{c^4}

.. math::

   h_{ij} = 2\gamma\frac{U}{c^2}\delta_{ij}
   + 2\epsilon\frac{U^2}{c^4}\delta_{ij}

.. math::

   h_{0i} = 0

where :math:`U=\mu/r` is the central-body potential, :math:`\mu=GM`,
:math:`r` is the distance from the central body, :math:`c` is the speed of light,
:math:`\delta_{ij}` is the Kronecker delta, and :math:`\beta,\gamma,\epsilon`
are the PPN parameters.

Parameters
----------
body : str
    Name of the central gravitating body used in the metric.
include_second_post_newtonian_order : bool, optional
    If true, include second post-Newtonian terms.

Returns
-------
SchwarzschildSpaceTimeMetricSettings
     Settings object for the Schwarzschild metric model. This can be assigned
     to :attr:`SpaceTimePropertiesSettings.metric_settings`.

        )doc" );

    m.def( "solar_system_metric_settings",
           &tss::solarSystemSpaceTimeMetricSettings,
           py::arg( "first_order_bodies" ),
           py::arg( "second_order_bodies" ) = std::vector< std::string >( ),
           py::arg( "spherical_harmonic_expansions" ) = std::map< std::string, std::pair< int, int > >( ),
           py::arg( "angular_momentum_bodies" ) = std::vector< std::string >( ),
           py::arg( "use_body_accelerations" ) = true,
           R"doc(

Create settings for a solar-system post-Newtonian metric.

The implemented metric uses the standard post-Newtonian decomposition
with scalar potential :math:`w` and vector potential :math:`\mathbf{w}`.
The formulation follows the IAU-relativistic framework conventions
summarized by :cite:t:`soffel2003`.

.. math::

    h_{00} = 2\frac{w}{c^2} - 2\beta\frac{w^2}{c^4}

.. math::

    h_{0i} = -2(\gamma + 1)\frac{w_i}{c^3}

.. math::

    h_{ij} = 2\gamma\frac{w}{c^2}\delta_{ij}
    + 2(\gamma^2 + \beta - 1)\frac{\left(w_{\mathrm{(2)}}\right)^2}{c^4}\delta_{ij}
    + 2(\gamma + 1)\frac{q_{ij}}{c^4}

In Tudat, these quantities are evaluated as

.. math::

    w = \sum_{a\in\mathcal{F}} w_a

.. math::

    w_a =
    \begin{cases}
    \mu_a/R_a, & \text{point-mass term}, \\
    U_a^{\mathrm{SH}}, & \text{if spherical-harmonic expansion is configured},
    \end{cases}

.. math::

    \mathbf{w} = \sum_{a\in\mathcal{F}} w_a\,\mathbf{v}_a
    - \sum_{k\in\mathcal{A}}
    \frac{G}{2R_k^3}\left(\mathbf{S}_k\times\mathbf{r}_k\right)

.. math::

    w_i = (\mathbf{w})_i

.. math::

    \left(w_{\mathrm{(2)}}\right)^2 = \sum_{b\in\mathcal{S}} \left(w_b\right)^2

.. math::

    q_{ij} = \sum_{b\in\mathcal{S}} q^{(b)}_{ij}

.. math::

    q^{(b)}_{ij} =
    \frac{\mu_b^2}{4R_b^2}
    \left(\frac{r_{b,i}r_{b,j}}{R_b^2} - \delta_{ij}\right)

where :math:`\mathcal{F}` is the set of bodies in ``first_order_bodies``,
:math:`\mathcal{S}` is the set of bodies in ``second_order_bodies``,
:math:`\mathcal{A}` is the set of bodies in ``angular_momentum_bodies``,
:math:`\mu_a=GM_a` is the gravitational parameter of body :math:`a`,
:math:`\mathbf{r}_a` is the relative position vector from body :math:`a` to the metric evaluation point,
:math:`R_a=\lVert\mathbf{r}_a\rVert`, :math:`\mathbf{v}_a` is the barycentric velocity of body :math:`a`,
:math:`\mathbf{S}_k` is the angular-momentum vector of body :math:`k`,
:math:`U_a^{\mathrm{SH}}` is the configured spherical-harmonic scalar potential contribution,
:math:`w_b` follows the same definition as :math:`w_a` but evaluated for body :math:`b`,
:math:`\delta_{ij}` is the Kronecker delta, :math:`c` is the speed of light,
:math:`G` is the gravitational constant, and :math:`\beta,\gamma` are the PPN parameters.

Parameters
----------
first_order_bodies : list[str]
    Bodies included in the first-order post-Newtonian terms.
second_order_bodies : list[str], optional
    Bodies included in second-order terms.
spherical_harmonic_expansions : dict[str, tuple[int, int]], optional
    Optional spherical-harmonic degree/order settings per body.
angular_momentum_bodies : list[str], optional
    Bodies for which angular-momentum terms are included.
use_body_accelerations : bool, optional
    Whether body-acceleration terms are included in the metric model.

Returns
-------
SolarSystemSpaceTimeMetricSettings
     Settings object for the solar-system metric model. This can be assigned
     to :attr:`SpaceTimePropertiesSettings.metric_settings`.

        )doc" );
}

}  // namespace space_time
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
