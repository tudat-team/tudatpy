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

std::shared_ptr< tr::PPNParameterSet > getGlobalPpnParameterSet( )
{
    return tr::ppnParameterSet;
}

void setGlobalPpnParameterSet( const std::shared_ptr< tr::PPNParameterSet >& ppn_parameter_set )
{
    tr::ppnParameterSet = ppn_parameter_set;
}

double getEquivalencePrincipleLpiViolationParameter( )
{
    return tr::equivalencePrincipleLpiViolationParameter;
}

void setEquivalencePrincipleLpiViolationParameter( const double value )
{
    tr::equivalencePrincipleLpiViolationParameter = value;
}

void createBaseMetric(
        const std::shared_ptr< tss::SpaceTimeMetricSettings >& space_time_metric_settings,
        const tss::SystemOfBodies& bodies )
{
    tss::createBaseMetric( space_time_metric_settings, bodies );
}

void expose_space_time_setup( py::module& m )
{
    py::class_< tr::PPNParameterSet, std::shared_ptr< tr::PPNParameterSet > >(
            m, "PPNParameterSet", R"doc(

        Container class for PPN parameters used in relativistic models.

        The first-order PPN parameters :math:`\gamma` and :math:`\beta` follow
        the conventions summarized by :cite:t:`willlrr2014`.
        The parameter :math:`\epsilon` is used here as a second-order
        post-Newtonian parameter.

     )doc" )
            .def( py::init( []( const double parameter_gamma,
                                const double parameter_beta,
                                const double parameter_epsilon )
                            {
                                return std::make_shared< tr::PPNParameterSet >(
                                            parameter_gamma, parameter_beta, 0.0, parameter_epsilon );
                            } ),
                  py::arg( "parameter_gamma" ),
                  py::arg( "parameter_beta" ),
                  py::arg( "parameter_epsilon" ) = 0.0 )
            .def_property(
                    "parameter_gamma",
                    &tr::PPNParameterSet::getParameterGamma,
                    &tr::PPNParameterSet::setParameterGamma,
                    R"doc(PPN parameter :math:`\gamma`.)doc" )
            .def_property(
                    "parameter_beta",
                    &tr::PPNParameterSet::getParameterBeta,
                    &tr::PPNParameterSet::setParameterBeta,
                    R"doc(PPN parameter :math:`\beta`.)doc" )
            .def_property(
                    "parameter_epsilon",
                    &tr::PPNParameterSet::getParameterEpsilon,
                    &tr::PPNParameterSet::setParameterEpsilon,
                    R"doc(Second-order post-Newtonian parameter :math:`\epsilon`.)doc" );

    py::enum_< tss::SpaceTimeMetricTypes >( m, "SpaceTimeMetricType" )
            .value( "schwardschild_metric", tss::schwardschild_metric )
            .value( "solar_system_metric", tss::solar_system_metric );

    py::class_< tss::SpaceTimeMetricSettings, std::shared_ptr< tss::SpaceTimeMetricSettings > >(
            m, "SpaceTimeMetricSettings", R"doc(

        Base class for space-time metric settings.

     )doc" );

    py::class_< tss::SchwardschildSpaceTimeMetricSettings,
                std::shared_ptr< tss::SchwardschildSpaceTimeMetricSettings >,
                tss::SpaceTimeMetricSettings >(
            m, "SchwardschildSpaceTimeMetricSettings", R"doc(

        Settings for a harmonic Schwarzschild metric.

     )doc" );

    py::class_< tss::SolarSystemSpaceTimeMetricSettings,
                std::shared_ptr< tss::SolarSystemSpaceTimeMetricSettings >,
                tss::SpaceTimeMetricSettings >(
            m, "SolarSystemSpaceTimeMetricSettings", R"doc(

        Settings for a solar-system post-Newtonian metric.

     )doc" );

    m.def(
            "schwarzschild_metric_settings",
            &tss::schwardschildSpaceTimeMetricSettings,
            py::arg( "body" ),
            py::arg( "include_second_post_newtonian_order" ) = false,
            py::arg( "ppn_parameter_set" ) = tr::ppnParameterSet,
            R"doc(

 Create settings for a harmonic Schwarzschild metric.

 The implemented covariant metric perturbation is

 .. math::

     h_{00}=2\frac{U}{c^2}-2\beta\frac{U^2}{c^4},\qquad
     h_{ij}=2\gamma\frac{U}{c^2}\delta_{ij}
     +2\epsilon\frac{U^2}{c^4}\delta_{ij},

 with :math:`h_{0i}=0` for this metric model.

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
 ppn_parameter_set : PPNParameterSet, optional
     PPN parameter set used in the metric. Defaults to the global PPN parameter set.

 Returns
 -------
 SchwardschildSpaceTimeMetricSettings
     Settings object for the Schwarzschild metric model.

        )doc" );

    m.def(
            "solar_system_metric_settings",
            &tss::solarSystemSpaceTimeMetricSettings,
            py::arg( "first_order_bodies" ),
            py::arg( "second_order_bodies" ) = std::vector< std::string >( ),
            py::arg( "spherical_harmonic_expansions" ) = std::map< std::string, std::pair< int, int > >( ),
            py::arg( "angular_momentum_bodies" ) = std::vector< std::string >( ),
            py::arg( "use_body_accelerations" ) = true,
            py::arg( "ppn_parameter_set" ) = tr::ppnParameterSet,
            R"doc(

 Create settings for a solar-system post-Newtonian metric.

 The implemented metric uses the standard post-Newtonian decomposition
 with scalar potential :math:`w` and vector potential :math:`\mathbf{w}`:

 .. math::

     h_{00}=2\frac{w}{c^2}-2\beta\frac{w^2}{c^4},\qquad
     h_{0i}=-2(\gamma+1)\frac{w_i}{c^3},

 .. math::

     h_{ij}=2\gamma\frac{w}{c^2}\delta_{ij}
     +2(\gamma^2+\beta-1)\frac{w^2_{\mathrm{(2)}}}{c^4}\delta_{ij}
     +2(\gamma+1)\frac{q_{ij}}{c^4},

 where :math:`w` is assembled from the configured first-order bodies,
 :math:`w^2_{\mathrm{(2)}}` contains squared-potential terms from bodies
 configured as second-order, :math:`q_{ij}` is the anisotropic matrix potential,
 :math:`\mathbf{w}` is the vector potential, :math:`w_i` its :math:`i`-th component,
 :math:`\delta_{ij}` is the Kronecker delta, :math:`c` is the speed of light,
 and :math:`\beta,\gamma` are the PPN parameters.

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
 ppn_parameter_set : PPNParameterSet, optional
     PPN parameter set used in the metric. Defaults to the global PPN parameter set.

 Returns
 -------
 SolarSystemSpaceTimeMetricSettings
     Settings object for the solar-system metric model.

        )doc" );

    m.def(
            "create_base_metric",
            &createBaseMetric,
            py::arg( "metric_settings" ),
            py::arg( "bodies" ),
            R"doc(

 Create and assign the global base metric from metric settings and bodies.

 This function creates the configured metric model and assigns it to Tudat's
 global ``baseMetric`` variable used by direct-from-metric relativistic-time propagation.

 Parameters
 ----------
 metric_settings : SpaceTimeMetricSettings
     Metric settings defining the metric model to create.
 bodies : SystemOfBodies
     System of bodies providing required environment models and state functions.

 Returns
 -------
 None
     The global base metric is updated in place.

        )doc" );

    m.def(
            "get_global_ppn_parameter_set",
            &getGlobalPpnParameterSet,
            R"doc(

 Retrieve Tudat's global PPN parameter set.

 This is the global PPN parameter set used by all Tudat functionality that
 requires PPN parameters, unless explicitly overridden in local settings.

 Returns
 -------
 PPNParameterSet
     Global PPN parameter set used by default in relativistic models.

        )doc" );

    m.def(
            "set_global_ppn_parameter_set",
            &setGlobalPpnParameterSet,
            py::arg( "ppn_parameter_set" ),
            R"doc(

 Replace Tudat's global PPN parameter set.

 The assigned value is returned by
 :func:`~get_global_ppn_parameter_set`.

 Parameters
 ----------
 ppn_parameter_set : PPNParameterSet
     New global PPN parameter set.

 Returns
 -------
 None

        )doc" );

    m.def(
            "get_equivalence_principle_lpi_violation_parameter",
            &getEquivalencePrincipleLpiViolationParameter,
            R"doc(

 Retrieve the global equivalence-principle LPI-violation parameter.

 This is the global value used by all Tudat functionality that requires the
 equivalence-principle local-position-invariance violation parameter.

 Returns
 -------
 float
     Global LPI-violation parameter value.

        )doc" );

    m.def(
            "set_equivalence_principle_lpi_violation_parameter",
            &setEquivalencePrincipleLpiViolationParameter,
            py::arg( "value" ),
            R"doc(

 Set the global equivalence-principle LPI-violation parameter.

 The assigned value is returned by
 :func:`~get_equivalence_principle_lpi_violation_parameter`.

 Parameters
 ----------
 value : float
     New value for the global LPI-violation parameter.

 Returns
 -------
 None

        )doc" );
}

}  // namespace space_time
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
