/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/docstrings.h"

#include "expose_cr3bp.h"

#include <tudat/astro/gravitation/librationPoint.h>
#include <tudat/astro/gravitation/jacobiEnergy.h>
#include <tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h>
#include <tudat/simulation/propagation_setup/createCR3BPPeriodicOrbits.h>
#include <tudat/simulation/propagation_setup/createCR3BPManifolds.h>
#include <tudat/math/integrators/createNumericalIntegrator.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace trf = tudat::root_finders;
namespace tg = tudat::gravitation;
namespace tni = tudat::numerical_integrators;
namespace tcr3bp = tudat::circular_restricted_three_body_problem;


namespace tudatpy {
namespace astro {
namespace cr3bp {

std::function< double( const Eigen::Vector6d ) > defaultNumericalContinuationFunction =
        std::bind( &tp::getDefaultPseudoArcLength, 1.0E-4, std::placeholders::_1 );

std::vector< std::shared_ptr< tp::PropagatedCR3BPPeriodicOrbitConditions > > createCR3BPPeriodicOrbitsThroughNumericalContinuationPy(
        std::vector< std::shared_ptr< tp::PropagatedCR3BPPeriodicOrbitConditions > >& periodicOrbits,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const tp::CR3BPPeriodicOrbitGenerationSettings periodicOrbitSettings,
        const std::function< double( const Eigen::Vector6d& ) > pseudoArcLengthFunction )
{
    tp::createCR3BPPeriodicOrbitsThroughNumericalContinuation(
                periodicOrbits, integratorSettings, periodicOrbitSettings, pseudoArcLengthFunction );
    return periodicOrbits;
}



std::vector< std::vector< std::map< double, Eigen::Vector6d > > >
computeManifoldsPy( const std::shared_ptr< tp::PropagatedCR3BPPeriodicOrbitConditions > periodicOrbitConditions,
                       const double eigenvectorDisplacementFromOrbit,
                       const int numberOfDeparturePoints,
                       const double maxEigenvalueDeviation,
                       const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > > integratorSettings ,
                       const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double> > manifoldIntegratorSettings )
{
        std::vector< std::vector< std::map< double, Eigen::Vector6d > > > manifolds;
        tp::computeManifolds(
                    manifolds, periodicOrbitConditions, eigenvectorDisplacementFromOrbit, numberOfDeparturePoints,
                    maxEigenvalueDeviation, integratorSettings, manifoldIntegratorSettings );
        return manifolds;
}


void expose_cr3bp(py::module &m)
{


    py::enum_<tp::CR3BPPeriodicOrbitTypes>(m, "CR3BPPeriodicOrbitTypes", get_docstring("CR3BPPeriodicOrbitTypes").c_str())
            .value("horizontal_lyapunov_orbit", tp::CR3BPPeriodicOrbitTypes::horizontal_lyapunov_orbit)
            .value("vertical_lyapunov_orbit", tp::CR3BPPeriodicOrbitTypes::vertical_lyapunov_orbit)
            .value("halo_orbit", tp::CR3BPPeriodicOrbitTypes::halo_orbit)
            .export_values();

    py::class_<tp::CR3BPPeriodicOrbitGenerationSettings,
            std::shared_ptr<tp::CR3BPPeriodicOrbitGenerationSettings> >(
                m, "CR3BPPeriodicOrbitGenerationSettings",
                get_docstring("CR3BPPeriodicOrbitGenerationSettings").c_str())
            .def(py::init<
                 const double,
                 const tp::CR3BPPeriodicOrbitTypes,
                 const int,
                 const int,
                 const double,
                 const double ,
                 const int,
                 const double>(),
                 py::arg("mass_parameter"),
                 py::arg("orbit_type"),
                 py::arg("libration_point_number"),
                 py::arg("maximum_number_differential_corrections"),
                 py::arg("maximum_position_deviation"),
                 py::arg("maximum_velocity_deviation"),
                 py::arg("maximum_numerical_continuation") = -1,
                 py::arg("maximum_eigenvalue_deviation") = TUDAT_NAN );


    py::class_<tp::CR3BPPeriodicOrbitConditions,
            std::shared_ptr<tp::CR3BPPeriodicOrbitConditions> >(
                m, "CR3BPPeriodicOrbitConditions",
                get_docstring("CR3BPPeriodicOrbitConditions").c_str())
            .def_property_readonly( "initial_state", &tp::CR3BPPeriodicOrbitConditions::getInitialState )
            .def_property_readonly( "orbital_period", &tp::CR3BPPeriodicOrbitConditions::getOrbitPeriod )
            .def_property_readonly( "mass_parameter", &tp::CR3BPPeriodicOrbitConditions::getMassParameter );

    py::class_<tp::PropagatedCR3BPPeriodicOrbitConditions,
            std::shared_ptr<tp::PropagatedCR3BPPeriodicOrbitConditions>,
            tp::CR3BPPeriodicOrbitConditions >(
                m, "PropagatedCR3BPPeriodicOrbitConditions",
                get_docstring("PropagatedCR3BPPeriodicOrbitConditions").c_str())
            .def_property_readonly( "monodromy_matrix", &tp::PropagatedCR3BPPeriodicOrbitConditions::getMonodromyMatrix )
            .def_property_readonly( "orbit_type", &tp::PropagatedCR3BPPeriodicOrbitConditions::getOrbitType )
            .def_property_readonly( "libration_point_number", &tp::PropagatedCR3BPPeriodicOrbitConditions::getLibrationPointNumber );

    py::class_<tp::GeneratedCR3BPPeriodicOrbitConditions,
            std::shared_ptr<tp::GeneratedCR3BPPeriodicOrbitConditions>,
            tp::PropagatedCR3BPPeriodicOrbitConditions >(
                m, "GeneratedCR3BPPeriodicOrbitConditions",
                get_docstring("GeneratedCR3BPPeriodicOrbitConditions").c_str());

    m.def("dimensionless_to_regular_time",
          &tcr3bp::convertDimensionlessTimeToDimensionalTime,
          py::arg( "dimensionless_time" ),
          py::arg( "primary_gravitational_parameter" ),
          py::arg( "secondary_gravitational_parameter" ),
          py::arg( "primary_secondary_distance" ),
          get_docstring("corotating_normalized_to_inertial_state").c_str() );

    m.def("regular_to_dimensionless_time",
          &tcr3bp::convertDimensionalTimeToDimensionlessTime,
          py::arg( "regular_time" ),
          py::arg( "primary_gravitational_parameter" ),
          py::arg( "secondary_gravitational_parameter" ),
          py::arg( "primary_secondary_distance" ),
          get_docstring("corotating_normalized_to_inertial_state").c_str() );

    m.def("corotating_normalized_to_inertial_state",
          &tcr3bp::convertCorotatingNormalizedToCartesianCoordinates,
          py::arg( "primary_gravitational_parameter" ),
          py::arg( "secondary_gravitational_parameter" ),
          py::arg( "primary_secondary_distance" ),
          py::arg( "normalized_corotating_state" ),
          py::arg( "regular_time" ),
          get_docstring("corotating_normalized_to_inertial_state").c_str() );

    m.def("inertial_to_corotating_normalized_state",
          &tcr3bp::convertCartesianToCorotatingNormalizedCoordinates,
          py::arg( "primary_gravitational_parameter" ),
          py::arg( "secondary_gravitational_parameter" ),
          py::arg( "primary_secondary_distance" ),
          py::arg( "inertial_state" ),
          py::arg( "time" ),
          get_docstring("inertial_to_corotating_normalized_state").c_str() );

    m.def("compute_jacobi_energy",
          &tg::computeJacobiEnergy,
          py::arg( "mass_parameter" ),
          py::arg( "state" ),
          get_docstring("computeJacobiEnergy").c_str() );

    m.def("compute_mass_parameter",
          &tcr3bp::computeMassParameter,
          py::arg( "primary_gravitational_parameter" ),
          py::arg( "secondary_gravitational_parameter" ),
          get_docstring("computeJacobiEnergy").c_str() );

    m.def("compute_libration_point_position",
          &tcr3bp::computeLibrationPointPosition,
          py::arg( "mass_parameter" ),
          py::arg( "libration_point_number" ),
          py::arg( "root_finder" )= std::make_shared< trf::NewtonRaphson< > >( 1.0e-14, 1000 ),
          get_docstring("compute_libration_point_position").c_str() );


    m.def("compute_cr3bp_state_derivative",
          &tp::computeCr3bpStateDerivative,
          py::arg("dimensionless_time"),
          py::arg("normalized_corotating_state"),
          py::arg("mass_parameter"),
          get_docstring("compute_cr3bp_state_derivative").c_str() );

   m.def("get_orbit_jacobi_energy",
         &tp::getOrbitJacobiEnergyHistory,
         py::arg( "orbit_state_history" ),
         py::arg( "mass_parameter" ),
         get_docstring("get_orbit_jacobi_energy").c_str() );

    m.def("propagate_periodic_orbit",
          &tp::propagatePeriodicOrbit,
          py::arg( "periodic_orbit_definition" ),
          py::arg( "integrator_settings" ),
          py::arg( "number_of_periods" ) = 1.0,
          get_docstring("propagate_periodic_orbit").c_str() );

    m.def("propagate_cr3bp",
          py::overload_cast<
          const std::shared_ptr< tni::IntegratorSettings< double > >,
          const double,
          const Eigen::Vector6d&,
          const double,
          const bool,
          const bool
          >( &tp::performCR3BPIntegration ),
          py::arg( "integrator_settings" ),
          py::arg( "mass_parameter" ),
          py::arg( "initial_normalized_corotating_state" ),
          py::arg( "final_time" ),
          py::arg( "propagate_to_exact_final_time" ) = false,
          py::arg( "save_full_results" ) = true,
          get_docstring("propagated_cr3bp").c_str() );

    m.def("propagated_cr3bp_with_state_transition",
          py::overload_cast<
          const std::shared_ptr< tni::IntegratorSettings< double > >,
          const double,
          const Eigen::Vector6d&,
          const double,
          const bool,
          const bool
          >( &tp::performCR3BPWithStmIntegration ),
          py::arg( "integrator_settings" ),
          py::arg( "mass_parameter" ),
          py::arg( "initial_normalized_corotating_state" ),
          py::arg( "final_time" ),
          py::arg( "propagate_to_exact_final_time" ) = false,
          py::arg( "save_full_results" ) = true,
          get_docstring("propagated_cr3bp").c_str() );

    m.def("default_numerical_continuaion_pseuco_arc_length",
          &tp::getDefaultPseudoArcLength,
          py::arg( "distance_increment" ),
          py::arg( "current_state" ),
          get_docstring("default_numerical_continuaion_pseuco_arc_length").c_str() );

    m.def("numerically_continue_libration_point_periodic_orbits",
          py::overload_cast<
          std::vector< std::shared_ptr< tp::PropagatedCR3BPPeriodicOrbitConditions > >&,
          const std::shared_ptr< tni::IntegratorSettings< double > >,
          const tp::CR3BPPeriodicOrbitGenerationSettings,
          const std::function< double( const Eigen::Vector6d& ) > >(
              &createCR3BPPeriodicOrbitsThroughNumericalContinuationPy ),
          py::arg( "periodic_orbits" ),
          py::arg( "integrator_settings" ),
          py::arg( "periodic_orbit_settings" ),
          py::arg( "pseudo_arc_length_function" ) = &defaultNumericalContinuationFunction,
          get_docstring("numerically_continue_libration_point_periodic_orbits").c_str() );

    m.def("create_earth_moon_periodic_orbit_initial_guess",
          &tp::richardsonApproximationEarthMoonLibrationPointPeriodicOrbit,
          py::arg( "mass_parameter" ),
          py::arg( "orbit_type" ),
          py::arg( "libration_point_number" ),
          py::arg( "guess_iteration" ),
          py::arg( "n" ) = 1.0,
          get_docstring("create_earth_moon_periodic_orbit_initial_guess").c_str() );

    m.def("create_libration_point_periodic_orbit",
          &tp::createCR3BPPeriodicOrbit,
          py::arg( "initial_normalized_corotating_state_guess" ),
          py::arg( "orbital_period_guess" ),
          py::arg( "periodic_orbit_settings" ),
          py::arg( "integrator_settings" ),
          get_docstring("create_libration_point_periodic_orbit").c_str() );

    m.def("create_manifolds",
          &computeManifoldsPy,
          py::arg( "periodic_orbit" ),
          py::arg( "initial_eigenvector_displacement" ),
          py::arg( "number_of_departure_points" ),
          py::arg( "maximum_eigenvalue_deviation" ),
          py::arg( "orbit_integrator_settings" ),
          py::arg( "manifold_integrator_settings" ),
          get_docstring("create_manifolds").c_str() );

}
} // namespace cr3bp
} // namespace astro
} // namespace tudatpy
