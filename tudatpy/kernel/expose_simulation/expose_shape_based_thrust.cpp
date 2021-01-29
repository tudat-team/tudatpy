/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_shape_based_thrust.h"

#include <tudat/astro/low_thrust/shape_based/hodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/getRecommendedBaseFunctionsHodographicShaping.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tsbm = tudat::shape_based_methods;
namespace tltt = tudat::low_thrust_trajectories;
namespace tni = tudat::numerical_integrators;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace shape_based_methods
{

inline std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getHodographicShapingAccelerationSettings(
        const Eigen::Vector6d& initialState,
        const Eigen::Vector6d& finalState,
        const double timeOfFlight,
        const double centralBodyGravitationalParameter,
        const int numberOfRevolutions,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        const Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsAxialVelocityFunction,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const double specificImpulse,
        const double timeOffset )
{
    HodographicShaping* shapingObject =
            new HodographicShaping(
                initialState, finalState, timeOfFlight, centralBodyGravitationalParameter,
                numberOfRevolutions, radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction );
    return shapingObject->getLowThrustAccelerationSettings( bodies, bodyToPropagate, specificImpulse, timeOffset );
}

inline Eigen::Vector6d getHodographicShapingStateAtEpoch(
        const Eigen::Vector6d& initialState,
        const Eigen::Vector6d& finalState,
        const double timeOfFlight,
        const double centralBodyGravitationalParameter,
        const int numberOfRevolutions,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        const Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsAxialVelocityFunction,
        const double epoch )
{
    return HodographicShaping(
                initialState, finalState, timeOfFlight, centralBodyGravitationalParameter,
                numberOfRevolutions, radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction ).
            getStateAtEpoch( epoch );
}

inline double getHodographicShapingDeltaV(
        const Eigen::Vector6d& initialState,
        const Eigen::Vector6d& finalState,
        const double timeOfFlight,
        const double centralBodyGravitationalParameter,
        const int numberOfRevolutions,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        const Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsAxialVelocityFunction )
{
    return HodographicShaping(
                initialState, finalState, timeOfFlight, centralBodyGravitationalParameter,
                numberOfRevolutions, radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction ).computeDeltaV( );
}




} // namespace shape_based_methods
} // namespace tudat

namespace tudatpy {


void expose_shape_based_thrust(py::module &m)
{
    py::class_<
            tltt::LowThrustLeg,
            std::shared_ptr<tltt::LowThrustLeg> >(m, "LowThrustLeg")
            .def( "get_thrust_acceleration_profile",
                  py::overload_cast<
                  std::vector< double >&,
                  std::function< double ( const double ) >,
                  std::shared_ptr<tni::IntegratorSettings< double > > >( &tltt::LowThrustLeg::getThrustAccelerationProfile ),
                  py::arg("output_times"),
                  py::arg("specific_impulse_function"),
                  py::arg("integrator_settings" ) = nullptr )
            .def( "get_trajectory",
                  py::overload_cast<
                  std::vector< double >& >( &tltt::LowThrustLeg::getTrajectory ),
                  py::arg("times") )
            .def( "get_state",
                  &tltt::LowThrustLeg::getStateAtEpoch,
                  py::arg("time") )
            .def( "get_low_thrust_acceleration_settings",
                  py::overload_cast<
                          const tss::SystemOfBodies&,
                          const std::string&,
                          const double,
                          const double >(
                      &tltt::LowThrustLeg::getLowThrustAccelerationSettings ),
                  py::arg("bodies"),
                  py::arg("body_to_propagate"),
                  py::arg("specific_impulse_function"),
                  py::arg("time_offset" ) )
            .def( "compute_delta_v",
                  &tltt::LowThrustLeg::computeDeltaV );

    py::class_<
            tsbm::ShapeBasedMethod,
            std::shared_ptr<tsbm::ShapeBasedMethod>,
            tltt::LowThrustLeg
            >(m, "ShapeBasedMethod");

    py::class_<
            tsbm::HodographicShaping,
            std::shared_ptr<tsbm::HodographicShaping>,
            tsbm::ShapeBasedMethod
            >(m, "HodographicShaping")
            .def(py::init<
                 const Eigen::Vector6d&,
                 const Eigen::Vector6d&,
                 const double,
                 const double,
                 const int,
                 const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
                 const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
                 const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
                 const Eigen::VectorXd&,
                 const Eigen::VectorXd&,
                 const Eigen::VectorXd&,
                 const double >(),
                 py::arg("initial_state"),
                 py::arg("final_state"),
                 py::arg("time_of_flight"),
                 py::arg("central_body_gravitational_parameter"),
                 py::arg("number_of_revolutions"),
                 py::arg("radial_velocity_functions"),
                 py::arg("normal_velocity_functions"),
                 py::arg("axial_velocity_functions"),
                 py::arg("radial_free_coefficients"),
                 py::arg("normal_free_coefficients"),
                 py::arg("axial_free_coefficients"),
                 py::arg("initial_mass") = TUDAT_NAN
            );

    py::class_<
            tsbm::BaseFunctionHodographicShaping,
            std::shared_ptr<tsbm::BaseFunctionHodographicShaping>
            >(m, "BaseFunctionHodographicShaping");

    m.def("recommended_radial_hodograph_functions",
          py::overload_cast< const double >(
              &tsbm::getRecommendedRadialVelocityBaseFunctions ),
          py::arg("time_of_flight") );

    m.def("recommended_normal_hodograph_functions",
          py::overload_cast< const double >(
              &tsbm::getRecommendedNormalBaseFunctions ),
          py::arg("time_of_flight") );


    m.def("recommended_axial_hodograph_functions",
          py::overload_cast< const double, const int >(
              &tsbm::getRecommendedAxialVelocityBaseFunctions ),
          py::arg("time_of_flight"),
          py::arg("number_of_revolutions") );


    m.def("hodograph_constant",
          &tsbm::hodographConstant );

    m.def("hodograph_sine",
          &tsbm::hodographSine,
          py::arg("frequency") );

    m.def("hodograph_cosine",
          &tsbm::hodographCosine,
          py::arg("frequency") );

    m.def("hodograph_exponential",
          &tsbm::hodographExponential,
          py::arg("exponent") );

    m.def("hodograph_scaled_exponential",
          &tsbm::hodographScaledExponential,
          py::arg("exponent"),
          py::arg("scale_factor"));

    m.def("hodograph_exponential_sine",
          &tsbm::hodographExponentialSine,
          py::arg("exponent"),
          py::arg("frequency") );

    m.def("hodograph_scaled_exponential_sine",
          &tsbm::hodographScaledExponentialSine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

    m.def("hodograph_exponential_cosine",
          &tsbm::hodographExponentialCosine,
          py::arg("exponent"),
          py::arg("frequency"));

    m.def("hodograph_scaled_exponential_cosine",
          &tsbm::hodographScaledExponentialCosine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

    m.def("hodograph_power",
          &tsbm::hodographPower,
          py::arg("exponent") );

    m.def("hodograph_scaled_power",
          &tsbm::hodographScaledPower,
          py::arg("exponent"),
          py::arg("scale_factor") );


    m.def("hodograph_power_sine",
          &tsbm::hodographPowerSine,
          py::arg("exponent"),
          py::arg("frequency") );

    m.def("hodograph_scaled_power_sine",
          &tsbm::hodographScaledPowerSine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

    m.def("hodograph_power_cosine",
          &tsbm::hodographPowerCosine,
          py::arg("exponent"),
          py::arg("frequency") );

    m.def("hodograph_scaled_power_cosine",
          &tsbm::hodographScaledPowerCosine,
          py::arg("exponent"),
          py::arg("frequency"),
          py::arg("scale_factor") );

    m.def("get_hodographic_shaping_acceleration_settings",
          &tsbm::getHodographicShapingAccelerationSettings,
          py::arg("initialState"),
          py::arg("finalState"),
          py::arg("timeOfFlight"),
          py::arg("centralBodyGravitationalParameter"),
          py::arg("numberOfRevolutions"),
          py::arg("radialVelocityFunctionComponents"),
          py::arg("normalVelocityFunctionComponents"),
          py::arg("axialVelocityFunctionComponents"),
          py::arg("freeCoefficientsRadialVelocityFunction"),
          py::arg("freeCoefficientsNormalVelocityFunction"),
          py::arg("freeCoefficientsAxialVelocityFunction"),
          py::arg("bodies"),
          py::arg("bodyToPropagate"),
          py::arg("specificImpulse"),
          py::arg("timeOffset") );

    m.def("get_hodographic_state_at_epoch",
          &tsbm::getHodographicShapingStateAtEpoch,
          py::arg("initialState"),
          py::arg("finalState"),
          py::arg("timeOfFlight"),
          py::arg("centralBodyGravitationalParameter"),
          py::arg("numberOfRevolutions"),
          py::arg("radialVelocityFunctionComponents"),
          py::arg("normalVelocityFunctionComponents"),
          py::arg("axialVelocityFunctionComponents"),
          py::arg("freeCoefficientsRadialVelocityFunction"),
          py::arg("freeCoefficientsNormalVelocityFunction"),
          py::arg("freeCoefficientsAxialVelocityFunction"),
          py::arg("epoch") );

    m.def("get_hodographic_delta_v",
          &tsbm::getHodographicShapingDeltaV,
          py::arg("initialState"),
          py::arg("finalState"),
          py::arg("timeOfFlight"),
          py::arg("centralBodyGravitationalParameter"),
          py::arg("numberOfRevolutions"),
          py::arg("radialVelocityFunctionComponents"),
          py::arg("normalVelocityFunctionComponents"),
          py::arg("axialVelocityFunctionComponents"),
          py::arg("freeCoefficientsRadialVelocityFunction"),
          py::arg("freeCoefficientsNormalVelocityFunction"),
          py::arg("freeCoefficientsAxialVelocityFunction") );

};

}// namespace tudatpy
