/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SHAPE_BASED_METHOD_H
#define TUDAT_SHAPE_BASED_METHOD_H

#include <tudat/simulation/simulation.h>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "tudat/astro/LowThrustTrajectories/lowThrustLeg.h"
#include "tudat/math/quadrature/createNumericalQuadrature.h"

using namespace tudat::low_thrust_trajectories;

namespace tudat
{
namespace shape_based_methods
{


// Base class for shape based methods.
class ShapeBasedMethod : public LowThrustLeg
{
public:

    //! Empty constructor.
    ShapeBasedMethod( const Eigen::Vector6d& stateAtDeparture,
                      const Eigen::Vector6d& stateAtArrival,
                      const double timeOfFlight,
                      const double initialBodyMass = TUDAT_NAN ) :
    LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, initialBodyMass, false ){ }

    //! Default destructor.
    virtual ~ShapeBasedMethod( ) { }

    //! Returns initial value of the independent variable.
    virtual double getInitialValueInpendentVariable( ) = 0;

    //! Returns final value of the independent variable.
    virtual double getFinalValueInpendentVariable( ) = 0;

    //! Returns state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory );

    Eigen::Vector3d computeCurrentThrustForce( double time,
                                          std::function< double ( const double ) > specificImpulseFunction,
                                          std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust profile.
    void getThrustForceProfile( std::vector< double >& epochsVector,
                           std::map< double, Eigen::VectorXd >& thrustProfile,
                           std::function< double ( const double ) > specificImpulseFunction,
                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
            const simulation_setup::NamedBodyMap& bodyMapTest,
            const std::string& bodyToPropagate,
            const std::string& centralBody,
            const std::function< double ( const double ) > specificImpulseFunction,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            const std::string& bodyToPropagate,
            const std::string& centralBody,
            const basic_astrodynamics::AccelerationMap& accelerationModelMap,
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


protected:

    //! Numerical quadrature settings, required to compute the time of flight and total deltaV.
    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_;


};


} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_SHAPE_BASED_METHOD_H
