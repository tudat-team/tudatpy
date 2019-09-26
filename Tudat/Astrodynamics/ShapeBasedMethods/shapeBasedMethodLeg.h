/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef SHAPEBASEDMETHOD_H
#define SHAPEBASEDMETHOD_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "Tudat/Astrodynamics/LowThrustDirectMethods/lowThrustLeg.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"

using namespace tudat::transfer_trajectories;

namespace tudat
{
namespace shape_based_methods
{


// Base class for shape based methods.
class ShapeBasedMethodLeg : public LowThrustLeg
{
public:

    //! Empty constructor.
    ShapeBasedMethodLeg( const Eigen::Vector6d& stateAtDeparture,
                      const Eigen::Vector6d& stateAtArrival,
                      const double timeOfFlight,
                      simulation_setup::NamedBodyMap& bodyMap,
                      const std::string bodyToPropagate,
                      const std::string centralBody,
                      std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings
                         = std::shared_ptr< numerical_integrators::IntegratorSettings< > >( ) ) :
    LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, bodyMap, bodyToPropagate, centralBody ),
    integratorSettings_( integratorSettings ){ }

    //! Default destructor.
    virtual ~ShapeBasedMethodLeg( ) { }

    //! Returns initial value of the independent variable.
    virtual double getInitialValueInpendentVariable( ) = 0;

    //! Returns final value of the independent variable.
    virtual double getFinalValueInpendentVariable( ) = 0;

//    //! Compute current thrust acceleration in cartesian coordinates.
//    virtual Eigen::Vector3d computeCurrentThrustAccelerationVector( const double independentVariable ) = 0;

//    //! Compute deltaV.
//    double computeDeltaV( );

//    //! Compute current cartesian state.
//    virtual Eigen::Vector6d computeCurrentStateVector( const double independentVariable ) = 0;

    //! Returns state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory );

//    //! Compute current mass of the spacecraft.
//    double computeCurrentMass( const double independentVariable,
//                               std::function< double ( const double ) > specificImpulseFunction,
//                               std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    //! Return mass profile.
//    void getMassProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& massProfile,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    //! Compute current thrust vector.
//    Eigen::Vector3d computeCurrentThrustAcceleration( double independentVariable );

//    //! Return thrust acceleration profile.
//    void getThrustAccelerationProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile );

//    //! Compute current thrust vector.
//    Eigen::Vector3d computeCurrentThrust(
//            double time,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    //! Return thrust profile.
//    void getThrustProfile( std::vector< double >& epochsVector,
//                           std::map< double, Eigen::VectorXd >& thrustProfile,
//                           std::function< double ( const double ) > specificImpulseFunction,
//                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap( std::function< double ( const double ) > specificImpulseFunction );

    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            basic_astrodynamics::AccelerationMap accelerationModelMap,
            std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


protected:

//    //! Compute current mass of the spacecraft between two epochs.
//    double computeCurrentMass(
//            const double timeInitialEpoch,
//            const double timeFinalEpoch,
//            const double massInitialEpoch,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Numerical quadrature settings, required to compute the time of flight and total deltaV.
    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_;

private:

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;




};


} // namespace shape_based_methods
} // namespace tudat

#endif // SHAPEBASEDMETHOD_H
