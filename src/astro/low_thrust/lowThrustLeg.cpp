/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/low_thrust/lowThrustLeg.h"

namespace tudat
{

namespace low_thrust_trajectories
{


//Eigen::Vector3d LowThrustLeg::getThrustDirection( const double timeOffset )
//{
//    // Compute current independent variable.
//    double currentIndependentVariable = convertTimeToAzimuth( currentTime - timeOffset  );

//    // Compute current direction of the acceleration vector.
//    Eigen::Vector3d currentAccelerationDirection = computeCurrentThrustAccelerationDirectionFromAzimuth(
//                currentIndependentVariable );

//    // Return direction of the low-thrust acceleration.
//    return currentAccelerationDirection;
//}


////! Compute current mass of the spacecraft between two epochs.
//double LowThrustLeg::computeCurrentMass(
//        const double timeInitialEpoch,
//        const double timeFinalEpoch,
//        const double massInitialEpoch,
//        const std::function< double ( const double ) > specificImpulseFunction,
//        const std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    simulation_setup::SystemOfBodies bodies;
//    std::string bodyToPropagate = "Vehicle";
//    bodies.createEmptyBody( bodyToPropagate );
//    bodies.at( bodyToPropagate )->setConstantBodyMass( massInitialEpoch );


//    // Retrieve acceleration map.
//    basic_astrodynamics::AccelerationMap accelerationMap;
//    accelerationMap[ bodyToPropagate ][ bodyToPropagate ].push_back( getLowThrustAccelerationModel(
//                                                                         bodies, bodyToPropagate,  specificImpulseFunction, integratorSettings ) );

//    // Create mass rate models
//    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
//    massRateModels[ bodyToPropagate ] = createMassRateModel(
//                bodyToPropagate, std::make_shared< simulation_setup::FromThrustMassRateSettings >( 1 ),
//                bodies, accelerationMap );

//    // Define mass propagator settings.
//    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettings =
//            std::make_shared< propagators::MassPropagatorSettings< double > >(
//                std::vector< std::string >{ bodyToPropagate }, massRateModels,
//                ( Eigen::Vector1d() << massInitialEpoch ).finished(),
//                std::make_shared< propagators::PropagationTimeTerminationSettings >( timeFinalEpoch, true ) );

//    // Re-initialise integrator settings.
//    integratorSettings->initialTime_ = timeInitialEpoch;
//    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

//    // Create dynamics simulation object.
//    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
//                bodies, integratorSettings, massPropagatorSettings, true, false, false );

//    // Propagate spacecraft mass.
//    return dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second[ 0 ];

//}

////! Compute current mass of the spacecraft.
//double LowThrustLeg::computeCurrentMass( const double independentVariable,
//                                         std::function< double ( const double ) > specificImpulseFunction,
//                                         std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    return computeCurrentMass( 0.0, independentVariable, initialMass_, specificImpulseFunction, integratorSettings );
//}

////! Return mass profile.
//void LowThrustLeg::getMassProfile(
//        const std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& massProfile,
//        const std::function< double ( const double ) > specificImpulseFunction,
//        const std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    massProfile.clear( );

//    if( std::isnan( initialMass_ ) )
//    {
//        throw std::runtime_error( "Error when getting mass profile, initial mass is NaN" );
//    }
//    double currentMass = initialMass_;

//    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector.at( i ) < epochsVector.at( i - 1 ) ) )
//        {
//            throw std::runtime_error( "Error when retrieving the mass profile of a low-thrust trajectory, "
//                                      "epochs are not provided in increasing order." );
//        }

//        if ( i == 0 )
//        {
//            currentMass = computeCurrentMass( 0.0, epochsVector.at( i ), currentMass, specificImpulseFunction, integratorSettings );
//            massProfile[ epochsVector.at( i ) ] = ( Eigen::Vector1d( ) << currentMass ).finished( );
//        }
//        else
//        {
//            currentMass = computeCurrentMass( epochsVector.at( i - 1 ), epochsVector.at( i ), currentMass, specificImpulseFunction, integratorSettings );
//            massProfile[ epochsVector.at( i ) ] = ( Eigen::Vector1d( ) << currentMass ).finished();
//        }
//    }

//}

////! Return thrust profile.
//void LowThrustLeg::getThrustForceProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& thrustProfile,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings)
//{
//    thrustProfile.clear( );

//    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector.at( i ) < epochsVector.at( i - 1 ) ) )
//        {
//            throw std::runtime_error( "Error when retrieving the thrust profile of a low-thrust trajectory, "
//                                      "epochs are not provided in increasing order." );
//        }

//        Eigen::Vector3d currentThrustVector = computeCurrentThrustForce( epochsVector.at( i ), specificImpulseFunction, integratorSettings );
//        thrustProfile[ epochsVector.at( i ) ] = currentThrustVector;

//    }
//}


////! Return thrust acceleration profile.
//void LowThrustLeg::getThrustAccelerationProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    thrustAccelerationProfile.clear();

//    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector.at( i ) < epochsVector.at( i - 1 ) ) )
//        {
//            throw std::runtime_error( "Error when retrieving the thrust profile of a low-thrust trajectory, "
//                                      "epochs are not provided in increasing order." );
//        }

//        Eigen::Vector3d currentThrustAccelerationVector =
//                computeCurrentThrustAccelerationFromAzimuth( epochsVector.at( i ), specificImpulseFunction, integratorSettings );
//        thrustAccelerationProfile[ epochsVector.at( i ) ] = currentThrustAccelerationVector;

//    }

//}




} // namespace low_thrust_trajectories

} // namespace tudat
