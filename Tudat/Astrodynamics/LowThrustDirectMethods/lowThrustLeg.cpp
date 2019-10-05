/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustDirectMethods/lowThrustLeg.h"

namespace tudat
{

namespace transfer_trajectories
{

//! Retrieve acceleration model (thrust).
std::shared_ptr< propulsion::ThrustAcceleration > LowThrustLeg::getLowThrustAccelerationModel(
        std::function< double( const double ) > specificImpulseFunction )
{

    std::shared_ptr< simulation_setup::Body > vehicle = bodyMap_[ bodyToPropagate_ ];

    // Define thrust magnitude function from the shaped trajectory.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {

        // Compute current independent variable.
        double currentIndependentVariable = convertTimeToIndependentVariable( currentTime );

        // Compute current acceleration.
        double currentAcceleration = computeCurrentThrustAccelerationMagnitude( currentIndependentVariable );

        // Compute current mass of the vehicle.
        double currentMass = vehicle->getBodyMass();

        // Compute and return magnitude of the low-thrust force.
        return currentAcceleration * currentMass;
    };

    // Define thrust magnitude settings from thrust magnitude function.
    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction );


    // Define thrust direction function from the shaped trajectory.
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction = [ = ]( const double currentTime )
    {
        // Compute current independent variable.
        double currentIndependentVariable = convertTimeToIndependentVariable( currentTime );

        // Compute current direction of the acceleration vector.
        Eigen::Vector3d currentAccelerationDirection = computeCurrentThrustAccelerationDirection( currentIndependentVariable );

        // Return direction of the low-thrust acceleration.
        return currentAccelerationDirection;
    };

    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            createThrustAcceleratioModel( thrustAccelerationSettings, bodyMap_, bodyToPropagate_ );

    return lowThrustAccelerationModel;

}


//! Compute current mass of the spacecraft between two epochs.
double LowThrustLeg::computeCurrentMass(
        const double timeInitialEpoch,
        const double timeFinalEpoch,
        const double massInitialEpoch,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    // Retrieve acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMap = retrieveLowThrustAccelerationMap( specificImpulseFunction ); // hybridMethodLeg_->getLowThrustTrajectoryAccelerationMap( );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationMap );

    // Define mass propagator settings.
    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
                ( Eigen::Vector1d() << massInitialEpoch ).finished(),
                std::make_shared< propagators::PropagationTimeTerminationSettings >( timeFinalEpoch, true ) );

    // Re-initialise integrator settings.
    integratorSettings->initialTime_ = timeInitialEpoch;
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

    // Create dynamics simulation object.
    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodyMap_, integratorSettings, massPropagatorSettings, true, false, false );

    // Propagate spacecraft mass.
    std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    double currentMass = propagatedMass.rbegin()->second[ 0 ];

    return currentMass;

}

//! Compute current mass of the spacecraft.
double LowThrustLeg::computeCurrentMass( const double independentVariable,
                           std::function< double ( const double ) > specificImpulseFunction,
                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    return computeCurrentMass( 0.0, independentVariable, initialMass_, specificImpulseFunction, integratorSettings );
}

//! Return mass profile.
void LowThrustLeg::getMassProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& massProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    massProfile.clear( );

    double currentMass = initialMass_;

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the mass profile of a hybrid trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        if ( i == 0 )
        {
            currentMass = computeCurrentMass( 0.0, epochsVector[ i ], currentMass, specificImpulseFunction, integratorSettings );
            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished( );
        }
        else
        {
            currentMass = computeCurrentMass( epochsVector[ i - 1 ], epochsVector[ i ], currentMass, specificImpulseFunction, integratorSettings );
            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished();
        }
    }

}

//! Compute current thrust vector.
Eigen::Vector3d LowThrustLeg::computeCurrentThrust( double time,
                                                    std::function< double ( const double ) > specificImpulseFunction,
                                                    std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    double independentVariable = convertTimeToIndependentVariable( time );
    return computeCurrentMass( 0.0, time, initialMass_, specificImpulseFunction, integratorSettings )
            * computeCurrentThrustAccelerationMagnitude( independentVariable ) * computeCurrentThrustAccelerationDirection( independentVariable );
}

//! Return thrust profile.
void LowThrustLeg::getThrustProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings)
{
    thrustProfile.clear( );

    // Retrieve corresponding mass profile.
    std::map< double, Eigen::VectorXd > massProfile;
    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );
    std::vector< double > massProfileVector;
    for ( std::map< double, Eigen::VectorXd >::iterator itr = massProfile.begin( ) ; itr != massProfile.end( ) ; itr++ )
    {
        massProfileVector.push_back( itr->second[ 0 ] );
    }

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectories, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector3d currentThrustVector = computeCurrentThrustAcceleration( epochsVector[ i ] ) * massProfileVector[ i ];
        thrustProfile[ epochsVector[ i ] ] = currentThrustVector;

    }
}


//! Return thrust acceleration profile.
void LowThrustLeg::getThrustAccelerationProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile )
{
    thrustAccelerationProfile.clear();

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a low-thrust trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector3d currentThrustAccelerationVector = computeCurrentThrustAcceleration( epochsVector[ i ] );
        thrustAccelerationProfile[ epochsVector[ i ] ] = currentThrustAccelerationVector;

    }

}


////! Full propagation.
//void LowThrustLeg::computeSemiAnalyticalAndFullPropagation(
//        std::function< double ( const double ) > specificImpulseFunction,
//        const std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings,
//        std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
//                std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
//        std::map< double, Eigen::VectorXd >& fullPropagationResults,
//        std::map< double, Eigen::VectorXd >& semiAnalyticalResults,
//        std::map< double, Eigen::VectorXd>& dependentVariablesHistory )





} // namespace transfer_trajectories

} // namespace tudat
