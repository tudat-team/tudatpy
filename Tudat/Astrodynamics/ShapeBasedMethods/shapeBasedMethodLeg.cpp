/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/ShapeBasedMethods/shapeBasedMethodLeg.h"

namespace tudat
{

namespace shape_based_methods
{


basic_astrodynamics::AccelerationMap ShapeBasedMethodLeg::retrieveLowThrustAccelerationMap(
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            getLowThrustAccelerationModel( specificImpulseFunction, integratorSettings );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationSettingsMap;
    accelerationSettingsMap[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = accelerationSettingsMap;

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

    accelerationModelMap[ bodyToPropagate_ ][ bodyToPropagate_ ].push_back( lowThrustAccelerationModel );

    return accelerationModelMap;

}


//! Returns state history.
void ShapeBasedMethodLeg::getTrajectory(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::Vector6d >& propagatedTrajectory )
{
    propagatedTrajectory.clear( );

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the mass profile of a shape-based trajectories, "
                                      "epochs are not provided in increasing order." );
        }

        double currentIndependentVariable = convertTimeToIndependentVariable( epochsVector[ i ] );
        propagatedTrajectory[ epochsVector[ i ] ] = computeCurrentStateVector( currentIndependentVariable );
    }
}


//! Compute current thrust vector.
Eigen::Vector3d ShapeBasedMethodLeg::computeCurrentThrust( double time,
                                                    std::function< double ( const double ) > specificImpulseFunction,
                                                    std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    double independentVariable = convertTimeToIndependentVariable( time );

    Eigen::Vector3d currentThrustVector = computeCurrentMass( 0.0, time, initialMass_, specificImpulseFunction, integratorSettings )
            * computeCurrentThrustAccelerationMagnitude( independentVariable, specificImpulseFunction, integratorSettings )
            * computeCurrentThrustAccelerationDirection( independentVariable, specificImpulseFunction, integratorSettings );

    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass_ );
    return currentThrustVector;
}


//! Return thrust profile.
void ShapeBasedMethodLeg::getThrustProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    thrustProfile.clear( );
    std::map< double, Eigen::VectorXd > massProfile;

    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectories, "
                                      "epochs are not provided in increasing order." );
        }

        double independentVariable = convertTimeToIndependentVariable( epochsVector[ i ] );

        double currentMass = massProfile[ epochsVector[ i ] ][ 0 ];
        thrustProfile[ epochsVector[ i ] ] = currentMass * computeCurrentThrustAccelerationMagnitude( independentVariable, specificImpulseFunction, integratorSettings )
                * computeCurrentThrustAccelerationDirection( independentVariable, specificImpulseFunction, integratorSettings );
    }
}


////! Compute current mass of the spacecraft.
//double ShapeBasedMethodLeg::computeCurrentMass(
//        const double independentVariable,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings )
//{
//    // Retrieve acceleration map.
//    basic_astrodynamics::AccelerationMap accelerationMap = retrieveLowThrustAccelerationMap( specificImpulseFunction );

//    // Create mass rate models
//    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
//    massRateModels[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
//                                                       bodyMap_, accelerationMap );

//    // Define mass propagator settings.
//    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettings =
//            std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
//                ( Eigen::Vector1d() << initialMass_ ).finished(),
//                std::make_shared< propagators::PropagationTimeTerminationSettings >( independentVariable, true ) );

//    // Re-initialise integrator settings.
//    integratorSettings->initialTime_ = 0.0;
//    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

//    // Create dynamics simulation object.
//    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
//                bodyMap_, integratorSettings, massPropagatorSettings, true, false, false );

//    // Propagate spacecraft mass.
//    std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//    double currentMass = propagatedMass.rbegin()->second[ 0 ];

//    return currentMass;

//}


////! Compute current mass of the spacecraft between two epochs.
//double ShapeBasedMethodLeg::computeCurrentMass(
//        const double timeInitialEpoch,
//        const double timeFinalEpoch,
//        const double massInitialEpoch,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    // Retrieve acceleration map.
//    basic_astrodynamics::AccelerationMap accelerationMap = retrieveLowThrustAccelerationMap( specificImpulseFunction );

//    // Create mass rate models
//    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
//    massRateModels[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
//                                                       bodyMap_, accelerationMap );

//    // Define mass propagator settings.
//    std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettings =
//            std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate_ }, massRateModels,
//                ( Eigen::Vector1d() << massInitialEpoch ).finished(),
//                std::make_shared< propagators::PropagationTimeTerminationSettings >( timeFinalEpoch, true ) );

//    // Re-initialise integrator settings.
//    integratorSettings->initialTime_ = timeInitialEpoch;
//    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

//    // Create dynamics simulation object.
//    propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
//                bodyMap_, integratorSettings, massPropagatorSettings, true, false, false );

//    // Propagate spacecraft mass.
//    std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//    double currentMass = propagatedMass.rbegin()->second[ 0 ];

//    return currentMass;
//}


////! Return mass profile.
//void ShapeBasedMethodLeg::getMassProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& massProfile,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    massProfile.clear( );

//    double currentMass = initialMass_;

//    for ( int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
//        {
//            throw std::runtime_error( "Error when retrieving the mass profile of a shape-based trajectories, "
//                                      "epochs are not provided in increasing order." );
//        }

//        if ( i == 0 )
//        {
//            currentMass = computeCurrentMass( 0.0, epochsVector[ i ], initialMass_, specificImpulseFunction, integratorSettings );
//            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished( );
//        }
//        else
//        {
//            currentMass = computeCurrentMass( epochsVector[ i - 1 ], epochsVector[ i ], currentMass, specificImpulseFunction, integratorSettings );
//            massProfile[ epochsVector[ i ] ] = ( Eigen::Vector1d( ) << currentMass ).finished();
//        }
//    }

//}


////! Compute current thrust vector.
//Eigen::Vector3d ShapeBasedMethodLeg::computeCurrentThrustAcceleration( double time )
//{
//    double independentVariable = convertTimeToIndependentVariable( time );
//    return computeCurrentThrustAccelerationMagnitude( independentVariable ) * computeCurrentThrustAccelerationDirection( independentVariable );
//}


////! Return thrust acceleration profile.
//void ShapeBasedMethodLeg::getThrustAccelerationProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile )
//{
//    thrustAccelerationProfile.clear();

//    for ( int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
//        {
//            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectories, "
//                                      "epochs are not provided in increasing order." );
//        }

//        Eigen::Vector3d currentThrustAccelerationVector = computeCurrentThrustAcceleration( epochsVector[ i ] );
//        thrustAccelerationProfile[ epochsVector[ i ] ] = currentThrustAccelerationVector;

//    }

//}

////! Compute current thrust vector.
//Eigen::Vector3d ShapeBasedMethodLeg::computeCurrentThrust( double time,
//                                                           std::function< double ( const double ) > specificImpulseFunction,
//                                                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//{
//    double independentVariable = convertTimeToIndependentVariable( time );
//    return computeCurrentMass( 0.0, time, initialMass_, specificImpulseFunction, integratorSettings )
//            * computeCurrentThrustAccelerationMagnitude( independentVariable ) * computeCurrentThrustAccelerationDirection( independentVariable );
//}


////! Return thrust profile.
//void ShapeBasedMethodLeg::getThrustProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& thrustProfile,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings)
//{
//    thrustProfile.clear( );

//    // Retrieve corresponding mass profile.
//    std::map< double, Eigen::VectorXd > massProfile;
//    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );
//    std::vector< double > massProfileVector;
//    for ( std::map< double, Eigen::VectorXd >::iterator itr = massProfile.begin( ) ; itr != massProfile.end( ) ; itr++ )
//    {
//        massProfileVector.push_back( itr->second[ 0 ] );
//    }

//    for ( int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
//        {
//            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectories, "
//                                      "epochs are not provided in increasing order." );
//        }

//        Eigen::Vector3d currentThrustVector = computeCurrentThrustAcceleration( epochsVector[ i ] ) * massProfileVector[ i ];
//        thrustProfile[ epochsVector[ i ] ] = currentThrustVector;

//    }
//}


//double ShapeBasedMethodLeg::computeDeltaV( )
//{
//    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of the azimuth angle.
//    std::function< double( const double ) > derivativeFunctionDeltaV = [ = ] ( const double currentIndependentVariable ){

//        double thrustAcceleration = computeCurrentThrustAccelerationMagnitude( currentIndependentVariable ) / physical_constants::ASTRONOMICAL_UNIT
//                * std::pow( physical_constants::JULIAN_YEAR, 2.0 ); // computeThrustAccelerationInSphericalCoordinates( currentAzimuthAngle ).norm()
////                * std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthAngle )
////                             * std::pow( radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ), 2.0 )
////                             / centralBodyGravitationalParameter_ );

//        return thrustAcceleration;

//    };

//    double finalValueIndependentValue = getFinalValueInpendentVariable( );

//    // Define numerical quadrature from quadratrure settings.
//    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
//            numerical_quadrature::createQuadrature( derivativeFunctionDeltaV, quadratureSettings_, finalValueIndependentValue );

//    // Return dimensional deltaV
//    return quadrature->getQuadrature( ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;
//}


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > ShapeBasedMethodLeg::createLowThrustTranslationalStatePropagatorSettings(
        basic_astrodynamics::AccelerationMap accelerationModelMap,
        std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );

    // Compute state vector at half of the time of flight.
    double independentVariableAtHalfTimeOfFlight = convertTimeToIndependentVariable( timeOfFlight_ / 2.0 );
    Eigen::Vector6d stateAtHalfOfTimeOfFlight = computeCurrentStateVector( independentVariableAtHalfTimeOfFlight );

    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
              stateAtHalfOfTimeOfFlight, terminationConditions.first, propagators::cowell, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
              stateAtHalfOfTimeOfFlight, terminationConditions.second, propagators::cowell, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}


} // namespace transfer_trajectories

} // namespace tudat
