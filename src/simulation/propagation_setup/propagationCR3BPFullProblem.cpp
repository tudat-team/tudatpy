/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "tudat/astro/gravitation/librationPoint.h"
#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "tudat/simulation/propagation_setup/propagationCR3BPFullProblem.h"

namespace tudat
{

namespace propagators
{


simulation_setup::BodyListSettings setupBodySettingsCR3BP(
        const double distancePrimarySecondary,
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& frameOrientation,
        const double primaryGravitationalParameter,
        const double secondaryGravitationalParameter )
{
    spice_interface::loadStandardSpiceKernels( );

    double gravitationalParameterPrimary = ( primaryGravitationalParameter == primaryGravitationalParameter ) ?
                primaryGravitationalParameter :
                simulation_setup::createGravityFieldModel(
                    simulation_setup::getDefaultGravityFieldSettings( namePrimaryBody, TUDAT_NAN, TUDAT_NAN ),
                    namePrimaryBody )->getGravitationalParameter( );
    double gravitationalParameterSecondary =( secondaryGravitationalParameter == secondaryGravitationalParameter ) ?
                secondaryGravitationalParameter :
                simulation_setup::createGravityFieldModel(
                    simulation_setup::getDefaultGravityFieldSettings( namePrimaryBody, TUDAT_NAN, TUDAT_NAN ),
                    nameSecondaryBody )->getGravitationalParameter( );
    double massParameter = circular_restricted_three_body_problem::computeMassParameter(
                gravitationalParameterPrimary, gravitationalParameterSecondary );

    // Initial state for the primary
    Eigen::Vector6d initialStateInKeplerianElementsPrimary = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElementsPrimary( orbital_element_conversions::semiMajorAxisIndex ) =
            massParameter * distancePrimarySecondary;
    initialStateInKeplerianElementsPrimary( orbital_element_conversions::trueAnomalyIndex ) = mathematical_constants::PI;

    // Initial state for the secondary
    Eigen::Vector6d initialStateInKeplerianElementsSecondary = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElementsSecondary( orbital_element_conversions::semiMajorAxisIndex ) =
            ( 1.0 - massParameter) * distancePrimarySecondary;

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( namePrimaryBody );
    bodiesToCreate.push_back( nameSecondaryBody );
    simulation_setup::BodyListSettings bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    // Compute effective gravitational parameters
    double sumGravitationalParameter = gravitationalParameterPrimary + gravitationalParameterSecondary;
    double distanceBarycenterPrimary = gravitationalParameterSecondary * distancePrimarySecondary / ( sumGravitationalParameter );
    double distanceBarycenterSecondary = distancePrimarySecondary - distanceBarycenterPrimary;
    double gravitationalParameterPrimaryTwoBodyProblem =
            std::pow( distanceBarycenterPrimary, 3 ) * sumGravitationalParameter /
            std::pow( distanceBarycenterPrimary + distanceBarycenterSecondary, 3 );
    double gravitationalParameterSecondaryTwoBodyProblem =
            std::pow( distanceBarycenterSecondary, 3 ) * sumGravitationalParameter /
            std::pow( distanceBarycenterPrimary + distanceBarycenterSecondary, 3 );


    // Define body ephemeris settings
    bodySettings.at( namePrimaryBody )->ephemerisSettings = std::make_shared< simulation_setup::KeplerEphemerisSettings >(
                initialStateInKeplerianElementsPrimary, 0.0, gravitationalParameterPrimaryTwoBodyProblem,
                "SSB", frameOrientation );
    bodySettings.at( nameSecondaryBody )->ephemerisSettings = std::make_shared< simulation_setup::KeplerEphemerisSettings >(
                initialStateInKeplerianElementsSecondary, 0.0, gravitationalParameterSecondaryTwoBodyProblem,
                "SSB", frameOrientation );
    for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
    {
        bodySettings.at( bodiesToCreate.at( j ) )->ephemerisSettings->resetFrameOrientation( frameOrientation );
        bodySettings.at( bodiesToCreate.at( j ) )->rotationModelSettings->resetOriginalFrame( frameOrientation );
    }

    return bodySettings;
}

//! Function to directly setup CR3BP bodies
simulation_setup::SystemOfBodies setupBodyMapCR3BP(
        const double distancePrimarySecondary,
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const std::string& frameOrientation,
        const double primaryGravitationalParameter,
        const double secondaryGravitationalParameter )
{
    // Create system of bodies
    simulation_setup::SystemOfBodies bodies = createSystemOfBodies(
                setupBodySettingsCR3BP( distancePrimarySecondary, namePrimaryBody, nameSecondaryBody, frameOrientation,
                                        primaryGravitationalParameter, secondaryGravitationalParameter ) );
    bodies.createEmptyBody( nameBodyToPropagate );
    bodies.at( nameBodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                      std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                      < double, Eigen::Vector6d > >( ), "SSB", frameOrientation ) );

    return bodies;

}


//! Setup CR3BP acceleration map.
basic_astrodynamics::AccelerationMap setupAccelerationMapCR3BP(
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const std::string& centralBody,
        const simulation_setup::SystemOfBodies& bodies )
{

    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ namePrimaryBody ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                   basic_astrodynamics::point_mass_gravity ) );
    bodyToPropagateAccelerations[ nameSecondaryBody ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                     basic_astrodynamics::point_mass_gravity ) );
    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ nameBodyToPropagate ] = bodyToPropagateAccelerations;

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, { nameBodyToPropagate }, { centralBody } );

    return accelerationModelMap;

}



//! Propagate CR3BP from CR3BP environment
void propagateCR3BPFromEnvironment(
        const double initialTime,
        const double finalPropagationTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP,
        std::map< double, Eigen::Vector6d >& stateHistory,
        const bool outputInNormalizedCoordinates )
{
    // CR3BP definition
    double gravitationalParameterPrimary;
    double gravitationalParameterSecondary;
    Eigen::Vector6d initialStatePrimary;
    Eigen::Vector6d initialStateSecondary;

    if( bodies.at( bodiesCR3BP.at( 0 ) )->getBodyMass( ) > bodies.at( bodiesCR3BP.at( 1 ) )->getBodyMass( ) )
    {
        gravitationalParameterPrimary = bodies.at( bodiesCR3BP.at( 0 ) )->getGravityFieldModel( )->getGravitationalParameter( );
        initialStatePrimary = bodies.at( bodiesCR3BP.at( 0 ) )->getEphemeris( )->getCartesianState( initialTime );
        gravitationalParameterSecondary = bodies.at( bodiesCR3BP.at( 1 ) )->getGravityFieldModel( )->getGravitationalParameter( );
        initialStateSecondary = bodies.at( bodiesCR3BP.at( 1 ) )->getEphemeris( )->getCartesianState( initialTime );
    }
    else
    {
        gravitationalParameterPrimary = bodies.at( bodiesCR3BP.at( 1 ) )->getGravityFieldModel( )->getGravitationalParameter( );
        initialStatePrimary = bodies.at( bodiesCR3BP.at( 1 ) )->getEphemeris( )->getCartesianState( initialTime );
        gravitationalParameterSecondary = bodies.at( bodiesCR3BP.at( 0 ) )->getGravityFieldModel( )->getGravitationalParameter( );
        initialStateSecondary = bodies.at( bodiesCR3BP.at( 0 ) )->getEphemeris( )->getCartesianState( initialTime );
    }


    double massParameter = circular_restricted_three_body_problem::computeMassParameter(
                gravitationalParameterPrimary, gravitationalParameterSecondary );
    double distanceBetweenPrimaries = ( initialStateSecondary.segment(0,3) - initialStatePrimary.segment(0,3) ).norm( );

    double normalizedInitialTime = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                initialTime, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);
    double normalizedFinalPropagationTime = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                finalPropagationTime, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);
    Eigen::Vector6d normalizedInitialState =
            circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
                gravitationalParameterSecondary, gravitationalParameterPrimary,
                distanceBetweenPrimaries, initialState, initialTime );



    // CR3BP propagation
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > CR3BPintegratorSettings = integratorSettings;
    double originalInitialTime = CR3BPintegratorSettings->initialTime_;
    double originalInitialTimeStep = CR3BPintegratorSettings->initialTimeStep_;

    CR3BPintegratorSettings->initialTime_ = normalizedInitialTime;
    CR3BPintegratorSettings->initialTimeStep_ = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                integratorSettings->initialTimeStep_, gravitationalParameterPrimary, gravitationalParameterSecondary,
                distanceBetweenPrimaries );

    stateHistory.clear( );
    if( outputInNormalizedCoordinates )
    {
        stateHistory = performCR3BPIntegration(
                    CR3BPintegratorSettings, massParameter, normalizedInitialState, normalizedFinalPropagationTime, true );
    }
    else
    {
        std::map< double, Eigen::Vector6d > normalizedStateHistory = performCR3BPIntegration(
                    CR3BPintegratorSettings, massParameter, normalizedInitialState, normalizedFinalPropagationTime, true );
        for( auto it = normalizedStateHistory.begin( ); it != normalizedStateHistory.end( ); it++ )
        {
            stateHistory[ circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(
                        it->first, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries ) ] =
                    circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates(
                        gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries,
                        it->second, it->first );
        }
    }

    CR3BPintegratorSettings->initialTime_ = originalInitialTime;
    CR3BPintegratorSettings->initialTimeStep_ = originalInitialTimeStep;
}

void propagateCR3BPAndFullDynamicsProblem(
        const double initialTime,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP,
        std::map< double, Eigen::Vector6d >& directPropagationResult,
        std::map< double, Eigen::Vector6d >& cr3bpPropagationResult,
        std::map< double, Eigen::VectorXd >& dependentVariableValues )
{
    // Propagate the full problem
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings );


    std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    utilities::castDynamicToFixedSizeEigenVectorMap< double, double, 6 >(
                stateHistory, directPropagationResult );
    dependentVariableValues = dynamicsSimulator.getDependentVariableHistory( );

    double finalPropagationTime = directPropagationResult.rbegin( )->first;

    cr3bpPropagationResult.clear( );
    propagateCR3BPFromEnvironment(
                initialTime, finalPropagationTime, propagatorSettings->getInitialStates( ), integratorSettings, bodies,
                bodiesCR3BP, cr3bpPropagationResult, false );
}

//! Propagate the CR3BP and the full dynamics problem
void propagateCR3BPAndFullDynamicsProblem(
        const double initialTime,
        const double finalTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP,
        std::map< double, Eigen::Vector6d >& directPropagationResult,
        std::map< double, Eigen::Vector6d >& cr3bpPropagationResult )
{
    std::shared_ptr< TranslationalStatePropagatorSettings< double> > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, initialState,
              std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );

    std::map< double, Eigen::VectorXd > dependentVariableValues;
    propagateCR3BPAndFullDynamicsProblem(
            initialTime, integratorSettings, propagatorSettings, bodies, bodiesCR3BP, directPropagationResult,
            cr3bpPropagationResult, dependentVariableValues );

}

//! Propagate the CR3BP and the full dynamics problem and compute the state difference at the end of the propagation.
Eigen::Vector6d getFinalStateDifferenceFullPropagationWrtCR3BP(
        const double initialTime,
        const double finalTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP )
{
    std::map< double, Eigen::Vector6d > directPropagationResult;
    std::map< double, Eigen::Vector6d > cr3bpPropagationResult;
    propagateCR3BPAndFullDynamicsProblem(
                initialTime, finalTime, initialState, integratorSettings, accelerationModelMap, bodiesToPropagate, centralBodies,
                bodies, bodiesCR3BP, directPropagationResult, cr3bpPropagationResult );

    Eigen::Vector6d finalPropagatedStateFullProblem = directPropagationResult.rbegin( )->second;
    Eigen::Vector6d finalPropagatedStateCR3BP = cr3bpPropagationResult.rbegin( )->second;

    // Difference between the two propagated states at finalTime
    return finalPropagatedStateCR3BP - finalPropagatedStateFullProblem;
}

}

}


