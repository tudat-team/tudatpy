/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/createStateDerivativeModel.h"
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"




//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "mainTestPropagationCR3BPfullProblem.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

namespace tudat
{

namespace propagators
{


//!Function to transform normalized co-rotating coordinates into cartesian ones
Eigen::Vector6d convertCorotatingNormalizedToCartesianCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& normalizedState,
        const double normalizedTime )
{
    Eigen::Vector3d normalizedPosition = normalizedState.segment( 0, 3 );
    Eigen::Vector3d normalizedVelocity = normalizedState.segment( 3, 3 );

    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.setZero( );
    Eigen::Matrix3d derivativeRotationMatrix;
    derivativeRotationMatrix.setZero( );

    rotationMatrix( 0, 0 ) = std::cos( normalizedTime );
    rotationMatrix( 0, 1 ) = - std::sin( normalizedTime );
    rotationMatrix( 1, 0 ) = std::sin( normalizedTime );
    rotationMatrix( 1, 1 ) = std::cos( normalizedTime );

    derivativeRotationMatrix( 0, 0 ) = - std::sin( normalizedTime );
    derivativeRotationMatrix( 0, 1 ) = - std::cos( normalizedTime );
    derivativeRotationMatrix( 1, 0 ) = std::cos( normalizedTime );
    derivativeRotationMatrix( 1, 1 ) = - std::sin( normalizedTime );

    Eigen::Vector6d inertialNormalizedState;
    inertialNormalizedState.segment( 0, 3 ) = rotationMatrix * normalizedPosition;
    inertialNormalizedState.segment( 3, 3 ) = derivativeRotationMatrix * normalizedPosition + rotationMatrix * normalizedVelocity;
    Eigen::Vector6d cartesianState = circular_restricted_three_body_problem::convertDimensionlessCartesianStateToDimensionalUnits(inertialNormalizedState, gravitationalParameterPrimary, gravitationalParameterSecondary, distancePrimarySecondary);

    return cartesianState;
}



//! Function to transform cartesian coordinates into co-rotating normalized ones
Eigen::Vector6d convertCartesianToCorotatingNormalizedCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& cartesianState,
        const double time )
{
    Eigen::Vector3d cartesianPosition = cartesianState.segment( 0, 3 );
    Eigen::Vector3d cartesianVelocity = cartesianState.segment( 3, 3 );

    double meanMotion = std::sqrt( ( gravitationalParameterPrimary + gravitationalParameterSecondary ) /
                                   std::pow( distancePrimarySecondary, 3 ) );

    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.setZero( );
    rotationMatrix( 0, 0 ) = std::cos( meanMotion * time );
    rotationMatrix( 0, 1 ) = std::sin( meanMotion * time );
    rotationMatrix( 1, 0 ) = -std::sin( meanMotion * time );
    rotationMatrix( 1, 1 ) = std::cos( meanMotion * time );

    Eigen::Matrix3d derivativeRotationMatrix;
    derivativeRotationMatrix.setZero( );
    derivativeRotationMatrix( 0, 0 ) = -std::sin( meanMotion * time );
    derivativeRotationMatrix( 0, 1 ) = std::cos( meanMotion * time );
    derivativeRotationMatrix( 1, 0 ) = -std::cos( meanMotion * time );
    derivativeRotationMatrix( 1, 1 ) = -std::sin( meanMotion * time );
    derivativeRotationMatrix = meanMotion * derivativeRotationMatrix;

    Eigen::Vector6d corotatingDimensionalState;
    corotatingDimensionalState.segment( 0, 3 ) = rotationMatrix * cartesianPosition;
    corotatingDimensionalState.segment( 3, 3 ) = derivativeRotationMatrix * cartesianPosition + rotationMatrix * cartesianVelocity;

    Eigen::Vector6d normalizedState = circular_restricted_three_body_problem::convertDimensionalCartesianStateToDimensionlessState(
                corotatingDimensionalState, gravitationalParameterPrimary, gravitationalParameterSecondary, distancePrimarySecondary );

    return normalizedState;
}



//! Function to directly setup CR3BP bodyMap
simulation_setup::NamedBodyMap setupBodyMapCR3BPBodyMap(
        const double distancePrimarySecondary,
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const double initialTime )
{

    spice_interface::loadStandardSpiceKernels( );


    double gravitationalParameterPrimary;
    double gravitationalParameterSecondary;

    // retrieve the gravitational parameter of the primary
    if ( namePrimaryBody == "Earth" || namePrimaryBody == "Mars" || namePrimaryBody == "Moon")
    {
        std::shared_ptr< simulation_setup::GravityFieldSettings > gravityFieldSettings = simulation_setup::getDefaultGravityFieldSettings( namePrimaryBody, 0.0, 1.0);
        std::shared_ptr< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings > modelGravityFieldSettings =
                std::dynamic_pointer_cast< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        gravitationalParameterPrimary = modelGravityFieldSettings->getGravitationalParameter( );
    }
    else
    {
        gravitationalParameterPrimary = spice_interface::getBodyGravitationalParameter( namePrimaryBody );
    }


    // retrieve the gravitational parameter of the secondary
    if ( nameSecondaryBody == "Earth" || nameSecondaryBody == "Mars" || nameSecondaryBody == "Moon")
    {
        std::shared_ptr< simulation_setup::GravityFieldSettings > gravityFieldSettings = simulation_setup::getDefaultGravityFieldSettings( nameSecondaryBody, 0.0, 1.0);
        std::shared_ptr< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings > modelGravityFieldSettings =
                std::dynamic_pointer_cast< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        gravitationalParameterSecondary = modelGravityFieldSettings->getGravitationalParameter( );
    }
    else
    {
        gravitationalParameterSecondary = spice_interface::getBodyGravitationalParameter( nameSecondaryBody );
    }



    double massParameter = circular_restricted_three_body_problem::computeMassParameter(
                gravitationalParameterPrimary, gravitationalParameterSecondary );

    std::string frameOrientation = "J2000";

    // Initial state for the primary
    Eigen::Vector6d initialStateInKeplerianElementsPrimary = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::semiMajorAxisIndex) = massParameter * distancePrimarySecondary;
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::trueAnomalyIndex) = mathematical_constants::PI;


    // Initial state for the secondary
    Eigen::Vector6d initialStateInKeplerianElementsSecondary = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElementsSecondary(orbital_element_conversions::semiMajorAxisIndex) =
            ( 1.0 - massParameter) * distancePrimarySecondary;


    // Create body objects.
    std::vector< std::string > bodiesToCreate;

    bodiesToCreate.push_back( namePrimaryBody );
    bodiesToCreate.push_back( nameSecondaryBody );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    double sumGravitationalParameter = gravitationalParameterPrimary + gravitationalParameterSecondary;
    double distanceBarycenterPrimary = gravitationalParameterSecondary * distancePrimarySecondary / ( sumGravitationalParameter );
    double distanceBarycenterSecondary = distancePrimarySecondary - distanceBarycenterPrimary;
    double gravitationalParameterPrimaryTwoBodyProblem =
            std::pow( distanceBarycenterPrimary, 3 ) * sumGravitationalParameter /
            std::pow( distanceBarycenterPrimary + distanceBarycenterSecondary, 3 );
    double gravitationalParameterSecondaryTwoBodyProblem =
            std::pow( distanceBarycenterSecondary, 3 ) * sumGravitationalParameter /
            std::pow( distanceBarycenterPrimary + distanceBarycenterSecondary, 3 );

    // Primary ephemeris
    bodySettings[ namePrimaryBody ]->ephemerisSettings = std::make_shared< simulation_setup::KeplerEphemerisSettings >(
                initialStateInKeplerianElementsPrimary, initialTime, gravitationalParameterPrimaryTwoBodyProblem, "SSB", frameOrientation );

    // Secondary ephemeris
    bodySettings[ nameSecondaryBody ]->ephemerisSettings = std::make_shared< simulation_setup::KeplerEphemerisSettings >(
                initialStateInKeplerianElementsSecondary, initialTime, gravitationalParameterSecondaryTwoBodyProblem, "SSB", frameOrientation );



    for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
    {
        bodySettings[ bodiesToCreate.at( j ) ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
        bodySettings[ bodiesToCreate.at( j ) ]->rotationModelSettings->resetOriginalFrame( frameOrientation );
    }


    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );


    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ nameBodyToPropagate ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                      std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                      < double, Eigen::Vector6d > >( ), "SSB", frameOrientation ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", frameOrientation );

    return bodyMap;

}


//! Function to directly setup CR3BP acceleration map
basic_astrodynamics::AccelerationMap setupAccelerationMapCR3BP(
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::NamedBodyMap& bodyMap )
{

    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ namePrimaryBody ].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                  basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ nameSecondaryBody ].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );
    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ nameBodyToPropagate ] = bodyToPropagateAccelerations;

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    return accelerationModelMap;

}



//! Function to simultaneously propagate the dynamics in the CR3BP and in the full dynamics problem
//! and compute the difference in state at the end of the propagation
Eigen::Vector6d propagateCR3BPandFullDynamicsProblem(
        const double initialTime,
        const double finalTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector < std::string >& bodiesCR3BP )
{
    // Create barycenter object to retrieve the position of the primary and secondary.
    bodyMap[ "TwoBodyBarycenter" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "TwoBodyBarycenter" ]->setEphemeris(std::make_shared< ephemerides::ConstantEphemeris >(
                                                     ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), "SSB", "J2000" ) );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_position_dependent_variable, bodiesCR3BP.at( 0 ), "TwoBodyBarycenter" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_velocity_dependent_variable, bodiesCR3BP.at( 0 ), "TwoBodyBarycenter" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_position_dependent_variable, bodiesCR3BP.at( 1 ), "TwoBodyBarycenter" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          relative_velocity_dependent_variable, bodiesCR3BP.at( 1 ), "TwoBodyBarycenter" ) );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );
    std::shared_ptr< TranslationalStatePropagatorSettings< double> > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, initialState,
              std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ), cowell,
              dependentVariablesToSave );

    // Propagate the full problem
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblem = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    double finalPropagationTime = stateHistoryFullProblem.rbegin( )->first;

    Eigen::Vector6d finalPropagatedStateFullProblem = stateHistoryFullProblem.rbegin( )->second.transpose( );
    Eigen::VectorXd initialStateBodies = dynamicsSimulator.getDependentVariableHistory( ).begin( )->second;


    // CR3BP definition
    double gravitationalParameterPrimary;
    double gravitationalParameterSecondary;
    Eigen::Vector6d initialStatePrimary;
    Eigen::Vector6d initialStateSecondary;

    if( bodyMap[ bodiesCR3BP.at( 0 ) ]->getBodyMass( ) > bodyMap[ bodiesCR3BP.at( 1 ) ]->getBodyMass( ) )
    {
        gravitationalParameterPrimary = bodyMap[ bodiesCR3BP.at( 0 ) ]->getGravityFieldModel( )->getGravitationalParameter( );
        initialStatePrimary = initialStateBodies.segment( 0, 6 );
        gravitationalParameterSecondary = bodyMap[ bodiesCR3BP.at( 1 ) ]->getGravityFieldModel( )->getGravitationalParameter( );
        initialStateSecondary = initialStateBodies.segment( 6, 6 );
    }
    else
    {
        gravitationalParameterPrimary = bodyMap[ bodiesCR3BP.at( 1 ) ]->getGravityFieldModel( )->getGravitationalParameter( );
        initialStatePrimary = initialStateBodies.segment( 6, 6 );
        gravitationalParameterSecondary = bodyMap[ bodiesCR3BP.at( 0 ) ]->getGravityFieldModel( )->getGravitationalParameter( );
        initialStateSecondary = initialStateBodies.segment( 0, 6 );
    }


    double massParameter = circular_restricted_three_body_problem::computeMassParameter(
                gravitationalParameterPrimary, gravitationalParameterSecondary );
    double distanceBetweenPrimaries = ( initialStateSecondary - initialStatePrimary ).norm( );

    double normalizedInitialTime = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                initialTime, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);
    double normalizedFinalPropagationTime = circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(
                finalPropagationTime, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);
    Eigen::Vector6d normalizedInitialState = convertCartesianToCorotatingNormalizedCoordinates(
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

    std::map< double, Eigen::Vector6d > stateHistoryCR3BP = performCR3BPIntegration(
                CR3BPintegratorSettings, massParameter, normalizedInitialState, normalizedFinalPropagationTime, true );
    Eigen::Vector6d normalizedFinalPropagatedStateCR3BP = stateHistoryCR3BP.rbegin( )->second;
    double normalizedFinalPropagationTimeCR3BP = stateHistoryCR3BP.rbegin( )->first;


    CR3BPintegratorSettings->initialTime_ = originalInitialTime;
    CR3BPintegratorSettings->initialTimeStep_ = originalInitialTimeStep;

    // Transformation to inertial coordinates
    Eigen::Vector6d finalPropagatedStateCR3BP = convertCorotatingNormalizedToCartesianCoordinates(
                gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries,
                normalizedFinalPropagatedStateCR3BP, normalizedFinalPropagationTimeCR3BP);

    // Difference between the two propagated states at finalTime
    Eigen::Vector6d finalStateDifference = finalPropagatedStateCR3BP - finalPropagatedStateFullProblem;


    return finalStateDifference;

}

}

}



int main( ){

    using namespace tudat;
    using namespace propagators;


    std::cout.precision( 20);

    double initialTime = 0.0;
    double finalTime = 120000000.0;

    std::vector < std::string > bodiesCR3BP;
    bodiesCR3BP.push_back("Sun");
    bodiesCR3BP.push_back("Earth");

    simulation_setup::NamedBodyMap bodyMap = setupBodyMapCR3BPBodyMap(
                physical_constants::ASTRONOMICAL_UNIT, "Sun", "Earth", "spacecraft", 0.0);

    // Spacecraft properties
    bodyMap[ "spacecraft" ]->setConstantBodyMass( 100.0);

    // Initialization of the spacecraft state
    Eigen::Vector6d initialState;
    initialState[0] = 2.991957413820000e+10;
    initialState[1] = 1.295555563704656e+11;
    initialState[2] = 0.0;
    initialState[3] = -2.579433850734350e+04;
    initialState[4] = 5.956947312313238e+03;
    initialState[5] = 0.0;


    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "spacecraft" );
    centralBodies.push_back( "SSB" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = setupAccelerationMapCR3BP(
                "Sun", "Earth", "spacecraft", bodiesToPropagate, centralBodies, bodyMap);



    // Create integrator settings
    const double fixedStepSize = 1000;

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);


    // calculate the difference between CR3BP and full problem
    Eigen::Vector6d stateDifference = propagateCR3BPandFullDynamicsProblem(initialTime, finalTime, initialState,
                                                                           integratorSettings, accelerationModelMap,
                                                                           bodiesToPropagate, centralBodies,
                                                                           bodyMap, bodiesCR3BP);

    std::cout << "stateDifference: " << stateDifference << "\n\n";

}


