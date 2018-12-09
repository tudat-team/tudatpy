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



using namespace tudat;


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


//!Function to transform normalized co-rotating coordinates into cartesian ones
Eigen::Vector6d convertCorotatingNormalizedToCartesianCoordinates(double gravitationalParameterPrimary, double gravitationalParameterSecondary, double distancePrimarySecondary, Eigen::Vector6d normalizedState, double normalizedTime)
{
    Eigen::Vector3d normalizedPosition = normalizedState.segment(0,3);
    Eigen::Vector3d normalizedVelocity = normalizedState.segment(3,3);

    Eigen::Matrix3d rotationMatrix;
    Eigen::Matrix3d derivativeRotationMatrix;

   rotationMatrix(0,0) = cos(normalizedTime);
   rotationMatrix(0,1) = - sin(normalizedTime);
   rotationMatrix(0,2) = 0.0;
   rotationMatrix(1,0) = sin(normalizedTime);
   rotationMatrix(1,1) = cos(normalizedTime);
   rotationMatrix(1,2) = 0.0;
   rotationMatrix(2,0) = 0.0;
   rotationMatrix(2,1) = 0.0;
   rotationMatrix(2,2) = 1.0;

   derivativeRotationMatrix(0,0) = - sin(normalizedTime);
   derivativeRotationMatrix(0,1) = - cos(normalizedTime);
   derivativeRotationMatrix(0,2) = 0;
   derivativeRotationMatrix(1,0) = cos(normalizedTime);
   derivativeRotationMatrix(1,1) = - sin(normalizedTime);
   derivativeRotationMatrix(1,2) = 0;
   derivativeRotationMatrix(2,0) = 0;
   derivativeRotationMatrix(2,1) = 0;
   derivativeRotationMatrix(2,2) = 0;

   Eigen::Vector6d inertialNormalizedState;
   inertialNormalizedState.segment(0,3) = rotationMatrix * normalizedPosition;
   inertialNormalizedState.segment(3,3) = derivativeRotationMatrix * normalizedPosition + rotationMatrix * normalizedVelocity;
   Eigen::Vector6d cartesianState = gravitation::circular_restricted_three_body_problem::convertDimensionlessCartesianStateToDimensionalUnits(inertialNormalizedState, gravitationalParameterPrimary, gravitationalParameterSecondary, distancePrimarySecondary);

   return cartesianState;
}



//! Function to transform cartesian coordinates into co-rotating normalized ones
Eigen::Vector6d convertCartesianToCorotatingNormalizedCoordinates(double gravitationalParameterPrimary, double gravitationalParameterSecondary, double distancePrimarySecondary, Eigen::Vector6d cartesianState, double time)
{
    Eigen::Vector3d cartesianPosition = cartesianState.segment(0,3);
    Eigen::Vector3d cartesianVelocity = cartesianState.segment(3,3);

    double n = std::sqrt(( gravitationalParameterPrimary + gravitationalParameterSecondary ) / std::pow(distancePrimarySecondary,3));

    Eigen::Matrix3d rotationMatrix;
    rotationMatrix(0,0) = cos( n * time );
    rotationMatrix(0,1) = sin( n * time );
    rotationMatrix(0,2) = 0.0;
    rotationMatrix(1,0) = - sin( n * time );
    rotationMatrix(1,1) = cos( n * time );
    rotationMatrix(1,2) = 0.0;
    rotationMatrix(2,0) = 0.0;
    rotationMatrix(2,1) = 0.0;
    rotationMatrix(2,2) = 1.0;

    Eigen::Matrix3d derivativeRotationMatrix;
    derivativeRotationMatrix(0,0) = - sin( n * time );
    derivativeRotationMatrix(0,1) = cos( n * time );
    derivativeRotationMatrix(0,2) = 0;
    derivativeRotationMatrix(1,0) = - cos( n * time );
    derivativeRotationMatrix(1,1) = - sin( n * time );
    derivativeRotationMatrix(1,2) = 0.0;
    derivativeRotationMatrix(2,0) = 0.0;
    derivativeRotationMatrix(2,1) = 0.0;
    derivativeRotationMatrix(2,2) = 0.0;
    derivativeRotationMatrix = n * derivativeRotationMatrix;

   Eigen::Vector6d corotatingDimensionalState;
   corotatingDimensionalState.segment(0,3) = rotationMatrix * cartesianPosition;
   corotatingDimensionalState.segment(3,3) = derivativeRotationMatrix * cartesianPosition + rotationMatrix * cartesianVelocity;

   Eigen::Vector6d normalizedState = gravitation::circular_restricted_three_body_problem::convertDimensionalCartesianStateToDimensionlessState(corotatingDimensionalState, gravitationalParameterPrimary, gravitationalParameterSecondary, distancePrimarySecondary);

   return normalizedState;
}



//! Function to directly setup CR3BP bodyMap
simulation_setup::NamedBodyMap setupBodyMapCR3BPBodyMap(
        double distancePrimarySecondary,
        std::string namePrimaryBody,
        std::string nameSecondaryBody,
        std::string nameBodyToPropagate,
        double initialTime
            ){

    spice_interface::loadStandardSpiceKernels( );


    double gravitationalParameterPrimary;
    double gravitationalParameterSecondary;

    // retrieve the gravitational parameter of the primary
    if (namePrimaryBody == "Earth" || namePrimaryBody == "Mars" || namePrimaryBody == "Moon"){
        std::shared_ptr< simulation_setup::GravityFieldSettings > gravityFieldSettings = simulation_setup::getDefaultGravityFieldSettings(namePrimaryBody, 0.0, 1.0);
        std::shared_ptr< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings > modelGravityFieldSettings =
                std::dynamic_pointer_cast< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        gravitationalParameterPrimary = modelGravityFieldSettings->getGravitationalParameter();
    }
    else{
        gravitationalParameterPrimary = spice_interface::getBodyGravitationalParameter(namePrimaryBody);
    }


    // retrieve the gravitational parameter of the secondary
    if (nameSecondaryBody == "Earth" || nameSecondaryBody == "Mars" || nameSecondaryBody == "Moon"){
        std::shared_ptr< simulation_setup::GravityFieldSettings > gravityFieldSettings = simulation_setup::getDefaultGravityFieldSettings(nameSecondaryBody, 0.0, 1.0);
        std::shared_ptr< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings > modelGravityFieldSettings =
                std::dynamic_pointer_cast< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        gravitationalParameterSecondary = modelGravityFieldSettings->getGravitationalParameter();
    }
    else{
        gravitationalParameterSecondary = spice_interface::getBodyGravitationalParameter(nameSecondaryBody);
    }



    double massParameter = gravitation::circular_restricted_three_body_problem::computeMassParameter(
        gravitationalParameterPrimary, gravitationalParameterSecondary );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "J2000";

    // Initial state for the primary
    Eigen::Vector6d initialStateInKeplerianElementsPrimary;
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::semiMajorAxisIndex) = massParameter * distancePrimarySecondary;
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::eccentricityIndex) = 0.0;
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::inclinationIndex) = 0.0;
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::argumentOfPeriapsisIndex) = 0.0;
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::longitudeOfAscendingNodeIndex) = 0.0;
    initialStateInKeplerianElementsPrimary(orbital_element_conversions::trueAnomalyIndex) = mathematical_constants::PI;


    // Initial state for the secondary
    Eigen::Vector6d initialStateInKeplerianElementsSecondary;
    initialStateInKeplerianElementsSecondary(orbital_element_conversions::semiMajorAxisIndex) = (1 - massParameter) * distancePrimarySecondary;
    initialStateInKeplerianElementsSecondary(orbital_element_conversions::eccentricityIndex) = 0.0;
    initialStateInKeplerianElementsSecondary(orbital_element_conversions::inclinationIndex) = 0.0;
    initialStateInKeplerianElementsSecondary(orbital_element_conversions::argumentOfPeriapsisIndex) = 0.0;
    initialStateInKeplerianElementsSecondary(orbital_element_conversions::longitudeOfAscendingNodeIndex) = 0.0;
    initialStateInKeplerianElementsSecondary(orbital_element_conversions::trueAnomalyIndex) = 0.0;


    // Create body objects.
    std::vector< std::string > bodiesToCreate;

    bodiesToCreate.push_back(namePrimaryBody);
    bodiesToCreate.push_back(nameSecondaryBody);

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
                simulation_setup::getDefaultBodySettings( bodiesToCreate);

    double sumGravitationalParameter = gravitationalParameterPrimary + gravitationalParameterSecondary;
    double distanceBarycenterPrimary = gravitationalParameterSecondary * distancePrimarySecondary / ( sumGravitationalParameter );
    double distanceBarycenterSecondary = distancePrimarySecondary - distanceBarycenterPrimary;
    double gravitationalParameterPrimaryTwoBodyProblem = std::pow(distanceBarycenterPrimary, 3) * sumGravitationalParameter / std::pow(distanceBarycenterPrimary + distanceBarycenterSecondary, 3);
    double gravitationalParameterSecondaryTwoBodyProblem = std::pow(distanceBarycenterSecondary, 3) * sumGravitationalParameter / std::pow(distanceBarycenterPrimary + distanceBarycenterSecondary, 3);

    // Primary ephemeris
    bodySettings[ namePrimaryBody ]->ephemerisSettings = std::make_shared< simulation_setup::KeplerEphemerisSettings >(
        initialStateInKeplerianElementsPrimary, initialTime, gravitationalParameterPrimaryTwoBodyProblem, "SSB", frameOrientation );

    // Secondary ephemeris
    bodySettings[ nameSecondaryBody ]->ephemerisSettings = std::make_shared< simulation_setup::KeplerEphemerisSettings >(
        initialStateInKeplerianElementsSecondary, initialTime, gravitationalParameterSecondaryTwoBodyProblem, "SSB", frameOrientation );



    for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
    {
        bodySettings[ bodiesToCreate.at( j ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( j ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }


    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );


    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ nameBodyToPropagate ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                std::shared_ptr< interpolators::OneDimensionalInterpolator
                < double, Eigen::Vector6d > >( ), "SSB", "J2000" ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    return bodyMap;

}


//! Function to directly setup CR3BP acceleration map
basic_astrodynamics::AccelerationMap setupAccelerationMapCR3BP(
        std::string namePrimaryBody,
        std::string nameSecondaryBody,
        std::string nameBodyToPropagate,
        std::vector< std::string > bodiesToPropagate,
        std::vector< std::string > centralBodies,
        simulation_setup::NamedBodyMap bodyMap
        ){

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
        double initialTime,
        double finalTime,
        Eigen::Vector6d initialState,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        basic_astrodynamics::AccelerationMap accelerationModelMap,
        std::vector< std::string > bodiesToPropagate,
        std::vector< std::string > centralBodies,
        simulation_setup::NamedBodyMap bodyMap,
        std::vector < std::string > bodiesCR3BP){


    using std::cout;
    cout.precision(15);


    // propagator definition (define new final propagation time longer than the expected final time for interpolation purposes)

    double initialIntegratorStepSize = integratorSettings->initialTimeStep_;
    double finalPropagationTime = finalTime + 3 * initialIntegratorStepSize;


    // Create barycenter object to retrieve the position of the primary and secondary.
    bodyMap[ "barycenter" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "barycenter" ]->setEphemeris(std::make_shared< ephemerides::ConstantEphemeris >(
                                              ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), "SSB", "J2000" ) );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    std::vector< std::shared_ptr< tudat::propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                                          tudat::propagators::relative_position_dependent_variable, bodiesCR3BP[0], "barycenter" ) );
    dependentVariablesList.push_back( std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                                          tudat::propagators::relative_velocity_dependent_variable, bodiesCR3BP[0], "barycenter" ) );
    dependentVariablesList.push_back( std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                                          tudat::propagators::relative_position_dependent_variable, bodiesCR3BP[1], "barycenter" ) );
    dependentVariablesList.push_back( std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                                          tudat::propagators::relative_velocity_dependent_variable, bodiesCR3BP[1], "barycenter" ) );
    std::shared_ptr< tudat::propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< tudat::propagators::DependentVariableSaveSettings >( dependentVariablesList );

    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double> > propagatorSettings =
        std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, initialState, finalPropagationTime, propagators::cowell, dependentVariablesToSave);



    // Create interpolator settings
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 );

    // Propagate the full problem
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator(bodyMap, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblem = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );



    // Interpolation of the full problem state history
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > interpolatorFullProblem =
            interpolators::createOneDimensionalInterpolator( stateHistoryFullProblem, interpolatorSettings );
    Eigen::Vector6d finalPropagatedStateFullProblem = interpolatorFullProblem->interpolate(finalTime);

    std::map< double, Eigen::VectorXd > stateHistoryBodies = dynamicsSimulator.getDependentVariableHistory();
    Eigen::VectorXd initialStateBodies = stateHistoryBodies[0];



    // CR3BP definition

    double gravitationalParameterPrimary;
    double gravitationalParameterSecondary;
    Eigen::Vector6d initialStatePrimary;
    Eigen::Vector6d initialStateSecondary;

    if (bodyMap[bodiesCR3BP[0]]->getBodyMass() > bodyMap[bodiesCR3BP[1]]->getBodyMass()){
        gravitationalParameterPrimary = bodyMap[bodiesCR3BP[0]]->getGravityFieldModel()->getGravitationalParameter();
        initialStatePrimary = initialStateBodies.segment(0,6);
        gravitationalParameterSecondary = bodyMap[bodiesCR3BP[1]]->getGravityFieldModel()->getGravitationalParameter();
        initialStateSecondary = initialStateBodies.segment(6,6);
    }

    else{
        gravitationalParameterPrimary = bodyMap[bodiesCR3BP[1]]->getGravityFieldModel()->getGravitationalParameter();
        initialStatePrimary = initialStateBodies.segment(6,6);
        gravitationalParameterSecondary = bodyMap[bodiesCR3BP[0]]->getGravityFieldModel()->getGravitationalParameter();
        initialStateSecondary = initialStateBodies.segment(0,6);
    }


    double massParameter = gravitation::circular_restricted_three_body_problem::computeMassParameter(gravitationalParameterPrimary, gravitationalParameterSecondary);
    double distanceBetweenPrimaries = std::sqrt( std::pow(initialStateSecondary[0] - initialStatePrimary[0], 2) + std::pow( initialStateSecondary[1] - initialStatePrimary[1], 2) + std::pow( initialStateSecondary[2] - initialStatePrimary[2], 2) );

    double normalizedInitialTime = gravitation::circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(initialTime, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);
    double normalizedFinalTime = gravitation::circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(finalTime, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);
    double normalizedFinalPropagationTime = gravitation::circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(finalPropagationTime, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);
    Eigen::Vector6d normalizedInitialState = convertCartesianToCorotatingNormalizedCoordinates(gravitationalParameterSecondary, gravitationalParameterPrimary, distanceBetweenPrimaries, initialState, initialTime);



    // CR3BP propagation
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > CR3BPintegratorSettings = integratorSettings;
    CR3BPintegratorSettings->initialTime_ = normalizedInitialTime;
    double initialTimeStepIntegrator = integratorSettings->initialTimeStep_;
    CR3BPintegratorSettings->initialTimeStep_ = gravitation::circular_restricted_three_body_problem::convertDimensionalTimeToDimensionlessTime(initialTimeStepIntegrator, gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries);

    std::map< double, Eigen::Vector6d > stateHistoryCR3BP = propagators::performCR3BPIntegration(CR3BPintegratorSettings, massParameter, normalizedInitialState, normalizedFinalPropagationTime );


    // Interpolation of the CR3BP state history
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > > interpolatorCR3BP =
            interpolators::createOneDimensionalInterpolator( stateHistoryCR3BP, interpolatorSettings );

    Eigen::Vector6d normalizedFinalPropagatedStateCR3BP = interpolatorCR3BP->interpolate(normalizedFinalTime);
    double normalizedFinalPropagationTimeCR3BP = normalizedFinalTime;



    // Transformation to inertial coordinates
    Eigen::Vector6d finalPropagatedStateCR3BP = convertCorotatingNormalizedToCartesianCoordinates(gravitationalParameterPrimary, gravitationalParameterSecondary, distanceBetweenPrimaries, normalizedFinalPropagatedStateCR3BP, normalizedFinalPropagationTimeCR3BP);


    // Difference between the two propagated states at finalTime
    Eigen::Vector6d finalStateDifference = finalPropagatedStateCR3BP - finalPropagatedStateFullProblem;


    return finalStateDifference;

}






int main(){

    std::cout.precision(20);

    double initialTime = 0.0;
    double finalTime = 120000000.0;

    std::vector < std::string > bodiesCR3BP;
    bodiesCR3BP.push_back("Sun");
    bodiesCR3BP.push_back("Earth");

    simulation_setup::NamedBodyMap bodyMap = setupBodyMapCR3BPBodyMap(physical_constants::ASTRONOMICAL_UNIT, "Sun", "Earth", "spacecraft", 0.0);

    // Spacecraft properties
    bodyMap[ "spacecraft" ]->setConstantBodyMass(100.0);

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

    basic_astrodynamics::AccelerationMap accelerationModelMap = setupAccelerationMapCR3BP("Sun", "Earth", "spacecraft", bodiesToPropagate, centralBodies, bodyMap);



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


