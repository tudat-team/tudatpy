/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <math.h>
#include <iostream>

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/LowThrustTrajectories/simsFlanagan.h"
#include "tudat/astro/LowThrustTrajectories/simsFlanaganModel.h"
#include "tudat/astro/LowThrustTrajectories/simsFlanaganOptimisationSetup.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "pagmo/algorithms/de1220.hpp"
#include "Problems/applicationOutput.h"
#include "tudat/astro/basic/celestialBodyConstants.h"

int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;
    using namespace tudat::simulation_setup;
    using namespace tudat::shape_based_methods;
    using namespace tudat::low_thrust_trajectories;
    using namespace shape_based_methods;


    spice_interface::loadStandardSpiceKernels( );


    double julianDate = 9264.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1000.0 * physical_constants::JULIAN_DAY;
    int numberOfRevolutions = 2;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    double maximumThrust = 5.0;
    double specificImpulse = 3000.0;
    double mass = 2800.0;
    int numberSegments = 500;


    std::string bodyToPropagate = "Borzi";
    std::string centralBody = "Sun";


    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ centralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ centralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ centralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ bodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( bodyToPropagate )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );



    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;


    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Get recommended base functions for the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Get recommended base functions for normal velocity composition function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Get recommended base functions for axial velocity composition function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    std::shared_ptr< shape_based_methods::HodographicShaping > hodographicShaping =
            std::make_shared< shape_based_methods::HodographicShaping >(
                cartesianStateDepartureBody, cartesianStateArrivalBody, timeOfFlight, numberOfRevolutions,
                bodyMap, bodyToPropagate, centralBody, radialVelocityFunctionComponents, normalVelocityFunctionComponents,
                axialVelocityFunctionComponents, freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );




    // Save results
    int numberOfSteps = 10000;
    double stepSize = timeOfFlight / static_cast< double >( numberOfSteps );

    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double time )
    {
        return specificImpulse;
    };

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize );

    std::vector< double > epochsToSaveResults;
    for ( int i = 0 ; i <= numberOfSteps ; i++ )
    {
        epochsToSaveResults.push_back( i * stepSize );
    }

    std::map< double, Eigen::Vector6d > hodographicShapingTrajectory;
    std::map< double, Eigen::VectorXd > hodographicShapingMassProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustAcceleration;
    hodographicShaping->getTrajectory( epochsToSaveResults, hodographicShapingTrajectory );
    hodographicShaping->getMassProfile( epochsToSaveResults, hodographicShapingMassProfile, specificImpulseFunction, integratorSettings );
    hodographicShaping->getThrustProfile( epochsToSaveResults, hodographicShapingThrustProfile, specificImpulseFunction, integratorSettings );
    hodographicShaping->getThrustAccelerationProfile( epochsToSaveResults, hodographicShapingThrustAcceleration,
                                                      specificImpulseFunction, integratorSettings );




    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    std::shared_ptr< ShapeBasedMethod > shapeBased = std::dynamic_pointer_cast< ShapeBasedMethod >( hodographicShaping );

    std::function< Eigen::Vector3d( const double ) > initialGuessThrustFromShaping =
            getInitialGuessFunctionFromShaping( shapeBased, numberSegments, timeOfFlight, specificImpulseFunction, integratorSettings );

    std::vector< double > initialGuessVector = convertToSimsFlanaganThrustModel( initialGuessThrustFromShaping, maximumThrust,
                                                                                 timeOfFlight, ( numberSegments + 1 ) / 2, numberSegments / 2 );

    std::map< double, Eigen::Vector3d > initialGuessThrustProfile;
    std::map< double, Eigen::Vector3d > initialGuessThrustAccelerationProfile;
    std::map< double, Eigen::Vector1d > initialGuessMassProfile;

    for ( unsigned int i = 0 ; i < epochsToSaveResults.size( ) ; i++ )
    {
        double currentTime = epochsToSaveResults[ i ];
        initialGuessThrustProfile[ currentTime ] =  initialGuessThrustFromShaping( currentTime );

        if ( i == 0 )
        {
            initialGuessMassProfile[ currentTime ] = ( Eigen::Vector1d( ) << mass ).finished( );
        }
        else
        {
            initialGuessMassProfile[ currentTime ] = ( Eigen::Vector1d( ) << initialGuessMassProfile[ epochsToSaveResults[ i - 1 ] ][ 0 ]
                    - initialGuessThrustProfile[ currentTime ].norm( ) /
                    ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) * stepSize ).finished( );
        }

        initialGuessThrustAccelerationProfile[ currentTime ] = initialGuessThrustProfile[ currentTime ] / initialGuessMassProfile[ currentTime ][ 0 ];

    }


    bodyMap[ bodyToPropagate ]->setConstantBodyMass( mass );

    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    std::shared_ptr< OptimisationSettings > optimisationSettings = std::make_shared< OptimisationSettings >(
                optimisationAlgorithm, 10, 1024, 1.0e-6, std::make_pair( initialGuessVector, 0.3 ) );

    SimsFlanagan simsFlanagan = SimsFlanagan(
                cartesianStateDepartureBody, cartesianStateArrivalBody, maximumThrust, specificImpulseFunction, numberSegments,
                timeOfFlight, bodyMap, bodyToPropagate, centralBody, optimisationSettings );


    std::map< double, Eigen::Vector6d > SimsFlanaganTrajectory;
    std::map< double, Eigen::VectorXd > SimsFlanaganMassProfile;
    std::map< double, Eigen::VectorXd > SimsFlanaganThrustProfile;
    std::map< double, Eigen::VectorXd > SimsFlanaganThrustAcceleration;
    simsFlanagan.getTrajectory( epochsToSaveResults, SimsFlanaganTrajectory );
    simsFlanagan.getMassProfile( epochsToSaveResults, SimsFlanaganMassProfile, specificImpulseFunction, integratorSettings );
    simsFlanagan.getThrustProfile( epochsToSaveResults, SimsFlanaganThrustProfile, specificImpulseFunction, integratorSettings );
    simsFlanagan.getThrustAccelerationProfile( epochsToSaveResults, SimsFlanaganThrustAcceleration, specificImpulseFunction, integratorSettings );


    input_output::writeDataMapToTextFile( hodographicShapingTrajectory,
                                          "hodographicShapingTrajectory.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingMassProfile,
                                          "hodographicShapingMassProfile.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingThrustProfile,
                                          "hodographicShapingThrustProfile.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingThrustAcceleration,
                                          "hodographicShapingThrustAcceleration.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( SimsFlanaganTrajectory,
                                          "SimsFlanaganTrajectory.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( SimsFlanaganMassProfile,
                                          "SimsFlanaganMassProfile.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( SimsFlanaganThrustProfile,
                                          "SimsFlanaganThrustProfile.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( SimsFlanaganThrustAcceleration,
                                          "SimsFlanaganThrustAcceleration.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( initialGuessThrustProfile,
                                          "initialGuessThrustProfile.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( initialGuessThrustAccelerationProfile,
                                          "initialGuessThrustAccelerationProfile.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( initialGuessMassProfile,
                                          "initialGuessMassProfile.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    std::cout << "DELTAV SIMS FLANAGAN: " << simsFlanagan.computeDeltaV( ) << "\n\n";
    std::cout << "DELTAV SHAPE BASED: " << hodographicShaping->computeDeltaV( ) << "\n\n";


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
