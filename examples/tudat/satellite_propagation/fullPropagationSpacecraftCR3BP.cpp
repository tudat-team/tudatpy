/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <Eigen/Core>

#include <tudat/astro/basic/physicalConstants.h>
#include <tudat/basics/testMacros.h>
#include <tudat/math/basic/mathematicalConstants.h>
#include "tudat/astro/basic/unitConversions.h"
#include <tudat/astro/basic/orbitalElementConversions.h>
#include "tudat/astro/gravitation/librationPoint.h"
#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "tudat/simulation/propagation/propagationCR3BPFullProblem.h"

#include <tudat/simulation/simulation.h>
#include <tudat/io/basicInputOutput.h>
#include <tudat/io/applicationOutput.h>
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"


int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;
    using namespace tudat::simulation_setup;

    spice_interface::loadStandardSpiceKernels( );


    // Global characteristics of the problem
    double distanceSunJupiter = 778.0e9;

    // Initialise the spacecraft state (B. Taylor, D. (1981). Horseshoe periodic orbits in the restricted problem of three bodies
    // for a sun-Jupiter mass ratio. Astronomy and Astrophysics. 103. 288-294.)
    Eigen::Vector6d initialState = Eigen::Vector6d::Zero();
    initialState[0] = - 7.992e11;
    initialState[4] =  -1.29e4;

    // Create integrator settings.
    double initialTime = 0.0;
    const double fixedStepSize = 100000.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize );

    // Create body map.
    std::vector < std::string > bodiesCR3BP;
    bodiesCR3BP.push_back( "Sun" );
    bodiesCR3BP.push_back( "Jupiter" );

    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Spacecraft" );
    centralBodies.push_back( "SSB" );

    // Define final time for the propagation.
    double gravitationalParameterSun = createGravityFieldModel(
                getDefaultGravityFieldSettings(
                    "Sun", TUDAT_NAN, TUDAT_NAN ), "Sun" )->getGravitationalParameter( );
    double gravitationalParameterJupiter = createGravityFieldModel(
                getDefaultGravityFieldSettings(
                    "Jupiter", TUDAT_NAN, TUDAT_NAN ), "Jupiter" )->getGravitationalParameter( );
    double finalTime = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(
                29.2386 * ( 2.0 * mathematical_constants::PI ), gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter);

    NamedBodyMap idealBodyMap = propagators::setupBodyMapCR3BP(
                distanceSunJupiter, "Sun", "Jupiter", "Spacecraft" );

    std::map< double, Eigen::Vector6d> fullPropagation;
    std::map< double, Eigen::Vector6d> cr3bpPropagation;

    /// Ideal case: full dynamics problem with the CR3BP assumtions
    {


        // Create acceleration map.
        basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapCR3BP(
                    "Sun", "Jupiter", bodiesToPropagate.at( 0 ), centralBodies.at( 0 ), idealBodyMap );

        // Calculate the difference between CR3BP and full problem.
        propagators::propagateCR3BPAndFullDynamicsProblem(
                    initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                    bodiesToPropagate, centralBodies, idealBodyMap, bodiesCR3BP, fullPropagation,
                    cr3bpPropagation );

        Eigen::Vector6d stateDifference =
                fullPropagation.rbegin( )->second - cr3bpPropagation.rbegin( )->second;

        std::cout << "state difference at final time: " << stateDifference << std::endl;
    }


    std::map< double, Eigen::Vector6d> fullPropagationPerturbedCase;
    std::map< double, Eigen::Vector6d> cr3bpPropagationPerturbedCase;

    /// Perturbed case
    {

        std::string frameOrigin = "SSB";
        std::string frameOrientation = "ECLIPJ2000";

        NamedBodyMap perturbedBodyMap;


        std::vector< std::string > additionalBodies = { "Earth", "Mars", "Venus", "Saturn" };
        for( unsigned int i = 0; i < additionalBodies.size( ); i++ )
        {

            perturbedBodyMap[ additionalBodies.at( i ) ] = std::make_shared< Body >( );
            perturbedBodyMap[ additionalBodies.at( i ) ]->setEphemeris(
                        std::make_shared< ephemerides::ApproximatePlanetPositions>(
                                additionalBodies.at( i ) ) );
            perturbedBodyMap[ additionalBodies.at( i ) ]->setGravityFieldModel(
                        createGravityFieldModel(
                            std::make_shared< CentralGravityFieldSettings >(
                                spice_interface::getBodyGravitationalParameter(
                                    additionalBodies.at( i ) ) ), additionalBodies.at( i ) ) );
        }

        perturbedBodyMap[ "Sun" ] = idealBodyMap[ "Sun" ];
        perturbedBodyMap[ "Jupiter" ] = idealBodyMap[ "Jupiter" ];


        // Create the body to be propagated.
        perturbedBodyMap[ "Spacecraft" ] = std::make_shared< Body >( );
        perturbedBodyMap[ "Spacecraft" ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                            < double, Eigen::Vector6d > >( ), "SSB", frameOrientation ) );

        setGlobalFrameBodyEphemerides( perturbedBodyMap, frameOrigin, frameOrientation );


        // Set of accelerations experienced by the spacecraft.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > bodyToPropagateAccelerations;
        bodyToPropagateAccelerations["Sun"].push_back(std::make_shared< AccelerationSettings >(
                                                          basic_astrodynamics::central_gravity ) );
        bodyToPropagateAccelerations["Jupiter"].push_back(std::make_shared< AccelerationSettings >(
                                                              basic_astrodynamics::central_gravity ) );
        bodyToPropagateAccelerations["Earth"].push_back(std::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::central_gravity ) );
        bodyToPropagateAccelerations["Mars"].push_back(std::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
        bodyToPropagateAccelerations["Venus"].push_back(std::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::central_gravity ) );
        bodyToPropagateAccelerations["Saturn"].push_back(std::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );

        SelectedAccelerationMap accelerationMap;
        accelerationMap[ "Spacecraft" ] = bodyToPropagateAccelerations;


        // Create the acceleration map.
        basic_astrodynamics::AccelerationMap accelerationModelMapPerturbedCase = createAccelerationModelsMap(
                    perturbedBodyMap, accelerationMap, bodiesToPropagate, centralBodies );


        // Calculate the difference between CR3BP and full problem.
        propagators::propagateCR3BPAndFullDynamicsProblem(
                    initialTime, finalTime, initialState, integratorSettings,
                    accelerationModelMapPerturbedCase,
                    bodiesToPropagate, centralBodies, perturbedBodyMap, bodiesCR3BP,
                    fullPropagationPerturbedCase,
                    cr3bpPropagationPerturbedCase );

        Eigen::Vector6d stateDifferencePerturbedCase =
                fullPropagationPerturbedCase.rbegin( )->second - cr3bpPropagationPerturbedCase.rbegin( )->second;

        std::cout << "state difference at final time for the perturbed case: " << stateDifferencePerturbedCase << std::endl;

    }

    /// Ouputs

    // Outputs for the ideal case
    {
        std::map< double, Eigen::Vector6d > fullPropagationNormalisedCoRotatingFrame;
        for( std::map< double, Eigen::Vector6d >::iterator itr = fullPropagation.begin( );
             itr != fullPropagation.end( ); itr++ ){
            fullPropagationNormalisedCoRotatingFrame[ itr->first ] = tudat::circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
                        gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
        }

        std::map< double, Eigen::Vector6d > cr3bpNormalisedCoRotatingFrame;
        for( std::map< double, Eigen::Vector6d >::iterator itr = cr3bpPropagation.begin( );
             itr != cr3bpPropagation.end( ); itr++ ){
            cr3bpNormalisedCoRotatingFrame[ itr->first ] = tudat::circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
                        gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
        }


        input_output::writeDataMapToTextFile( fullPropagation,
                                              "fullProblemPropagation.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( fullPropagationNormalisedCoRotatingFrame,
                                              "fullProblemPropagationNormalisedCoRotatingFrame.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( cr3bpPropagation,
                                              "CR3BPsolution.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( cr3bpNormalisedCoRotatingFrame,
                                              "CR3BPnormalisedCoRotatingFrame.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );
    }

    // Outputs for the perturbed case
    {
        std::map< double, Eigen::Vector6d > fullPropagationNormalisedCoRotatingFramePerturbedCase;
        for( std::map< double, Eigen::Vector6d >::iterator itr = fullPropagationPerturbedCase.begin( );
             itr != fullPropagationPerturbedCase.end( ); itr++ ){
            fullPropagationNormalisedCoRotatingFramePerturbedCase[ itr->first ] = tudat::circular_restricted_three_body_problem::
                    convertCartesianToCorotatingNormalizedCoordinates(gravitationalParameterSun, gravitationalParameterJupiter,
                                                                      distanceSunJupiter, itr->second, itr->first);
        }

        std::map< double, Eigen::Vector6d > cr3bpNormalisedCoRotatingFramePerturbedCase;
        for( std::map< double, Eigen::Vector6d >::iterator itr = cr3bpPropagationPerturbedCase.begin( );
             itr != cr3bpPropagationPerturbedCase.end( ); itr++ ){
            cr3bpNormalisedCoRotatingFramePerturbedCase[ itr->first ] = tudat::circular_restricted_three_body_problem::
                    convertCartesianToCorotatingNormalizedCoordinates(
                        gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
        }


        input_output::writeDataMapToTextFile( fullPropagationPerturbedCase,
                                              "fullProblemPropagationPerturbedCase.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( fullPropagationNormalisedCoRotatingFramePerturbedCase,
                                              "fullProblemPropagationNormalisedCoRotatingFramePerturbedCase.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( cr3bpPropagationPerturbedCase,
                                              "CR3BPsolutionPerturbedCase.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( cr3bpNormalisedCoRotatingFramePerturbedCase,
                                              "CR3BPnormalisedCoRotatingFramePerturbedCase.dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
