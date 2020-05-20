/*    Copyright (c) 2010-2019, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include <iostream>
#include <fstream>
#include <functional>

#include <boost/filesystem.hpp>
#include "Problems/applicationOutput.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"

#include "tudat/simulation/simulation.h"
#include "tudat/astro/LowThrustTrajectories/lowThrustOptimisationSetup.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/sphericalShaping.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "tudat/astro/LowThrustTrajectories/lowThrustLegSettings.h"
#include "tudat/astro/LowThrustTrajectories/lowThrustLeg.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/hodographicShapingOptimisationSetup.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/getRecommendedBaseFunctionsHodographicShaping.h"
#include "tudat/simulation/optimisationSettings.h"

using namespace tudat;
using namespace tudat::shape_based_methods;
using namespace tudat::numerical_integrators;
using namespace tudat::simulation_setup;

std::vector< std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > > getShapingBasisFunctions(
        const double timeOfFlight, const int numberOfRevolutions )
{
    Eigen::VectorXd dummyVector;

    // Get recommended base functions for the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    getRecommendedRadialVelocityBaseFunctions(
                radialVelocityFunctionComponents, dummyVector, timeOfFlight );

    // Get recommended base functions for the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    getRecommendedNormalAxialBaseFunctions(
                normalVelocityFunctionComponents, dummyVector, timeOfFlight );

    // Get recommended base functions for the axial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    getRecommendedAxialVelocityBaseFunctions(
                axialVelocityFunctionComponents, dummyVector, timeOfFlight, numberOfRevolutions );

    {
        double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
        double scaleFactor = 1.0 / timeOfFlight;

        std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
                std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );
        std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
                std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );

        // Add two additional base functions
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );

    }

    {
        double scaleFactor = 1.0 / timeOfFlight;

        // Create base function settings for the components of the axial velocity composite function.
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 3.0, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 4.0, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 5.0, scaleFactor );


        // Set components for the axial velocity function.
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, firstAxialVelocityBaseFunctionSettings ) );
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondAxialVelocityBaseFunctionSettings ) );
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, thirdAxialVelocityBaseFunctionSettings ) );

    }

    return { radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents };
}

//! Execute  main
int main( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 123 );

    tudat::spice_interface::loadStandardSpiceKernels( );

    // Ephemeris functions of bodies.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    std::function< Eigen::Vector6d( const double ) > departureStateFunction = [ = ]( const double currentTime )
    { return pointerToDepartureBodyEphemeris->getCartesianState( currentTime ); };
    std::function< Eigen::Vector6d( const double ) > arrivalStateFunction = [ = ]( const double currentTime )
    { return pointerToArrivalBodyEphemeris->getCartesianState( currentTime ); };


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////        GRID SEARCH FOR HODOGRAPHIC SHAPING LOWEST-ORDER SOLUTION            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define bounds for departure date and time-of-flight.
    std::pair< double, double > departureTimeBounds =
            std::make_pair( 7304.5 * physical_constants::JULIAN_DAY, 13225.5 * physical_constants::JULIAN_DAY  );
    std::pair< double, double > timeOfFlightBounds =
            std::make_pair( 300.0 * physical_constants::JULIAN_DAY, 2000.0 * physical_constants::JULIAN_DAY );

    // Define lower and upper bounds for the radial velocity free coefficients.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 7, 0.0 ) );
    bounds[ 0 ][ 0 ] =  departureTimeBounds.first;
    bounds[ 1 ][ 0 ] = departureTimeBounds.second;
    bounds[ 0 ][ 1 ] = timeOfFlightBounds.first;
    bounds[ 1 ][ 1 ] = timeOfFlightBounds.second;
    bounds[ 0 ][ 2 ] = - 5000.0;
    bounds[ 1 ][ 2 ] = 5000.0;
    bounds[ 0 ][ 3 ] = -5000;
    bounds[ 1 ][ 3 ] = 5000.0;
    bounds[ 0 ][ 4 ] = - 5000.0;
    bounds[ 1 ][ 4 ] = 5000.0;
    bounds[ 0 ][ 5 ] = -5000;
    bounds[ 1 ][ 5 ] = 5000.0;
    bounds[ 0 ][ 6 ] = -5000;
    bounds[ 1 ][ 6 ] = 5000.0;

    for( int useMultiObjective = 0; useMultiObjective < 2; useMultiObjective++ )
    {
        for( int revolutions = 2; revolutions < 3; revolutions++ )
        {
            // Create object to compute the problem fitness
            problem prob{ HodographicShapingOptimisationProblem(
                            departureStateFunction, arrivalStateFunction, spice_interface::getBodyGravitationalParameter( "Sun" ),
                            revolutions, std::bind( &getShapingBasisFunctions, std::placeholders::_1, revolutions ), bounds,
                            useMultiObjective, 2000.0 ) };

            //sade, gaco, sga, de
            algorithm algo;
            if( !useMultiObjective )
            {
                algo = algorithm{ simulated_annealing( ) };
            }
            else
            {
                algo = algorithm{ nsga2( ) };
            }

            // Create an island with 1000 individuals
            island isl{algo, prob, 2000 };

            // Evolve for 512 generations
            for( int i = 0 ; i < 101; i++ )
            {
                isl.evolve( );
                while( isl.status( ) != pagmo::evolve_status::idle &&
                       isl.status( ) != pagmo::evolve_status::idle_error )
                {
                    isl.wait( );
                }

                if( i % 25 == 0 )
                {
                    if( !useMultiObjective )
                    {
                        std::cout<<"Iteration: "<<" "<<i<<"; Best Delta V: "<<isl.get_population( ).champion_f( ).at( 0 )<<std::endl;

                        printPopulationToFile( isl.get_population( ).get_x( ), "hodograph_single_objective_" + std::to_string( i / 25 ), false );
                        printPopulationToFile( isl.get_population( ).get_f( ), "hodograph_single_objective_" + std::to_string( i / 25 ), true );
                    }
                    else
                    {
                        std::cout<<"Iteration: "<<" "<<i<<std::endl;

                        printPopulationToFile( isl.get_population( ).get_x( ), "hodograph_multi_objective_" + std::to_string( i / 25 ), false );
                        printPopulationToFile( isl.get_population( ).get_f( ), "hodograph_multi_objective_" + std::to_string( i / 25 ), true );
                    }
                }

            }


            if( !useMultiObjective )
            {
                std::cout<<"Final best Delta V: "<<isl.get_population( ).champion_f( ).at( 0 )<<std::endl;

                std::vector< double > bestPopulation = isl.get_population( ).champion_x( );

                Eigen::VectorXd radialFreeParameters = Eigen::VectorXd::Zero( 2 );
                Eigen::VectorXd normalFreeParameters = Eigen::VectorXd::Zero( 0 );
                Eigen::VectorXd axialFreeParameters = Eigen::VectorXd::Zero( 3 );

                radialFreeParameters << bestPopulation.at( 2 ), bestPopulation.at( 3 );
                axialFreeParameters << bestPopulation.at( 4 ), bestPopulation.at( 5 ), bestPopulation.at( 6 );

                double timeOfFlight = bestPopulation.at( 1 );
                double derpartureTime = bestPopulation.at( 0 );
                double arrivalTime = derpartureTime + timeOfFlight;

                double initialMass = 2000.0;
                double specificImpulse = 3000.0;


                std::vector< std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > > shapingFunctions = getShapingBasisFunctions(
                            timeOfFlight, 2 );

                std::shared_ptr< HodographicShaping > hodographicShaping =
                        std::make_shared< HodographicShaping >(
                            departureStateFunction( derpartureTime ), arrivalStateFunction( arrivalTime ), timeOfFlight,
                            spice_interface::getBodyGravitationalParameter( "Sun" ), 2,
                            shapingFunctions.at( 0 ), shapingFunctions.at( 1 ), shapingFunctions.at( 2 ),
                            radialFreeParameters, normalFreeParameters, axialFreeParameters, initialMass );

                // Save results
                int numberOfSteps = 1000;
                double stepSize = timeOfFlight / static_cast< double >( numberOfSteps );
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

                hodographicShaping->getTrajectory(
                            epochsToSaveResults, hodographicShapingTrajectory );
                hodographicShaping->getMassProfile(
                            epochsToSaveResults, hodographicShapingMassProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
                hodographicShaping->getThrustForceProfile(
                            epochsToSaveResults, hodographicShapingThrustProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
                hodographicShaping->getCylindricalThrustAccelerationProfile(
                            epochsToSaveResults, hodographicShapingThrustAcceleration );

                input_output::writeDataMapToTextFile(
                            hodographicShapingTrajectory, "hodographicShapingOptimalTrajectory.dat", tudat_pagmo_applications::getOutputPath( ) );

                input_output::writeDataMapToTextFile(
                            hodographicShapingMassProfile, "hodographicShapingOptimalMassProfile.dat", tudat_pagmo_applications::getOutputPath( ) );

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustProfile, "hodographicShapingOptimalThrustProfile.dat", tudat_pagmo_applications::getOutputPath( ) );

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustAcceleration, "hodographicShapingOptimalThrustAcceleration.dat", tudat_pagmo_applications::getOutputPath( ) );
            }
        }
    }
}
