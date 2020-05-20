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

using namespace tudat;

//! Execute  main
int main( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 123 );

    tudat::spice_interface::loadStandardSpiceKernels( );

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
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
            std::make_pair( 7304.5 * physical_constants::JULIAN_DAY, 10225.5 * physical_constants::JULIAN_DAY  );
    std::pair< double, double > timeOfFlightBounds =
            std::make_pair( 500.0 * physical_constants::JULIAN_DAY, 2000.0 * physical_constants::JULIAN_DAY );

    // Initialize free coefficients vectors
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    std::map< int, Eigen::Vector4d > hodographicShapingResultsLowOrder;

    int numberCases = 0;

    // for-loop parsing the time-of-flight values, ranging from 500 to 2000 days, with a time-step of 5 days.
    for ( int i = 0 ; i <= ( timeOfFlightBounds.second - timeOfFlightBounds.first ) / ( 5.0 * physical_constants::JULIAN_DAY ) ; i++  )
    {
        double currentTOF = timeOfFlightBounds.first + i * 5.0 * physical_constants::JULIAN_DAY;

        // Get recommended base functions for the radial velocity composite function.
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
        shape_based_methods::getRecommendedRadialVelocityBaseFunctions(
                    radialVelocityFunctionComponents, freeCoefficientsRadialVelocityFunction, currentTOF );

        // Get recommended base functions for the normal velocity composite function.
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
        shape_based_methods::getRecommendedNormalAxialBaseFunctions(
                    normalVelocityFunctionComponents, freeCoefficientsNormalVelocityFunction, currentTOF );

        // for-loop parsing the departure date values, ranging from 7304 MJD to 10225 MJD (with 401 steps)
        for ( int j = 0 ; j <= 400; j++ )
        {
            double currentDepartureDate = departureTimeBounds.first + j * ( departureTimeBounds.second - departureTimeBounds.first ) / 400.0;

            // Compute states at departure and arrival.
            Eigen::Vector6d cartesianStateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( currentDepartureDate );
            Eigen::Vector6d cartesianStateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( currentDepartureDate + currentTOF );

            int bestNumberOfRevolutions;
            double currentBestDeltaV;

            // Parse shaped trajectories with numbers of revolutions between 0 and 5.
            for ( int currentNumberOfRevolutions = 0 ; currentNumberOfRevolutions <= 5 ; currentNumberOfRevolutions++ )
            {
                // Get recommended base functions for the axial velocity composite function.
                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
                shape_based_methods::getRecommendedAxialVelocityBaseFunctions(
                            axialVelocityFunctionComponents, freeCoefficientsAxialVelocityFunction, currentTOF,
                            currentNumberOfRevolutions );

                // Create hodographically shaped trajectory.
                tudat::shape_based_methods::HodographicShaping hodographicShaping = shape_based_methods::HodographicShaping(
                            cartesianStateAtDeparture, cartesianStateAtArrival, currentTOF,
                            spice_interface::getBodyGravitationalParameter( "Sun" ), currentNumberOfRevolutions,
                            radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                            freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
                            freeCoefficientsAxialVelocityFunction );

                // Save trajectory with the lowest deltaV.
                if ( currentNumberOfRevolutions == 0 )
                {
                    bestNumberOfRevolutions = 0;
                    currentBestDeltaV = hodographicShaping.computeDeltaV( );
                }
                else
                {
                    if ( hodographicShaping.computeDeltaV( ) < currentBestDeltaV )
                    {
                        currentBestDeltaV = hodographicShaping.computeDeltaV( );
                        bestNumberOfRevolutions = currentNumberOfRevolutions;
                    }
                }
            }

            // Save results.
            Eigen::Vector4d outputVector =
                    ( Eigen::Vector4d( ) << currentTOF / physical_constants::JULIAN_DAY,
                      currentDepartureDate / physical_constants::JULIAN_DAY, currentBestDeltaV, bestNumberOfRevolutions ).finished( );
            numberCases++;
            hodographicShapingResultsLowOrder[ numberCases ] = outputVector;

        }
    }

    input_output::writeDataMapToTextFile( hodographicShapingResultsLowOrder,
                                          "hodographicShapingLowOrder.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////        RESTRICTED GRID SEARCH FOR HODOGRAPHIC SHAPING HIGH-ORDER SOLUTION             ///////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    numberCases = 0;

    // Define lower and upper bounds for the radial velocity free coefficients.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 2, 0.0 ) );
    bounds[ 0 ][ 0 ] = - 600.0;
    bounds[ 1 ][ 0 ] = 800.0;
    bounds[ 0 ][ 1 ] = 0.0;
    bounds[ 1 ][ 1 ] = 1500.0;

    // Set fixed number of revolutions.
    int numberOfRevolutions = 1;

    std::map< int, Eigen::Vector4d > hodographicShapingResultsHigherOrder;
    std::map< int, Eigen::Vector4d > hodographicShapingResultsLowOrderOneRevolution;

    // for-loop parsing the time-of-flight values, ranging from 500 to 900 days, with a time-step of 20 days.
    for ( int i = 0 ; i <= ( 900.0 * physical_constants::JULIAN_DAY - timeOfFlightBounds.first ) / ( 20 * physical_constants::JULIAN_DAY ) ; i++  )
    {
        double currentTOF = timeOfFlightBounds.first + i * 20.0 * physical_constants::JULIAN_DAY;

        double frequency = 2.0 * mathematical_constants::PI / currentTOF;
        double scaleFactor = 1.0 / currentTOF;

        // Define settings for the two additional base functions for the radial velocity composite function.
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );

        // Get recommended base functions for the radial velocity composite function, and add two additional base functions
        // (introducing two degrees of freedom in the trajectory design problem).
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
        shape_based_methods::getRecommendedRadialVelocityBaseFunctions( radialVelocityFunctionComponents, freeCoefficientsRadialVelocityFunction, currentTOF );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );

        // Get recommended base functions for the normal velocity composite function.
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
        shape_based_methods::getRecommendedNormalAxialBaseFunctions( normalVelocityFunctionComponents, freeCoefficientsNormalVelocityFunction, currentTOF );


        // for-loop parsing departure dates ranging from 7304 MJD to 7379 MJD (with a time-step of 15 days).
        for ( int j = 0 ; j <= ( 7379.5 * physical_constants::JULIAN_DAY - departureTimeBounds.first ) / ( 15 * physical_constants::JULIAN_DAY ); j++ )
        {
            double currentDepartureDate = departureTimeBounds.first +
                    j * 15.0 * physical_constants::JULIAN_DAY;

            // Compute states at departure and arrival.
            Eigen::Vector6d cartesianStateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( currentDepartureDate );
            Eigen::Vector6d cartesianStateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( currentDepartureDate + currentTOF );


            // Get recommended base functions for the axial velocity composite function.
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
            shape_based_methods::getRecommendedAxialVelocityBaseFunctions( axialVelocityFunctionComponents, freeCoefficientsAxialVelocityFunction,
                                                                           currentTOF, numberOfRevolutions );


            // Create hodographic shaping optimisation problem.
            problem prob{ shape_based_methods::FixedTimeHodographicShapingOptimisationProblem(
                            cartesianStateAtDeparture, cartesianStateAtArrival, currentTOF,
                            spice_interface::getBodyGravitationalParameter( "Sun" ), numberOfRevolutions,
                            radialVelocityFunctionComponents,
                            normalVelocityFunctionComponents, axialVelocityFunctionComponents, bounds ) };

            // Perform optimisation.
            algorithm algo{ pagmo::sga( ) };

            // Create an island with 1024 individuals
            island isl{ algo, prob, 25 };

            // Evolve for 100 generations
            for( int i = 0 ; i < 10; i++ )
            {
                isl.evolve( );
                while( isl.status( ) != pagmo::evolve_status::idle &&
                       isl.status( ) != pagmo::evolve_status::idle_error )
                {
                    isl.wait( );
                }
                isl.wait_check( ); // Raises errors

            }

            // Save high-order shaping solution.
            double currentBestDeltaV = isl.get_population( ).champion_f( )[ 0 ];

            Eigen::Vector4d outputVector = ( Eigen::Vector4d( ) << currentTOF / physical_constants::JULIAN_DAY,
                                             currentDepartureDate / physical_constants::JULIAN_DAY, currentBestDeltaV, 1 ).finished( );

            hodographicShapingResultsHigherOrder[ numberCases ] = outputVector;


            // Compute corresponding low-order hodographic shaping solution.
            Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
            Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );
            Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

            // Compute low-order hodographically shaped trajectory (number of revolutions set to 1).
            tudat::shape_based_methods::HodographicShaping hodographicShapingLowOrderOneRevolution = shape_based_methods::HodographicShaping(
                        cartesianStateAtDeparture, cartesianStateAtArrival, currentTOF,
                        spice_interface::getBodyGravitationalParameter( "Sun" ), numberOfRevolutions,
                        radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                        freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction );

            // Save low-order shaping solution.
            outputVector = ( Eigen::Vector4d( ) << currentTOF / physical_constants::JULIAN_DAY,
                             currentDepartureDate / physical_constants::JULIAN_DAY, hodographicShapingLowOrderOneRevolution.computeDeltaV( ), 1 ).finished( );
            hodographicShapingResultsLowOrderOneRevolution[ numberCases ] = outputVector;

            numberCases++;

        }
    }

    input_output::writeDataMapToTextFile( hodographicShapingResultsLowOrderOneRevolution,
                                          "hodographicShapingLowOrderOneRevolution.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingResultsHigherOrder,
                                          "hodographicShapingHigherOrder.dat",
                                          tudat_pagmo_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    return 0;

}
