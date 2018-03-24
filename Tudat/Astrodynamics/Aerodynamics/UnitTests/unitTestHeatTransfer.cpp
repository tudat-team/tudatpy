/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/equilibriumWallTemperature.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"

namespace tudat
{
namespace unit_tests
{

using namespace aerodynamics;

BOOST_AUTO_TEST_SUITE( test_fay_riddell_heat_flux )

//! Test if equilibrium temperature is properly computed for variety of cases.
BOOST_AUTO_TEST_CASE( testEquilibriumTemperature )
{
    double airSpeed, noseRadius, wallEmissivity, equilibriumWallTemperature,
            adiabaticWallTemperature;
    double airDensity = 1.0E-5;
    double airTemperature = 300.0;
    for( unsigned int i = 0; i <= 10; i++ )
    {
        airSpeed = 1.0E3 + static_cast< double >( i ) * 1.0E3;
        double machNumber = airSpeed / 300.0;

        for( unsigned int j = 0; j <= 10; j++ )
        {
            noseRadius = 0.001 + static_cast< double >( j ) * 0.1;
            for( unsigned int k = 0; k < 10; k++ )
            {
                wallEmissivity = 0.1 + static_cast< double >( k ) * 0.1;
                adiabaticWallTemperature
                        = computeAdiabaticWallTemperature( airTemperature , machNumber );

                boost::function< double( const double ) > heatTransferFunction = boost::bind(
                            &computeFayRiddellHeatFlux, airDensity, airSpeed, airTemperature, noseRadius, _1 );

                equilibriumWallTemperature =
                        computeEquilibiumWallTemperature( heatTransferFunction, wallEmissivity, adiabaticWallTemperature );

                BOOST_CHECK_CLOSE_FRACTION(
                            wallEmissivity * electro_magnetism::computeBlackbodyRadiationIntensity( equilibriumWallTemperature ),
                            heatTransferFunction( equilibriumWallTemperature ), 1.0E-12 );

            }
        }
    }

}

//! Test if Fay-Ridell heat flux functions produce consistent results for variety of cases.
BOOST_AUTO_TEST_CASE( testFayRiddellHeatFluxConsistency )
{
    double airSpeed, noseRadius, equilibriumWallTemperature,
            adiabaticWallTemperature;
    double airDensity = 1.0E-5;
    double airTemperature = 300.0;
    double wallEmissivity = 0.7;
    double heatFlux1, heatFlux2, heatFlux3;
    for( unsigned int i = 0; i <= 10; i++ )
    {
        airSpeed = 1.0E3 + static_cast< double >( i ) * 1.0E3;
        double machNumber = airSpeed / 300.0;

        for( unsigned int j = 0; j <= 10; j++ )
        {
            noseRadius = 0.001 + static_cast< double >( j ) * 0.1;

            adiabaticWallTemperature
                    = computeAdiabaticWallTemperature( airTemperature , machNumber );

            boost::function< double( const double ) > heatTransferFunction = boost::bind(
                        &computeFayRiddellHeatFlux, airDensity, airSpeed, airTemperature, noseRadius, _1 );

            heatFlux1 = computeEquilibriumFayRiddellHeatFlux(
                        airDensity, airSpeed, airTemperature, machNumber, noseRadius, wallEmissivity );
            heatFlux2 = computeEquilibriumHeatflux(
                        heatTransferFunction, wallEmissivity, adiabaticWallTemperature );

            equilibriumWallTemperature =
                    computeEquilibiumWallTemperature( heatTransferFunction, wallEmissivity, adiabaticWallTemperature );
            heatFlux3 = computeFayRiddellHeatFlux(
                        airDensity, airSpeed, airTemperature, noseRadius, equilibriumWallTemperature );

            BOOST_CHECK_CLOSE_FRACTION(
                        heatFlux1, heatFlux2, 4.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION(
                        heatFlux1, heatFlux3, 4.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION(
                        heatFlux2, heatFlux3, 4.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }
}

//! Test if Fay-Ridell heat flux functions produce correct results for limit cases
BOOST_AUTO_TEST_CASE( testFayRiddellHeatFluxFunctions )
{
    double airSpeed = 6.0E3;
    double noseRadius = 0.1;
    double airDensity = 1.0E-5;
    double airTemperature = 300.0;

    double computedHeatFlux  = computeFayRiddellHeatFlux(
                airDensity, airSpeed, airTemperature, noseRadius, airTemperature );

    double expectedHeatFlux = 0.5 * FAY_RIDDEL_HEAT_FLUX_CONSTANT * std::pow(
                airSpeed, 3.0 ) * std::pow( airDensity, 0.5 ) * std::pow( noseRadius, -0.5 );
    BOOST_CHECK_CLOSE_FRACTION(
                computedHeatFlux, expectedHeatFlux, 4.0 * std::numeric_limits< double >::epsilon( ) );

    computedHeatFlux  = computeFayRiddellHeatFlux(
                    0.0, airSpeed, airTemperature, noseRadius, airTemperature );
    BOOST_CHECK_SMALL( std::fabs( computedHeatFlux ), std::numeric_limits< double >::epsilon( ) );

    computedHeatFlux  = computeFayRiddellHeatFlux(
                    airDensity, 0.0, airTemperature, noseRadius, airTemperature );
    BOOST_CHECK_SMALL( std::fabs( computedHeatFlux ), std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat

