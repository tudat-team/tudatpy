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
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Gravitation/tabulatedGravityFieldVariations.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_gravity_field_variations )

using namespace tudat::gravitation;
using namespace tudat::spice_interface;

//! Function to get nominal gravity field coefficients for Jupiter
/*!
 * Function to get nominal gravity field coefficients for Jupiter (returned by reference)
 * \param cosineCoefficients Jupiter's cosine coefficients
 * \param sineCoefficients Jupiter's sine coefficients
 */
void getNominalJupiterGravityField(
        Eigen::MatrixXd& cosineCoefficients, Eigen::MatrixXd& sineCoefficients )
{
    // Define unnormalized coefficients.
    double jupiterJ2 = 14.736E-3;
    double jupiterJ3 = 1.4E-6;
    double jupiterJ4 = -587.0E-6;
    double jupiterJ6 = 31.0E-6;
    double jupiterc22 = -0.03E-6;
    double jupiters22 = -0.007E-6;

    cosineCoefficients.setZero( 6, 6 );
    sineCoefficients.setZero( 6, 6 );

    // Normalize coefficients and set in matrices.
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) =  -jupiterJ2 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    cosineCoefficients( 2, 2 ) =  jupiterc22 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
    cosineCoefficients( 3, 0 ) =  -jupiterJ3 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 0 );
    cosineCoefficients( 4, 0 ) =  -jupiterJ4 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 0 );
    cosineCoefficients( 5, 0 ) =  -jupiterJ6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 6, 0 );
    sineCoefficients( 2, 2 ) =  jupiters22/ basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
}

//! Function to get tabulated gravity field variations
/*!
 * Function to get tabulated gravity field variations, based on arbitrary values, purely for test
 * purposes.
 * \param cosineCoefficientCorrections Map of corrections to cosine coefficients as function of
 * time.
 * \param sineCoefficientCorrections Map of corrections to sine coefficients as function of
 * time.
 */
void getTabulatedGravityFieldVariationValues(
        std::map< double, Eigen::MatrixXd >& cosineCoefficientCorrections,
        std::map< double, Eigen::MatrixXd >& sineCoefficientCorrections )
{

    // Define settings for coefficient variations
    double angularFrequency = 2.0 * mathematical_constants::PI / 86400.0;
    double amplitude = 1.0E-7;
    double startTime = 0.99E7;
    double endTime = 1.01E7;
    double timeStep = 3600.0;

    // Iterate over time steps
    double currentTime = startTime;
    while( currentTime < endTime )
    {
        // Initialize correctiuons to zero
        cosineCoefficientCorrections[ currentTime ].setZero( 4, 5 );
        sineCoefficientCorrections[ currentTime ].setZero( 4, 5 );

        // Set corrections for current time step.
        for( unsigned int i = 1; i < 5; i++ )
        {
            for( unsigned int j = 0; j <= i; j++ )
            {
                cosineCoefficientCorrections[ currentTime ]( i - 1, j ) =
                        amplitude *
                        std::cos( currentTime * angularFrequency + static_cast< double >( i * j ) );
                if( j != 0 )
                {
                    sineCoefficientCorrections[ currentTime ]( i - 1, j ) =
                            amplitude *
                            std::cos( currentTime * angularFrequency -
                                      static_cast< double >( i * j ) );
                }
            }
        }
        currentTime += timeStep;
    }
}

//! Function to get predefined tabulated gravity field variations.
/*!
 * Function to get predefined tabulated gravity field variations (non-physical values, purelu
 * for test purposes) using getTabulatedGravityFieldVariationValues function
 * \return Predefined tabulated gravity field variations.
 */
std::shared_ptr< TabulatedGravityFieldVariations > getTabulatedGravityFieldVariations( )
{
    // Get coefficient tables.
    std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections;
    std::map< double, Eigen::MatrixXd > sineCoefficientCorrections;
    getTabulatedGravityFieldVariationValues(
                cosineCoefficientCorrections, sineCoefficientCorrections );

    // Create correction object.
    return std::make_shared< TabulatedGravityFieldVariations >(
                cosineCoefficientCorrections, sineCoefficientCorrections, 1, 0 );
}

//! Function to get predefined gravity field variations object.
/*!
 *  Function to get predefined gravity field variations object, consisting of dummy tabulated
 *  variations from getTabulatedGravityFieldVariationValues and degree 2 tidal variations.
 * \return Predefined gravity field variations object.
 */
std::shared_ptr< GravityFieldVariationsSet > getTestGravityFieldVariations( )
{
    // Define bodies raising rides.
    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Io" );
    deformingBodies.push_back( "Europa" );

    // Retrieve required data of bodies raising tides.
    std::vector< std::function< Eigen::Vector6d( const double ) > >
            deformingBodyStateFunctions;
    std::vector< std::function< double( ) > > deformingBodyMasses;
    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
    {
        deformingBodyStateFunctions.push_back(
                    std::bind( &getBodyCartesianStateAtEpoch,
                                 deformingBodies.at( i ), "SSB", "J2000",
                                 "None", std::placeholders::_1 ) );
        deformingBodyMasses.push_back(
                    std::bind( &getBodyGravitationalParameter,
                                 deformingBodies.at( i ) ) );
    }

    // Define Love numbers ( constant for degree 2 only)
    std::complex< double > constantLoveNumber( 0.5, 0.5E-3 );
    std::vector< std::complex< double > > constantSingleDegreeLoveNumber( 3, constantLoveNumber );
    std::vector< std::vector< std::complex< double > > >
        loveNumbers( 1, constantSingleDegreeLoveNumber );

    // Set up gravity field variation of Jupiter due to Galilean moons.
    std::shared_ptr< GravityFieldVariations > solidBodyGravityFieldVariations =
            std::make_shared< BasicSolidBodyTideGravityFieldVariations >(
                std::bind( &getBodyCartesianStateAtEpoch,
                             "Jupiter", "SSB", "J2000", "None", std::placeholders::_1 ),
                std::bind( &computeRotationQuaternionBetweenFrames,
                             "J2000", "IAU_Jupiter", std::placeholders::_1 ),
                deformingBodyStateFunctions,
                getAverageRadius( "Jupiter" ),
                std::bind( &getBodyGravitationalParameter, "Jupiter" ),
                deformingBodyMasses, loveNumbers, deformingBodies );

    // Get tabulated gravity field variations.
    std::shared_ptr< GravityFieldVariations > tabulatedGravityFieldVariations =
            getTabulatedGravityFieldVariations( );

    // Create and return full gravity field variations object.
    return std::make_shared< GravityFieldVariationsSet >(
                std::vector< std::shared_ptr< GravityFieldVariations > >{ solidBodyGravityFieldVariations, tabulatedGravityFieldVariations },
                std::vector< BodyDeformationTypes >{ basic_solid_body, tabulated_variation },
                std::vector< std::string >{ "BasicTidal", "Tabulated" } );
}

BOOST_AUTO_TEST_CASE( testGravityFieldVariations )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_jup310_small.bsp" } );

    // Define properties of nominal field
    double gravitationalParameter = getBodyGravitationalParameter( "Jupiter" );
    double referenceRadius = getAverageRadius( "Jupiter" );
    Eigen::MatrixXd nominalCosineCoefficients;
    Eigen::MatrixXd nominalSineCoefficients;
    getNominalJupiterGravityField( nominalCosineCoefficients, nominalSineCoefficients );

    // Create gravity field corrections object and retrieve corrections
    std::vector< std::shared_ptr< GravityFieldVariations > > gravityFieldVariationsList =
            getTestGravityFieldVariations( )->getVariationObjects( );

    // Define data structures for storing expected gravity field variations.
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directGravityFieldVariations;
    Eigen::Matrix< double, 5, 5 > expectedCosineCoefficientsCorrections =
            Eigen::Matrix< double, 5, 5 >::Zero( );
    Eigen::Matrix< double, 5, 5 > expectedSineCoefficientsCorrections =
            Eigen::Matrix< double, 5, 5 >::Zero( );

    // Calculate gravity field variations directly.
    int minimumDegree, minimumOrder, numberOfDegrees, numberOfOrders;
    double testTime = 1.0E7;
    for( unsigned int i = 0; i < gravityFieldVariationsList.size( ); i++ )
    {
        minimumDegree = gravityFieldVariationsList.at( i )->getMinimumDegree( );
        numberOfDegrees = gravityFieldVariationsList.at( i )->getNumberOfDegrees( );

        minimumOrder = gravityFieldVariationsList.at( i )->getMinimumOrder( );
        numberOfOrders = gravityFieldVariationsList.at( i )->getNumberOfOrders( );

        directGravityFieldVariations =
                gravityFieldVariationsList.at( i )->
                calculateSphericalHarmonicsCorrections( testTime );
        expectedCosineCoefficientsCorrections.block(
                    minimumDegree, minimumOrder, numberOfDegrees, numberOfOrders ) +=
                directGravityFieldVariations.first;
        expectedSineCoefficientsCorrections.block(
                    minimumDegree, minimumOrder, numberOfDegrees, numberOfOrders ) +=
                directGravityFieldVariations.second;
    }

    // Create time-varying gravity field.
    std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
            std::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                gravitationalParameter, referenceRadius, nominalCosineCoefficients,
                nominalSineCoefficients, getTestGravityFieldVariations( ) );
    timeDependentGravityField->update( 2.0 * testTime );

    // Calculate variations for current test time.
    timeDependentGravityField->update( testTime );

    // Get corrections at current time step from timeDependentGravityField
    Eigen::MatrixXd perturbedCosineCoefficients =
            timeDependentGravityField->getCosineCoefficients( );
    Eigen::MatrixXd perturbedSineCoefficients =
            timeDependentGravityField->getSineCoefficients( );
    Eigen::MatrixXd calculatedCosineCoefficientCorrections =
            ( perturbedCosineCoefficients - nominalCosineCoefficients ).block( 0, 0, 5, 5 );
    Eigen::MatrixXd calculatedSineCoefficientCorrections =
            ( perturbedSineCoefficients - nominalSineCoefficients ).block( 0, 0, 5, 5 );

    // Compare corrections against expected variations.
    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            BOOST_CHECK_SMALL( calculatedCosineCoefficientCorrections( i, j ) -
                               expectedCosineCoefficientsCorrections( i, j ), 1.0E-18 );
            BOOST_CHECK_SMALL( calculatedSineCoefficientCorrections( i, j ) -
                               expectedSineCoefficientsCorrections( i, j ), 1.0E-18 );
        }
    }

    // Test calculated tidal corrections against manual corrections directly from Cartesian states
    // of perturbing bodies.
    {
        // Define Love numbers
        // Define Love numbers ( constant for degree 2 only)
        std::complex< double > constantLoveNumber( 0.5, 0.5E-3 );
        std::vector< std::complex< double > > constantSingleDegreeLoveNumber( 3, constantLoveNumber );
        std::vector< std::vector< std::complex< double > > >
            loveNumbers( 1, constantSingleDegreeLoveNumber );

        // Manually calculate tidal corrections directly.
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directIoTide =
                calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                    loveNumbers, getBodyGravitationalParameter( "Io" )
                    / getBodyGravitationalParameter( "Jupiter" ),
                    getAverageRadius( "Jupiter" ),  spice_interface::getBodyCartesianPositionAtEpoch(
                        "Io", "Jupiter", "IAU_Jupiter", "None", testTime ), 2, 2 );
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directEuropaTide =
                calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                    loveNumbers, getBodyGravitationalParameter( "Europa" )
                    / getBodyGravitationalParameter( "Jupiter" ),
                    getAverageRadius( "Jupiter" ),  spice_interface::getBodyCartesianPositionAtEpoch(
                        "Europa", "Jupiter", "IAU_Jupiter", "None", testTime ), 2, 2 );

        // Calculate tidal corrections from objects
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > tidalCorrectionsFromObject =
                ( ( timeDependentGravityField->getGravityFieldVariationsSet( )
                    ->getGravityFieldVariation( basic_solid_body ) ).second
                  ->calculateSphericalHarmonicsCorrections( testTime ) );

        // Compare results.
        Eigen::MatrixXd directCosineCorrections = directIoTide.first + directEuropaTide.first;
        Eigen::MatrixXd directSineCorrections = directIoTide.second + directEuropaTide.second;

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( directCosineCorrections( 2, i )
                               - tidalCorrectionsFromObject.first( 0, i ), 1.0E-20 );
            BOOST_CHECK_SMALL( directSineCorrections( 2, i )
                               - tidalCorrectionsFromObject.second( 0, i ), 1.0E-20 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

