/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <sstream>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/tabulatedGravityFieldVariations.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_gravity_field_variations )

using namespace tudat::gravitation;
using namespace tudat::spice_interface;

// Compare coefficient rates with a central difference of the model's own variation values.
void checkGravityFieldVariationDerivative( const std::shared_ptr< GravityFieldVariations >& variation,
                                           const double time,
                                           const double timeStep,
                                           const double relativeTolerance = 1.0e-8 )
{
    const auto rates = variation->calculateSphericalHarmonicsCorrectionsTimeDerivative( time );
    const auto before = variation->calculateSphericalHarmonicsCorrections( time - timeStep );
    const auto after = variation->calculateSphericalHarmonicsCorrections( time + timeStep );
    const Eigen::MatrixXd numericalCosineRate = ( after.first - before.first ) / ( 2.0 * timeStep );
    const Eigen::MatrixXd numericalSineRate = ( after.second - before.second ) / ( 2.0 * timeStep );

    // Both coefficient blocks must agree with the numerical rates, including zero sine entries.
    BOOST_CHECK_SMALL( ( rates.first - numericalCosineRate ).norm( ), relativeTolerance * numericalCosineRate.norm( ) + 1.0e-25 );
    BOOST_CHECK_SMALL( ( rates.second - numericalSineRate ).norm( ), relativeTolerance * numericalSineRate.norm( ) + 1.0e-25 );
}

//! Function to get nominal gravity field coefficients for Jupiter
/*!
 * Function to get nominal gravity field coefficients for Jupiter (returned by reference)
 * \param cosineCoefficients Jupiter's cosine coefficients
 * \param sineCoefficients Jupiter's sine coefficients
 */
void getNominalJupiterGravityField( Eigen::MatrixXd& cosineCoefficients, Eigen::MatrixXd& sineCoefficients )
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
    cosineCoefficients( 2, 0 ) = -jupiterJ2 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    cosineCoefficients( 2, 2 ) = jupiterc22 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
    cosineCoefficients( 3, 0 ) = -jupiterJ3 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 0 );
    cosineCoefficients( 4, 0 ) = -jupiterJ4 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 0 );
    cosineCoefficients( 5, 0 ) = -jupiterJ6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 6, 0 );
    sineCoefficients( 2, 2 ) = jupiters22 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
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
void getTabulatedGravityFieldVariationValues( std::map< double, Eigen::MatrixXd >& cosineCoefficientCorrections,
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
                        amplitude * std::cos( currentTime * angularFrequency + static_cast< double >( i * j ) );
                if( j != 0 )
                {
                    sineCoefficientCorrections[ currentTime ]( i - 1, j ) =
                            amplitude * std::cos( currentTime * angularFrequency - static_cast< double >( i * j ) );
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
    getTabulatedGravityFieldVariationValues( cosineCoefficientCorrections, sineCoefficientCorrections );

    // Create correction object.
    return std::make_shared< TabulatedGravityFieldVariations >( cosineCoefficientCorrections, sineCoefficientCorrections, 1, 0 );
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
    std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyStateFunctions;
    std::vector< std::function< double( ) > > deformingBodyMasses;
    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
    {
        deformingBodyStateFunctions.push_back(
                std::bind( &getBodyCartesianStateAtEpoch, deformingBodies.at( i ), "SSB", "J2000", "None", std::placeholders::_1 ) );
        deformingBodyMasses.push_back( std::bind( &getBodyGravitationalParameter, deformingBodies.at( i ) ) );
    }

    // Define Love numbers ( constant for degree 2 only)
    std::complex< double > constantLoveNumber( 0.5, 0.5E-3 );
    std::vector< std::complex< double > > constantSingleDegreeLoveNumber( 3, constantLoveNumber );
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    loveNumbers[ 2 ] = constantSingleDegreeLoveNumber;

    // Set up gravity field variation of Jupiter due to Galilean moons.
    std::shared_ptr< GravityFieldVariations > solidBodyGravityFieldVariations =
            std::make_shared< BasicSolidBodyTideGravityFieldVariations >(
                    std::bind( &getBodyCartesianStateAtEpoch, "Jupiter", "SSB", "J2000", "None", std::placeholders::_1 ),
                    std::bind( &computeRotationQuaternionBetweenFrames, "J2000", "IAU_Jupiter", std::placeholders::_1 ),
                    deformingBodyStateFunctions,
                    getAverageRadius( "Jupiter" ),
                    std::bind( &getBodyGravitationalParameter, "Jupiter" ),
                    deformingBodyMasses,
                    loveNumbers,
                    deformingBodies );

    // Get tabulated gravity field variations.
    std::shared_ptr< GravityFieldVariations > tabulatedGravityFieldVariations = getTabulatedGravityFieldVariations( );

    // Create and return full gravity field variations object.
    return std::make_shared< GravityFieldVariationsSet >(
            std::vector< std::shared_ptr< GravityFieldVariations > >{ solidBodyGravityFieldVariations, tabulatedGravityFieldVariations },
            std::vector< BodyDeformationTypes >{ basic_solid_body, tabulated_variation },
            std::vector< std::string >{ "BasicTidal", "Tabulated" } );
}

// Compare combined tidal/tabulated values with direct corrections; check tabulated rates against
// finite differences in several intervals and after resetting the table, plus the nonlinear fallback.
BOOST_AUTO_TEST_CASE( testGravityFieldVariations )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

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
    Eigen::Matrix< double, 5, 5 > expectedCosineCoefficientsCorrections = Eigen::Matrix< double, 5, 5 >::Zero( );
    Eigen::Matrix< double, 5, 5 > expectedSineCoefficientsCorrections = Eigen::Matrix< double, 5, 5 >::Zero( );

    // Calculate gravity field variations directly.
    int minimumDegree, minimumOrder, numberOfDegrees, numberOfOrders;
    double testTime = 1.0E7;
    for( unsigned int i = 0; i < gravityFieldVariationsList.size( ); i++ )
    {
        minimumDegree = gravityFieldVariationsList.at( i )->getMinimumDegree( );
        numberOfDegrees = gravityFieldVariationsList.at( i )->getNumberOfDegrees( );

        minimumOrder = gravityFieldVariationsList.at( i )->getMinimumOrder( );
        numberOfOrders = gravityFieldVariationsList.at( i )->getNumberOfOrders( );

        directGravityFieldVariations = gravityFieldVariationsList.at( i )->calculateSphericalHarmonicsCorrections( testTime );
        expectedCosineCoefficientsCorrections.block( minimumDegree, minimumOrder, numberOfDegrees, numberOfOrders ) +=
                directGravityFieldVariations.first;
        expectedSineCoefficientsCorrections.block( minimumDegree, minimumOrder, numberOfDegrees, numberOfOrders ) +=
                directGravityFieldVariations.second;
    }

    // Create time-varying gravity field.
    std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
            std::make_shared< TimeDependentSphericalHarmonicsGravityField >( gravitationalParameter,
                                                                             referenceRadius,
                                                                             nominalCosineCoefficients,
                                                                             nominalSineCoefficients,
                                                                             getTestGravityFieldVariations( ) );
    // Reject indices outside the nominal coefficient blocks before attempting a matrix write.
    BOOST_CHECK_THROW( timeDependentGravityField->setNominalCosineCoefficient( nominalCosineCoefficients.rows( ), 0, 0.0 ),
                       std::runtime_error );
    BOOST_CHECK_THROW( timeDependentGravityField->setNominalSineCoefficient( 0, nominalSineCoefficients.cols( ), 0.0 ),
                       std::runtime_error );
    BOOST_CHECK_THROW( timeDependentGravityField->setNominalCosineCoefficient( -1, 0, 0.0 ), std::runtime_error );
    BOOST_CHECK_THROW( timeDependentGravityField->setNominalSineCoefficient( 0, -1, 0.0 ), std::runtime_error );
    timeDependentGravityField->update( 2.0 * testTime );

    // Calculate variations for current test time.
    timeDependentGravityField->update( testTime );

    // Get corrections at current time step from timeDependentGravityField
    Eigen::MatrixXd perturbedCosineCoefficients = timeDependentGravityField->getCosineCoefficients( );
    Eigen::MatrixXd perturbedSineCoefficients = timeDependentGravityField->getSineCoefficients( );
    Eigen::MatrixXd calculatedCosineCoefficientCorrections =
            ( perturbedCosineCoefficients - nominalCosineCoefficients ).block( 0, 0, 5, 5 );
    Eigen::MatrixXd calculatedSineCoefficientCorrections = ( perturbedSineCoefficients - nominalSineCoefficients ).block( 0, 0, 5, 5 );

    // Compare corrections against expected variations.
    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            BOOST_CHECK_SMALL( calculatedCosineCoefficientCorrections( i, j ) - expectedCosineCoefficientsCorrections( i, j ), 1.0E-18 );
            BOOST_CHECK_SMALL( calculatedSineCoefficientCorrections( i, j ) - expectedSineCoefficientsCorrections( i, j ), 1.0E-18 );
        }
    }

    // Test calculated tidal corrections against manual corrections directly from Cartesian states
    // of perturbing bodies.
    {
        // Define Love numbers
        // Define Love numbers ( constant for degree 2 only)
        std::complex< double > constantLoveNumber( 0.5, 0.5E-3 );
        std::vector< std::complex< double > > constantSingleDegreeLoveNumber( 3, constantLoveNumber );
        std::map< int, std::vector< std::complex< double > > > loveNumbers;
        loveNumbers[ 2 ] = constantSingleDegreeLoveNumber;

        // Manually calculate tidal corrections directly.
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directIoTide = calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                loveNumbers,
                getBodyGravitationalParameter( "Io" ) / getBodyGravitationalParameter( "Jupiter" ),
                getAverageRadius( "Jupiter" ),
                spice_interface::getBodyCartesianPositionAtEpoch( "Io", "Jupiter", "IAU_Jupiter", "None", testTime ),
                2,
                2 );
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directEuropaTide = calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                loveNumbers,
                getBodyGravitationalParameter( "Europa" ) / getBodyGravitationalParameter( "Jupiter" ),
                getAverageRadius( "Jupiter" ),
                spice_interface::getBodyCartesianPositionAtEpoch( "Europa", "Jupiter", "IAU_Jupiter", "None", testTime ),
                2,
                2 );

        // Calculate tidal corrections from objects
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > tidalCorrectionsFromObject =
                ( ( timeDependentGravityField->getGravityFieldVariationsSet( )->getGravityFieldVariation( basic_solid_body ) )
                          .second->calculateSphericalHarmonicsCorrections( testTime ) );

        // Compare results.
        Eigen::MatrixXd directCosineCorrections = directIoTide.first + directEuropaTide.first;
        Eigen::MatrixXd directSineCorrections = directIoTide.second + directEuropaTide.second;

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( directCosineCorrections( 2, i ) - tidalCorrectionsFromObject.first( 0, i ), 1.0E-19 );
            BOOST_CHECK_SMALL( directSineCorrections( 2, i ) - tidalCorrectionsFromObject.second( 0, i ), 1.0E-19 );
        }
    }

    const auto tabulatedVariation = getTabulatedGravityFieldVariations( );
    auto cosineTable = tabulatedVariation->getCosineCoefficientCorrections( );
    auto sineTable = tabulatedVariation->getSineCoefficientCorrections( );
    const double firstEpoch = cosineTable.begin( )->first;
    const double tableStep = 3600.0;
    // Interior points in three different intervals check the slope without crossing a knot.
    for( const double intervalOffset : { 0.25, 1.5, 2.75 } )
    {
        checkGravityFieldVariationDerivative( tabulatedVariation, firstEpoch + intervalOffset * tableStep, 10.0, 1.0e-10 );
    }

    for( auto& entry : cosineTable )
    {
        entry.second *= 2.0;
    }
    for( auto& entry : sineTable )
    {
        entry.second *= -3.0;
    }
    tabulatedVariation->resetCoefficientInterpolator( cosineTable, sineTable );
    // Replacing the table must also replace the cached slopes for both coefficient blocks.
    checkGravityFieldVariationDerivative( tabulatedVariation, firstEpoch + 1.5 * tableStep, 10.0, 1.0e-10 );

    TabulatedGravityFieldVariations nonlinearVariation(
            cosineTable,
            sineTable,
            1,
            0,
            std::make_shared< interpolators::InterpolatorSettings >( interpolators::cubic_spline_interpolator ) );
    std::ostringstream warnings;
    std::streambuf* originalBuffer = std::cerr.rdbuf( warnings.rdbuf( ) );
    const auto firstRates = nonlinearVariation.calculateSphericalHarmonicsCorrectionsTimeDerivative( firstEpoch );
    const std::string firstWarning = warnings.str( );
    const auto secondRates = nonlinearVariation.calculateSphericalHarmonicsCorrectionsTimeDerivative( firstEpoch + tableStep );
    std::cerr.rdbuf( originalBuffer );

    // Unsupported interpolation returns zero rates and emits a warning only on the first call.
    BOOST_CHECK( !firstWarning.empty( ) );
    BOOST_CHECK_EQUAL( warnings.str( ), firstWarning );
    BOOST_CHECK_SMALL( firstRates.first.norm( ) + firstRates.second.norm( ) + secondRates.first.norm( ) + secondRates.second.norm( ),
                       1.0e-25 );
}

// Keep the reference SPICE positions, with exactly consistent linear position/velocity functions.
// The kernels' separately supplied velocities are not exact derivatives of their interpolated positions.
std::function< Eigen::Vector6d( double ) > getGravityVariationTestStateFunction( const std::string& body )
{
    const double referenceTime = 1.0e7;
    const Eigen::Vector6d referenceState = getBodyCartesianStateAtEpoch( body, "SSB", "J2000", "None", referenceTime );
    return [ referenceTime, referenceState ]( const double time ) {
        Eigen::Vector6d state = referenceState;
        state.head< 3 >( ) += ( time - referenceTime ) * referenceState.tail< 3 >( );
        return state;
    };
}

std::shared_ptr< BasicSolidBodyTideGravityFieldVariations > getBasicGravityFieldVariation( const bool includeHigherDegrees = false )
{
    // Define bodies raising rides.
    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Io" );
    deformingBodies.push_back( "Europa" );

    // Retrieve required data of bodies raising tides.
    std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyStateFunctions;
    std::vector< std::function< double( ) > > deformingBodyMasses;
    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
    {
        deformingBodyStateFunctions.push_back( getGravityVariationTestStateFunction( deformingBodies.at( i ) ) );
        deformingBodyMasses.push_back( std::bind( &getBodyGravitationalParameter, deformingBodies.at( i ) ) );
    }

    // Define Love numbers ( constant for degree 2 only)
    std::vector< std::complex< double > > degreeTwoLoveNumber = { std::complex< double >( 0.1, 0.0 ),
                                                                  std::complex< double >( 0.2, 0.0 ),
                                                                  std::complex< double >( 0.3, 0.0 ) };
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    loveNumbers[ 2 ] = degreeTwoLoveNumber;

    std::map< int, std::vector< double > > meanCosineForcing, meanSineForcing;
    if( includeHigherDegrees )
    {
        loveNumbers[ 3 ] = { 0.1, 0.2, 0.3, 0.4 };
        loveNumbers[ 4 ] = { 0.1, std::complex< double >( 0.2, 0.03 ), 0.3, 0.4 };
        meanCosineForcing = generateZeroMeanTermsFromReference( loveNumbers );
        meanSineForcing = generateZeroMeanTermsFromReference( loveNumbers );
        meanCosineForcing[ 4 ][ 1 ] = 1.0e-12;
        meanSineForcing[ 4 ][ 1 ] = -2.0e-12;
    }

    // Set up gravity field variation of Jupiter due to Galilean moons.
    const auto variation = std::make_shared< BasicSolidBodyTideGravityFieldVariations >(
            getGravityVariationTestStateFunction( "Jupiter" ),
            std::bind( &computeRotationQuaternionBetweenFrames, "J2000", "IAU_Jupiter", std::placeholders::_1 ),
            deformingBodyStateFunctions,
            getAverageRadius( "Jupiter" ),
            std::bind( &getBodyGravitationalParameter, "Jupiter" ),
            deformingBodyMasses,
            loveNumbers,
            deformingBodies,
            meanCosineForcing,
            meanSineForcing );
    variation->resetRotationDerivativeFunction(
            std::bind( &computeRotationMatrixDerivativeBetweenFrames, "J2000", "IAU_Jupiter", std::placeholders::_1 ) );
    return variation;
}

std::shared_ptr< ModeCoupledSolidBodyTideGravityFieldVariations > getModeCoupledGravityFieldVariation( const int index )
{
    // Define bodies raising rides.
    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Io" );
    deformingBodies.push_back( "Europa" );

    // Retrieve required data of bodies raising tides.
    std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyStateFunctions;
    std::vector< std::function< double( ) > > deformingBodyMasses;
    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
    {
        deformingBodyStateFunctions.push_back( getGravityVariationTestStateFunction( deformingBodies.at( i ) ) );
        deformingBodyMasses.push_back( std::bind( &getBodyGravitationalParameter, deformingBodies.at( i ) ) );
    }

    // Define Love numbers ( constant for degree 2 only)
    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > loveNumbers;
    if( index == 0 )
    {
        loveNumbers[ std::make_pair( 2, 0 ) ][ std::make_pair( 2, 0 ) ] = 0.1;
        loveNumbers[ std::make_pair( 2, 1 ) ][ std::make_pair( 2, 1 ) ] = 0.2;
        loveNumbers[ std::make_pair( 2, 2 ) ][ std::make_pair( 2, 2 ) ] = 0.3;
    }
    else if( index == 1 )
    {
        loveNumbers[ std::make_pair( 2, 0 ) ][ std::make_pair( 2, 2 ) ] = 0.1;
        loveNumbers[ std::make_pair( 2, 2 ) ][ std::make_pair( 2, 1 ) ] = 0.3;
        loveNumbers[ std::make_pair( 2, 1 ) ][ std::make_pair( 2, 0 ) ] = 0.2;
    }
    else if( index == 2 )
    {
        loveNumbers[ std::make_pair( 2, 0 ) ][ std::make_pair( 3, 3 ) ] = 0.1;
        loveNumbers[ std::make_pair( 2, 2 ) ][ std::make_pair( 3, 2 ) ] = 0.3;
        loveNumbers[ std::make_pair( 2, 1 ) ][ std::make_pair( 3, 1 ) ] = 0.2;
    }

    // Set up gravity field variation of Jupiter due to Galilean moons.
    const auto variation = std::make_shared< ModeCoupledSolidBodyTideGravityFieldVariations >(
            getGravityVariationTestStateFunction( "Jupiter" ),
            std::bind( &computeRotationQuaternionBetweenFrames, "J2000", "IAU_Jupiter", std::placeholders::_1 ),
            deformingBodyStateFunctions,
            getAverageRadius( "Jupiter" ),
            std::bind( &getBodyGravitationalParameter, "Jupiter" ),
            deformingBodyMasses,
            loveNumbers,
            deformingBodies );
    variation->resetRotationDerivativeFunction(
            std::bind( &computeRotationMatrixDerivativeBetweenFrames, "J2000", "IAU_Jupiter", std::placeholders::_1 ) );
    return variation;
}

// Compare the existing basic/coupled tide values, and check their analytical rates against
// central differences at two epochs, including body rotation and forcing/response index changes.
BOOST_AUTO_TEST_CASE( testModeCoupledGravityFieldVariations )
{
    std::cout << std::endl << std::endl << std::endl;
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::shared_ptr< ModeCoupledSolidBodyTideGravityFieldVariations > modeCoupledControlVariationModel =
            getModeCoupledGravityFieldVariation( 0 );
    std::shared_ptr< ModeCoupledSolidBodyTideGravityFieldVariations > modeCoupledOrderOffsetVariationModel =
            getModeCoupledGravityFieldVariation( 1 );
    std::shared_ptr< ModeCoupledSolidBodyTideGravityFieldVariations > modeCoupledDegreeOrderOffsetVariationModel =
            getModeCoupledGravityFieldVariation( 2 );
    std::shared_ptr< BasicSolidBodyTideGravityFieldVariations > basicVariationModel = getBasicGravityFieldVariation( );

    double testTime = 1.0E7;
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > modeCoupledControlVariations =
            modeCoupledControlVariationModel->calculateBasicSphericalHarmonicsCorrections( testTime );
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > basicVariations =
            basicVariationModel->calculateBasicSphericalHarmonicsCorrections( testTime );
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > modeCoupledOrderOffsetVariations =
            modeCoupledOrderOffsetVariationModel->calculateBasicSphericalHarmonicsCorrections( testTime );
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > modeCoupledDegreeOrderOffsetVariations =
            modeCoupledDegreeOrderOffsetVariationModel->calculateBasicSphericalHarmonicsCorrections( testTime );

    std::cout << modeCoupledDegreeOrderOffsetVariations.first << std::endl
              << std::endl
              << modeCoupledDegreeOrderOffsetVariations.second << std::endl
              << std::endl
              << std::endl;
    std::cout << modeCoupledOrderOffsetVariations.first << std::endl
              << std::endl
              << modeCoupledOrderOffsetVariations.second << std::endl
              << std::endl
              << std::endl;
    std::cout << modeCoupledControlVariations.first << std::endl
              << std::endl
              << modeCoupledControlVariations.second << std::endl
              << std::endl
              << std::endl;
    std::cout << basicVariations.first << std::endl << std::endl << basicVariations.second << std::endl << std::endl << std::endl;

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( modeCoupledControlVariations.first( 0, i ) - basicVariations.first( 0, i ), 1.0E-19 );
        BOOST_CHECK_SMALL( modeCoupledControlVariations.second( 0, i ) - basicVariations.second( 0, i ), 1.0E-19 );
    }

    BOOST_CHECK_SMALL( modeCoupledOrderOffsetVariations.first( 0, 0 ) - modeCoupledControlVariations.first( 0, 1 ), 1.0E-19 );
    BOOST_CHECK_SMALL( modeCoupledOrderOffsetVariations.first( 0, 1 ) - modeCoupledControlVariations.first( 0, 2 ), 1.0E-19 );
    BOOST_CHECK_SMALL( modeCoupledOrderOffsetVariations.first( 0, 2 ) - modeCoupledControlVariations.first( 0, 0 ), 1.0E-19 );
    BOOST_CHECK_SMALL( modeCoupledOrderOffsetVariations.second( 0, 1 ) - modeCoupledControlVariations.second( 0, 2 ), 1.0E-19 );

    BOOST_CHECK_SMALL( modeCoupledDegreeOrderOffsetVariations.first( 1, 1 ) - basicVariations.first( 0, 1 ), 1.0E-19 );
    BOOST_CHECK_SMALL( modeCoupledDegreeOrderOffsetVariations.first( 1, 2 ) - basicVariations.first( 0, 2 ), 1.0E-19 );
    BOOST_CHECK_SMALL( modeCoupledDegreeOrderOffsetVariations.first( 1, 3 ) - basicVariations.first( 0, 0 ), 1.0E-19 );
    BOOST_CHECK_SMALL( modeCoupledDegreeOrderOffsetVariations.second( 1, 1 ) - basicVariations.second( 0, 1 ), 1.0E-19 );
    BOOST_CHECK_SMALL( modeCoupledDegreeOrderOffsetVariations.second( 1, 2 ) - basicVariations.second( 0, 2 ), 1.0E-19 );

    // Each existing mapping must differentiate its own cosine/sine variations correctly.
    for( const double time : { testTime, testTime + 3600.0 } )
    {
        checkGravityFieldVariationDerivative( basicVariationModel, time, 1.0, 5.0e-8 );
        checkGravityFieldVariationDerivative( modeCoupledControlVariationModel, time, 1.0, 5.0e-8 );
        checkGravityFieldVariationDerivative( modeCoupledOrderOffsetVariationModel, time, 1.0, 5.0e-8 );
        checkGravityFieldVariationDerivative( modeCoupledDegreeOrderOffsetVariationModel, time, 1.0, 5.0e-8 );
    }

    // Changing a Love number must immediately change the rate, including its imaginary component.
    auto loveNumbers = basicVariationModel->getLoveNumbersOfDegree( 2 );
    loveNumbers[ 1 ] = std::complex< double >( 0.2, 0.03 );
    basicVariationModel->resetLoveNumbersOfDegree( loveNumbers, 2 );
    checkGravityFieldVariationDerivative( basicVariationModel, testTime, 1.0, 5.0e-8 );

    // Degrees three/four exercise the general recurrence, complex Love numbers, and constant mean offsets.
    // Agreement with the differentiated values also checks the corrected degree-four odd-order signs.
    const auto higherDegreeVariation = getBasicGravityFieldVariation( true );
    checkGravityFieldVariationDerivative( higherDegreeVariation, testTime, 1.0, 5.0e-8 );
    auto degreeFourLoveNumbers = higherDegreeVariation->getLoveNumbersOfDegree( 4 );
    degreeFourLoveNumbers.push_back( 0.5 );
    higherDegreeVariation->resetLoveNumbersOfDegree( degreeFourLoveNumbers, 4 );
    // Enabling another order must resize both the forcing buffers and the constant mean-forcing arrays.
    checkGravityFieldVariationDerivative( higherDegreeVariation, testTime + 3600.0, 1.0, 5.0e-8 );
}

void getPeriodicGravityFieldVariationSettings( std::vector< Eigen::MatrixXd >& cosineShAmplitudesCosineTime,
                                               std::vector< Eigen::MatrixXd >& cosineShAmplitudesSineTime,
                                               std::vector< Eigen::MatrixXd >& sineShAmplitudesCosineTime,
                                               std::vector< Eigen::MatrixXd >& sineShAmplitudesSineTime,
                                               std::vector< double >& frequencies,
                                               double& referenceEpoch,
                                               int& minimumDegree,
                                               int& minimumOrder,
                                               const int type = 0 )
{
    cosineShAmplitudesCosineTime.clear( );
    cosineShAmplitudesSineTime.clear( );
    sineShAmplitudesCosineTime.clear( );
    sineShAmplitudesSineTime.clear( );
    frequencies.clear( );

    cosineShAmplitudesCosineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 2.0E-8, -2.0E-8, 4.2E-9, 4.0E-8, -3.7E-9, 2.9E-7 ).finished( ) );
    cosineShAmplitudesCosineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 5.0E-8, -6.0E-8, 2.2E-9, 4.0E-7, -3.4E-7, 8.9E-7 ).finished( ) );
    cosineShAmplitudesCosineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 8.0E-8, 2.0E-7, 7.2E-8, 7.5E-7, -1.9E-7, 2.5E-7 ).finished( ) );

    cosineShAmplitudesSineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 3.0E-8, -6.0E-8, 2.2E-9, 9.0E-8, -6.7E-9, 3.9E-7 ).finished( ) );
    cosineShAmplitudesSineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 8.0E-8, -9.0E-8, 4.2E-9, 2.0E-7, -6.4E-7, 4.9E-7 ).finished( ) );
    cosineShAmplitudesSineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 3.0E-8, 6.0E-7, 8.2E-8, 1.5E-7, 5.9E-7, 9.5E-7 ).finished( ) );

    sineShAmplitudesCosineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 0.0, 6.0E-8, 5.2E-9, 0.0, 7.4E-8, 5.6E-7 ).finished( ) );
    sineShAmplitudesCosineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 0.0, 5.5E-8, 9.2E-9, 0.0, 4.4E-8, 1.1E-7 ).finished( ) );
    sineShAmplitudesCosineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 0.0, 7.8E-7, 4.2E-8, 0.0, 8.5E-7, -2.5E-8 ).finished( ) );

    sineShAmplitudesSineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 0.0, 5.0E-8, -4.2E-9, 0.0, 4.4E-8, -1.6E-7 ).finished( ) );
    sineShAmplitudesSineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 0.0, -4.5E-8, 6.2E-9, 0.0, 6.4E-8, 2.1E-7 ).finished( ) );
    sineShAmplitudesSineTime.push_back( ( Eigen::MatrixXd( 2, 3 ) << 0.0, 5.8E-7, 6.2E-8, 0.0, 2.5E-7, 5.5E-8 ).finished( ) );

    frequencies.resize( 3 );

    frequencies = { 2.0E-5, 7.0E-5, 4.0E-6 };
    referenceEpoch = 5.0E8;

    if( type == 0 )
    {
        minimumDegree = 2;
        minimumOrder = 0;
    }
    else if( type == 1 )
    {
        minimumDegree = 3;
        minimumOrder = 0;
    }
    else if( type == 2 )
    {
        minimumDegree = 3;
        minimumOrder = 1;
    }
}

// Check periodic values for the existing frequency/block settings, and compare their rates with
// central differences before, at, and after the reference epoch.
BOOST_AUTO_TEST_CASE( testPeriodicGravityFieldVariations )
{
    using namespace tudat::simulation_setup;

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define properties of nominal field
    double gravitationalParameter = getBodyGravitationalParameter( "Jupiter" );
    double referenceRadius = getAverageRadius( "Jupiter" );
    Eigen::MatrixXd nominalCosineCoefficients;
    Eigen::MatrixXd nominalSineCoefficients;
    getNominalJupiterGravityField( nominalCosineCoefficients, nominalSineCoefficients );
    {
        double testTime = 1.0E7;
        for( int k = 0; k < 4; k++ )
        {
            std::vector< Eigen::MatrixXd > cosineShAmplitudesCosineTime;
            std::vector< Eigen::MatrixXd > cosineShAmplitudesSineTime;
            std::vector< Eigen::MatrixXd > sineShAmplitudesCosineTime;
            std::vector< Eigen::MatrixXd > sineShAmplitudesSineTime;
            std::vector< double > frequencies;
            double referenceEpoch;
            int minimumDegree;
            int minimumOrder;
            int testIndex = ( k < 3 ) ? k : 0;
            getPeriodicGravityFieldVariationSettings( cosineShAmplitudesCosineTime,
                                                      cosineShAmplitudesSineTime,
                                                      sineShAmplitudesCosineTime,
                                                      sineShAmplitudesSineTime,
                                                      frequencies,
                                                      referenceEpoch,
                                                      minimumDegree,
                                                      minimumOrder,
                                                      testIndex );

            std::shared_ptr< PeriodicGravityFieldVariationsSettings > variationSettings;
            if( k < 3 )
            {
                variationSettings = std::dynamic_pointer_cast< PeriodicGravityFieldVariationsSettings >(
                        periodicGravityFieldVariationsSettings( cosineShAmplitudesCosineTime,
                                                                cosineShAmplitudesSineTime,
                                                                sineShAmplitudesCosineTime,
                                                                sineShAmplitudesSineTime,
                                                                frequencies,
                                                                referenceEpoch,
                                                                minimumDegree,
                                                                minimumOrder ) );
            }
            else
            {
                variationSettings = std::dynamic_pointer_cast< PeriodicGravityFieldVariationsSettings >(
                        periodicGravityFieldVariationsSettingsSingleFrequency( cosineShAmplitudesCosineTime.at( 0 ),
                                                                               cosineShAmplitudesSineTime.at( 0 ),
                                                                               sineShAmplitudesCosineTime.at( 0 ),
                                                                               sineShAmplitudesSineTime.at( 0 ),
                                                                               frequencies.at( 0 ),
                                                                               referenceEpoch,
                                                                               minimumDegree,
                                                                               minimumOrder ) );
            }

            // Create time-varying gravity field.
            std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                    std::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                            gravitationalParameter, referenceRadius, nominalCosineCoefficients, nominalSineCoefficients );

            SystemOfBodies dummySystem;
            dummySystem.createEmptyBody( "Jupiter" );
            dummySystem.getBody( "Jupiter" )->setGravityFieldModel( timeDependentGravityField );
            std::vector< std::shared_ptr< GravityFieldVariationSettings > > variationSettingsList;
            variationSettingsList.push_back( variationSettings );

            std::shared_ptr< gravitation::GravityFieldVariationsSet > variations =
                    createGravityFieldModelVariationsSet( "Jupiter", dummySystem, variationSettingsList );
            timeDependentGravityField->updateCorrectionFunctions( );

            // Exercise the same single-frequency and multiple-frequency models used for the value checks.
            for( const double derivativeTime : { referenceEpoch - 2.0e4, referenceEpoch, referenceEpoch + 3.0e4 } )
            {
                checkGravityFieldVariationDerivative( variations->getVariationObjects( ).at( 0 ), derivativeTime, 1.0 );
            }

            if( k == 0 )
            {
                const auto periodicVariation =
                        std::dynamic_pointer_cast< PeriodicGravityFieldVariations >( variations->getVariationObjects( ).at( 0 ) );
                auto invalidAmplitudes = cosineShAmplitudesCosineTime;
                invalidAmplitudes.pop_back( );
                // Reject inconsistent counts and block sizes when resetting, before they reach the evaluation loop.
                BOOST_CHECK_THROW( periodicVariation->resetCosineShAmplitudesCosineTime( invalidAmplitudes ), std::runtime_error );
                invalidAmplitudes = cosineShAmplitudesCosineTime;
                invalidAmplitudes.front( ) = Eigen::MatrixXd::Zero( 1, 1 );
                BOOST_CHECK_THROW( periodicVariation->resetCosineShAmplitudesCosineTime( invalidAmplitudes ), std::runtime_error );
                // A rejected reset must leave the original variations and their rates consistent.
                checkGravityFieldVariationDerivative( periodicVariation, referenceEpoch + 3.0e4, 1.0 );
            }

            if( k < 3 )
            {
                timeDependentGravityField->update( 2.0 * testTime );

                // Calculate variations for current test time.
                timeDependentGravityField->update( testTime );

                // Get corrections at current time step from timeDependentGravityField
                Eigen::MatrixXd perturbedCosineCoefficients = timeDependentGravityField->getCosineCoefficients( );
                Eigen::MatrixXd perturbedSineCoefficients = timeDependentGravityField->getSineCoefficients( );
                Eigen::MatrixXd calculatedCosineCoefficientCorrections =
                        ( perturbedCosineCoefficients - nominalCosineCoefficients ).block( 0, 0, 5, 5 );
                Eigen::MatrixXd calculatedSineCoefficientCorrections =
                        ( perturbedSineCoefficients - nominalSineCoefficients ).block( 0, 0, 5, 5 );

                Eigen::MatrixXd expectedCosineCoefficientsCorrections = Eigen::MatrixXd::Zero( 5, 5 );
                Eigen::MatrixXd expectedSineCoefficientsCorrections = Eigen::MatrixXd::Zero( 5, 5 );

                for( unsigned int j = 0; j < cosineShAmplitudesCosineTime.size( ); j++ )
                {
                    expectedCosineCoefficientsCorrections.block( minimumDegree, minimumOrder, 2, 3 ) +=
                            cosineShAmplitudesCosineTime.at( j ) * std::cos( frequencies.at( j ) * ( testTime - referenceEpoch ) ) +
                            cosineShAmplitudesSineTime.at( j ) * std::sin( frequencies.at( j ) * ( testTime - referenceEpoch ) );
                    expectedSineCoefficientsCorrections.block( minimumDegree, minimumOrder, 2, 3 ) +=
                            sineShAmplitudesCosineTime.at( j ) * std::cos( frequencies.at( j ) * ( testTime - referenceEpoch ) ) +
                            sineShAmplitudesSineTime.at( j ) * std::sin( frequencies.at( j ) * ( testTime - referenceEpoch ) );
                }

                // Compare corrections against expected variations.
                for( unsigned int i = 0; i < 5; i++ )
                {
                    for( unsigned int j = 0; j < 5; j++ )
                    {
                        BOOST_CHECK_SMALL( calculatedCosineCoefficientCorrections( i, j ) - expectedCosineCoefficientsCorrections( i, j ),
                                           1.0E-18 );
                        BOOST_CHECK_SMALL( calculatedSineCoefficientCorrections( i, j ) - expectedSineCoefficientsCorrections( i, j ),
                                           1.0E-18 );
                    }
                }
            }
            else
            {
                timeDependentGravityField->update( referenceEpoch );

                // Get corrections at current time step from timeDependentGravityField
                Eigen::MatrixXd perturbedCosineCoefficients = timeDependentGravityField->getCosineCoefficients( );
                Eigen::MatrixXd perturbedSineCoefficients = timeDependentGravityField->getSineCoefficients( );
                Eigen::MatrixXd calculatedCosineCoefficientCorrections =
                        ( perturbedCosineCoefficients - nominalCosineCoefficients ).block( 0, 0, 5, 5 );
                Eigen::MatrixXd calculatedSineCoefficientCorrections =
                        ( perturbedSineCoefficients - nominalSineCoefficients ).block( 0, 0, 5, 5 );

                Eigen::MatrixXd expectedCosineCoefficientsCorrections = Eigen::MatrixXd::Zero( 5, 5 );
                Eigen::MatrixXd expectedSineCoefficientsCorrections = Eigen::MatrixXd::Zero( 5, 5 );

                expectedCosineCoefficientsCorrections.block( minimumDegree, minimumOrder, 2, 3 ) = cosineShAmplitudesCosineTime.at( 0 );
                expectedSineCoefficientsCorrections.block( minimumDegree, minimumOrder, 2, 3 ) = sineShAmplitudesCosineTime.at( 0 );

                // Compare corrections against expected variations.
                for( unsigned int i = 0; i < 5; i++ )
                {
                    for( unsigned int j = 0; j < 5; j++ )
                    {
                        BOOST_CHECK_SMALL( calculatedCosineCoefficientCorrections( i, j ) - expectedCosineCoefficientsCorrections( i, j ),
                                           1.0E-18 );
                        BOOST_CHECK_SMALL( calculatedSineCoefficientCorrections( i, j ) - expectedSineCoefficientsCorrections( i, j ),
                                           1.0E-18 );
                    }
                }
            }
        }
    }
}

void getPolynomialGravityFieldVariationSettings( std::map< int, Eigen::MatrixXd >& cosineAmplitudes,
                                                 std::map< int, Eigen::MatrixXd >& sineAmplitudes,
                                                 double& referenceEpoch,
                                                 int& minimumDegree,
                                                 int& minimumOrder,
                                                 const int type = 0 )
{
    cosineAmplitudes[ 1 ] =
            ( Eigen::MatrixXd( 2, 3 ) << 2.0E-8, -2.0E-8, 4.2E-9, 4.0E-8, -3.7E-9, 2.9E-7 ).finished( ) / physical_constants::JULIAN_YEAR;
    sineAmplitudes[ 1 ] =
            ( Eigen::MatrixXd( 2, 3 ) << 0.0, 6.0E-8, 5.2E-9, 0.0, 7.4E-8, 5.6E-7 ).finished( ) / physical_constants::JULIAN_YEAR;
    referenceEpoch = 5.0E8;
    minimumDegree = 2;
    minimumOrder = 0;
}

// Check polynomial values and rates for constant, linear, quadratic, and cubic terms, including
// the reference epoch where constant terms must contribute zero rate.
BOOST_AUTO_TEST_CASE( testPolynomialGravityFieldVariations )
{
    using namespace tudat::simulation_setup;

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define properties of nominal field
    double gravitationalParameter = getBodyGravitationalParameter( "Jupiter" );
    double referenceRadius = getAverageRadius( "Jupiter" );
    Eigen::MatrixXd nominalCosineCoefficients;
    Eigen::MatrixXd nominalSineCoefficients;
    getNominalJupiterGravityField( nominalCosineCoefficients, nominalSineCoefficients );
    double testTime = 5.0 * physical_constants::JULIAN_YEAR;
    //    for( int k = 0; k < 4; k++ )
    {
        std::map< int, Eigen::MatrixXd > cosineAmplitudes;
        std::map< int, Eigen::MatrixXd > sineAmplitudes;
        double referenceEpoch;
        int minimumDegree;
        int minimumOrder;
        getPolynomialGravityFieldVariationSettings( cosineAmplitudes, sineAmplitudes, referenceEpoch, minimumDegree, minimumOrder, 0 );

        cosineAmplitudes[ 0 ] = cosineAmplitudes.at( 1 ) * physical_constants::JULIAN_YEAR;
        sineAmplitudes[ 0 ] = sineAmplitudes.at( 1 ) * physical_constants::JULIAN_YEAR;
        cosineAmplitudes[ 2 ] = cosineAmplitudes.at( 1 ) / physical_constants::JULIAN_YEAR;
        sineAmplitudes[ 3 ] = sineAmplitudes.at( 1 ) / ( physical_constants::JULIAN_YEAR * physical_constants::JULIAN_YEAR );

        std::shared_ptr< PolynomialGravityFieldVariationsSettings > variationSettings =
                std::dynamic_pointer_cast< PolynomialGravityFieldVariationsSettings >( polynomialGravityFieldVariationsSettings(
                        cosineAmplitudes, sineAmplitudes, referenceEpoch, minimumDegree, minimumOrder ) );

        // Create time-varying gravity field.
        std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                std::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                        gravitationalParameter, referenceRadius, nominalCosineCoefficients, nominalSineCoefficients );

        SystemOfBodies dummySystem;
        dummySystem.createEmptyBody( "Jupiter" );
        dummySystem.getBody( "Jupiter" )->setGravityFieldModel( timeDependentGravityField );
        std::vector< std::shared_ptr< GravityFieldVariationSettings > > variationSettingsList;
        variationSettingsList.push_back( variationSettings );

        std::shared_ptr< gravitation::GravityFieldVariationsSet > variations =
                createGravityFieldModelVariationsSet( "Jupiter", dummySystem, variationSettingsList );
        timeDependentGravityField->updateCorrectionFunctions( );

        const auto polynomialVariation = variations->getVariationObjects( ).at( 0 );
        // Compare both rate blocks with finite differences on either side of, and at, the reference epoch.
        for( const double timeOffset : { -0.25 * physical_constants::JULIAN_YEAR, 0.0, 0.5 * physical_constants::JULIAN_YEAR } )
        {
            checkGravityFieldVariationDerivative( polynomialVariation, referenceEpoch + timeOffset, 100.0 );
        }
        const auto referenceRates = polynomialVariation->calculateSphericalHarmonicsCorrectionsTimeDerivative( referenceEpoch );
        // At the reference epoch, only the linear terms contribute; constants must not introduce NaNs.
        BOOST_CHECK_SMALL( ( referenceRates.first - cosineAmplitudes.at( 1 ) ).norm( ), 1.0e-25 );
        BOOST_CHECK_SMALL( ( referenceRates.second - sineAmplitudes.at( 1 ) ).norm( ), 1.0e-25 );

        timeDependentGravityField->update( 2.0 * testTime );

        // Calculate variations for current test time.
        timeDependentGravityField->update( testTime );

        // Get corrections at current time step from timeDependentGravityField
        Eigen::MatrixXd perturbedCosineCoefficients = timeDependentGravityField->getCosineCoefficients( );
        Eigen::MatrixXd perturbedSineCoefficients = timeDependentGravityField->getSineCoefficients( );
        Eigen::MatrixXd calculatedCosineCoefficientCorrections =
                ( perturbedCosineCoefficients - nominalCosineCoefficients ).block( 0, 0, 5, 5 );
        Eigen::MatrixXd calculatedSineCoefficientCorrections = ( perturbedSineCoefficients - nominalSineCoefficients ).block( 0, 0, 5, 5 );

        Eigen::MatrixXd expectedCosineCoefficientsCorrections = Eigen::MatrixXd::Zero( 5, 5 );
        Eigen::MatrixXd expectedSineCoefficientsCorrections = Eigen::MatrixXd::Zero( 5, 5 );

        for( auto it : cosineAmplitudes )
        {
            expectedCosineCoefficientsCorrections.block( minimumDegree, minimumOrder, 2, 3 ) +=
                    it.second * std::pow( testTime - referenceEpoch, it.first );
        }

        for( auto it : sineAmplitudes )
        {
            expectedSineCoefficientsCorrections.block( minimumDegree, minimumOrder, 2, 3 ) +=
                    it.second * std::pow( testTime - referenceEpoch, it.first );
        }

        // Compare corrections against expected variations.
        for( unsigned int i = 0; i < 5; i++ )
        {
            for( unsigned int j = 0; j < 5; j++ )
            {
                BOOST_CHECK_SMALL( calculatedCosineCoefficientCorrections( i, j ) - expectedCosineCoefficientsCorrections( i, j ),
                                   1.0E-18 );
                BOOST_CHECK_SMALL( calculatedSineCoefficientCorrections( i, j ) - expectedSineCoefficientsCorrections( i, j ), 1.0E-18 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
