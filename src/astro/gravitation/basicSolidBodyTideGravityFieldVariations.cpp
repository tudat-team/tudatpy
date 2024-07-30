/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{

namespace gravitation
{

//! Function to calculate solid body tide gravity field variation due to single body at single
//! degree and order.
std::complex< double > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
        const std::complex< double > loveNumber, const double massRatio,
        const double radiusRatioPowerN,
        const double amplitude, const std::complex< double > tideArgument,
        const int degree, const int order )
{
    // Calculate and return corrections.
    return loveNumber / ( 2.0 * static_cast< double >( degree ) + 1.0 ) *
            massRatio * radiusRatioPowerN *
            amplitude * basic_mathematics::calculateLegendreGeodesyNormalizationFactor(
                degree, order ) * std::exp( -tideArgument );
}

//! Function to calculate solid body tide gravity field variations due to single body at single degree and order directly
//! from perturbing body's Cartesian state.
std::complex< double > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
        const std::complex< double > loveNumber, const double massRatio,
        const double referenceRadius, const Eigen::Vector3d& relativeBodyFixedPosition, const int degree, const int order )
{
    // Calculate spherical position of perturbing body.
    Eigen::Vector3d sphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                relativeBodyFixedPosition );

    // Calculate amplitude and argument of tide
    double tideAmplitude = basic_mathematics::computeLegendrePolynomialExplicit(
                degree, order, std::sin( mathematical_constants::PI / 2.0 - sphericalPosition.y( ) ) ) *
            basic_mathematics::calculateLegendreGeodesyNormalizationFactor(
                            degree, order );
    std::complex< double > tideArgument = static_cast< double >( order ) * std::complex< double >( 0.0, sphericalPosition( 2 ) );

    // Calculate tidal corrections.
    return loveNumber / ( 2.0 * static_cast< double >( degree ) + 1.0 ) *
            massRatio * std::pow( referenceRadius / sphericalPosition( 0 ), degree + 1 ) *
            tideAmplitude * std::exp( -tideArgument );
}

//! Function to calculate solid body tide gravity field variations due to single body at a set of degrees and orders
//! from perturbing body's Cartesian state.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
        const std::map< int, std::vector< std::complex< double > > > loveNumbers, const double massRatio,
        const double referenceRadius, const Eigen::Vector3d& relativeBodyFixedPosition,
        const int maximumDegree, const int maximumOrder )
{
    // Initialize results.
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );

    std::complex< double > currentCorrections;

    // Iterate over all requested degrees and orders
    for( auto loveNumberIt : loveNumbers )
    {
        int n = loveNumberIt.first;
        for( int m = 0; ( m <= maximumOrder && m <= n ); m++ )
        {
            // Calculate and set corrections at current degree and order.
            currentCorrections = calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                       loveNumbers.at( n ).at( m ), massRatio, referenceRadius, relativeBodyFixedPosition, n, m );
            cosineCorrections( n, m ) = currentCorrections.real( );
            if( m > 0 )
            {
                sineCorrections( n, m ) = -currentCorrections.imag( );
            }
        }
    }

    return std::make_pair( cosineCorrections, sineCorrections );
}


//! Sets current properties (mass state) of body involved in tidal deformation.
void SolidBodyTideGravityFieldVariations::setBodyGeometryParameters(
        const int bodyIndex, const double evaluationTime )
{
    // Calculate current state and orientation of deformed body.
    if( bodyIndex == 0 )
    {
        deformedBodyPosition = std::move( deformedBodyStateFunction_( evaluationTime ) ).segment( 0, 3 );
        toDeformedBodyFrameRotation = deformedBodyOrientationFunction_( evaluationTime );
    }

    // Calculate current state of body causing deformation.
    Eigen::Vector3d relativeDeformingBodyPosition = toDeformedBodyFrameRotation * (
                std::move( deformingBodyStateFunctions_[ bodyIndex ]( evaluationTime ) ).segment( 0, 3 ) -
            deformedBodyPosition );
    Eigen::Vector3d  relativeDeformingBodySphericalPosition =coordinate_conversions::
            convertCartesianToSpherical( relativeDeformingBodyPosition );

    // Set geometric parameters of body causing deformation.
    radiusRatio = deformedBodyReferenceRadius_ / relativeDeformingBodySphericalPosition.x( );
    iLongitude = mathematical_constants::COMPLEX_I * relativeDeformingBodySphericalPosition.z( );
    sineOfLatitude = std::sin( mathematical_constants::PI / 2.0 -
                               relativeDeformingBodySphericalPosition.y( ) );
}

//! Function for calculating spherical harmonic coefficient corrections.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > SolidBodyTideGravityFieldVariations::
calculateBasicSphericalHarmonicsCorrections(
        const double time )
{
    // Initialize corrections to zero.
    Eigen::MatrixXd cTermCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sTermCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );

    // Iterate over all bodies causing deformation and calculate and add associated corrections
    for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
    {
        setBodyGeometryParameters( i, time );

        // Calculate properties of currently considered body
        massRatio = deformingBodyMasses_[ i ]( ) / deformedBodyMass_( );

        // Calculate all correction functions.
        for( unsigned int j = 0; j < correctionFunctions.size( ); j++ )
        {
            correctionFunctions[ j ]( cTermCorrections, sTermCorrections );
        }
    }

    return std::make_pair( cTermCorrections, sTermCorrections );
}

//! Calculates basic solid body gravity field corrections due to single body.
void BasicSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections(
        Eigen::MatrixXd& cTermCorrections,
        Eigen::MatrixXd& sTermCorrections )
{
//    // Initialize power of radiusRatio^(N+1) (calculation starts at N=2)

    currentCosineCorrections_.setZero( );
    currentSineCorrections_.setZero( );

    //    // Iterate over all love
    std::complex< double > stokesCoefficientCorrection( 0.0, 0.0 );

    for( auto loveNumberIt : loveNumbers_ )
    {
        unsigned int n = static_cast< unsigned int >( loveNumberIt.first );
        radiusRatioPower = basic_mathematics::raiseToIntegerPower( radiusRatio, n + 1 );

        for( unsigned int m = 0; ( m <= n && m < loveNumberIt.second.size( ) ); m++ )
        {
            updateTidalAmplitudeAndArgument( n, m );

            // Calculate and add coefficients.
            stokesCoefficientCorrection =
                    calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                        loveNumbers_[ n ][ m ], massRatio, radiusRatioPower, tideAmplitude,
                    tideArgument, n, m );


            currentCosineCorrections_( n - 2, m ) += stokesCoefficientCorrection.real( );
            if( m != 0 )
            {
                currentSineCorrections_( n - 2, m ) -= stokesCoefficientCorrection.imag( );
            }
        }
    }

    cTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_  - minimumOrder_ + 1 ) +=
            currentCosineCorrections_;
    sTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_  - minimumOrder_ + 1 ) +=
            currentSineCorrections_;

}


int getModeCoupledMaximumResponseDegree(
    const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > >& loveNumbers )
{
    int maximumDegree = 0;
    for( auto it : loveNumbers )
    {
        for( auto it2 : it.second )
        {
            if( it2.first.first > maximumDegree )
            {
                maximumDegree =  it2.first.first;
            }
        }
    }
    return maximumDegree;
}

int getModeCoupledMaximumResponseOrder(
    const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > >& loveNumbers )
{
    int maximumOrder = 0;
    for( auto it : loveNumbers )
    {
        for( auto it2 : it.second )
        {
            if( it2.first.second > maximumOrder )
            {
                maximumOrder =  it2.first.second;
            }
        }
    }
    return maximumOrder;
}


//! Calculates basic solid body gravity field corrections due to single body.
void ModeCoupledSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections(
    Eigen::MatrixXd& cTermCorrections,
    Eigen::MatrixXd& sTermCorrections )
{
//    // Initialize power of radiusRatio^(N+1) (calculation starts at N=2)

    currentCosineCorrections_.setZero( );
    currentSineCorrections_.setZero( );

    //    // Iterate over all love
    std::complex< double > stokesCoefficientCorrection( 0.0, 0.0 );

    for( auto loveNumberIt : loveNumbers_ )
    {
        // Retrieve forcing degree/order
        std::pair< int, int > forcingDegreeOrder = loveNumberIt.first;
        int n = forcingDegreeOrder.first;
        int m = forcingDegreeOrder.second;

        // Compute forcing quantities
        radiusRatioPower = basic_mathematics::raiseToIntegerPower( radiusRatio, n + 1 );
        updateTidalAmplitudeAndArgument( n, m );

        // Compute response for unity Love umber
        std::complex< double > unityLoveNumberStokesCoefficientCorrection =
                calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                    std::complex< double >( 1.0, 0.0 ), massRatio, radiusRatioPower, tideAmplitude,
                    tideArgument, n, m );

        // Compute response at required degrees/orders with required love numbers
        for( auto responseIt : loveNumberIt.second )
        {
            int nResponse = responseIt.first.first;
            int mResponse = responseIt.first.second;
            stokesCoefficientCorrection = std::complex< double >( responseIt.second, 0.0 ) *
                unityLoveNumberStokesCoefficientCorrection;
            currentCosineCorrections_( nResponse - 2, mResponse ) += stokesCoefficientCorrection.real( );
            if ( mResponse != 0 )
            {
                currentSineCorrections_( nResponse - 2, mResponse ) -= stokesCoefficientCorrection.imag( );
            }
        }
    }

    cTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_  - minimumOrder_ + 1 ) +=
        currentCosineCorrections_;
    sTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_  - minimumOrder_ + 1 ) +=
        currentSineCorrections_;

}

}

}
