/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"

namespace tudat
{

namespace gravitation
{

//! Function to create constant complex Love number list for a range of degrees and orders.
std::vector< std::vector< std::complex< double > > > getFullLoveNumbersVector(
        const std::complex< double > constantLoveNumber, const int maximumDegree, const int maximumOrder )
{
    // Define list of Love numbers
    std::vector< std::vector< std::complex< double > > > fullLoveNumbersVector;
    fullLoveNumbersVector.resize( maximumDegree - 1 );

    // Iterate over all degrees and orders.
    for( unsigned int i = 0; i <= static_cast< unsigned int >( maximumDegree ) - 2;i++ )
    {
        for( unsigned int j = 0; j <= ( i + 2 ); j++ )
        {
            // Set current Love numbers.
            if( static_cast< int >( j ) <= maximumOrder  )
            {
                fullLoveNumbersVector[ i ].push_back( constantLoveNumber );
            }
            else
            {
                fullLoveNumbersVector[ i ].push_back( std::complex< double >( 0.0, 0.0 ) );
            }
        }
    }

    return fullLoveNumbersVector;
}

//! Function to create constant real Love number list for a range of degrees and orders.
std::vector< std::vector< std::complex< double > > > getFullLoveNumbersVector(
        const double constantLoveNumber, const int maximumDegree, const int maximumOrder )
{
    return getFullLoveNumbersVector( std::complex< double >( constantLoveNumber, 0.0 ), maximumDegree, maximumOrder );
}

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
        const std::vector< std::vector< std::complex< double > > > loveNumbers, const double massRatio,
        const double referenceRadius, const Eigen::Vector3d& relativeBodyFixedPosition,
        const int maximumDegree, const int maximumOrder )
{
    // Initialize results.
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );

    std::complex< double > currentCorrections;

    // Iterate over all requested degrees and orders
    for( int n = 2; n <= maximumDegree; n++ )
    {
        for( int m = 0; ( m <= maximumOrder && m <= n ); m++ )
        {
            // Calculate and set corrections at current degree and order.
            currentCorrections = calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                       loveNumbers.at( n - 2 ).at( m ), massRatio, referenceRadius, relativeBodyFixedPosition, n, m );
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
void BasicSolidBodyTideGravityFieldVariations::setBodyGeometryParameters(
        const int bodyIndex, const double evaluationTime )
{
    // Calculate current state and orientation of deformed body.
    if( bodyIndex == 0 )
    {
        deformedBodyPosition = deformedBodyStateFunction_( evaluationTime ).segment( 0, 3 );
        toDeformedBodyFrameRotation = deformedBodyOrientationFunction_( evaluationTime );
    }

    // Calculate current state of body causing deformation.
    Eigen::Vector3d relativeDeformingBodyPosition = toDeformedBodyFrameRotation * (
                deformingBodyStateFunctions_[ bodyIndex ]( evaluationTime ).segment( 0, 3 ) -
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
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > BasicSolidBodyTideGravityFieldVariations::
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
    radiusRatioPower = radiusRatio * radiusRatio * radiusRatio;

    currentCosineCorrections_.setZero( );
    currentSineCorrections_.setZero( );

    //    // Iterate over all love
    std::complex< double > stokesCoefficientCorrection( 0.0, 0.0 );

    for( unsigned int n = 2; n < loveNumbers_.size( ) + 2; n++ )
    {
        for( unsigned int m = 0; ( m <= n && m < loveNumbers_.at( n - 2 ).size( ) ); m++ )
        {
            updateTidalAmplitudeAndArgument( n, m );

            // Calculate and add coefficients.
            stokesCoefficientCorrection =
                    calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                        loveNumbers_[ n - 2 ][ m ], massRatio, radiusRatioPower, tideAmplitude,
                    tideArgument, n, m );


            currentCosineCorrections_( n - 2, m ) += stokesCoefficientCorrection.real( );
            if( m != 0 )
            {
                currentSineCorrections_( n - 2, m ) -= stokesCoefficientCorrection.imag( );
            }
        }

        // Increment radius ratio power.
        radiusRatioPower *= radiusRatio;
    }

    cTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_  - minimumOrder_ + 1 ) +=
            currentCosineCorrections_;
    sTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_  - minimumOrder_ + 1 ) +=
            currentSineCorrections_;

}

}

}
