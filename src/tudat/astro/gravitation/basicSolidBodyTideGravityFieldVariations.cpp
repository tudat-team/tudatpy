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
std::complex< double > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude( const std::complex< double > loveNumber,
                                                                                          const double massRatio,
                                                                                          const double radiusRatioPowerN,
                                                                                          const double amplitude,
                                                                                          const std::complex< double > tideArgument,
                                                                                          const int degree,
                                                                                          const int order,
                                                                                          const double meanCosineForcing,
                                                                                          const double meanSineForcing )
{
    std::complex< double > tidalForcing = 1 / ( 2.0 * static_cast< double >( degree ) + 1.0 ) * massRatio * radiusRatioPowerN * amplitude *
            basic_mathematics::calculateLegendreGeodesyNormalizationFactor( degree, order ) * std::exp( -tideArgument );

    // subtract mean tidal forcing terms
    tidalForcing.real( tidalForcing.real( ) - meanCosineForcing );
    tidalForcing.imag( tidalForcing.imag( ) + meanSineForcing );  // sine correspons to -i, thus plus.

    // Calculate and return corrections.
    return loveNumber * tidalForcing;
}

//! Function to calculate solid body tide gravity field variations due to single body at single degree and order directly
//! from perturbing body's Cartesian state.
std::complex< double > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude( const std::complex< double > loveNumber,
                                                                                          const double massRatio,
                                                                                          const double referenceRadius,
                                                                                          const Eigen::Vector3d& relativeBodyFixedPosition,
                                                                                          const int degree,
                                                                                          const int order,
                                                                                          const double meanCosineForcing,
                                                                                          const double meanSineForcing )
{
    // Calculate spherical position of perturbing body.
    Eigen::Vector3d sphericalPosition = coordinate_conversions::convertCartesianToSpherical( relativeBodyFixedPosition );

    // Calculate amplitude and argument of tide
    double tideAmplitude = basic_mathematics::computeLegendrePolynomialExplicit(
                                   degree, order, std::sin( mathematical_constants::PI / 2.0 - sphericalPosition.y( ) ) ) *
            basic_mathematics::calculateLegendreGeodesyNormalizationFactor( degree, order );
    std::complex< double > tideArgument = static_cast< double >( order ) * std::complex< double >( 0.0, sphericalPosition( 2 ) );

    std::complex< double > tidalForcing = 1 / ( 2.0 * static_cast< double >( degree ) + 1.0 ) * massRatio *
            std::pow( referenceRadius / sphericalPosition( 0 ), degree + 1 ) * tideAmplitude * std::exp( -tideArgument );

    // subtract mean tidal forcing terms
    tidalForcing.real( tidalForcing.real( ) - meanCosineForcing );
    tidalForcing.imag( tidalForcing.imag( ) + meanSineForcing );  // sine correspons to -i, thus plus.

    // Calculate tidal corrections.
    return loveNumber * tidalForcing;
}

std::map< int, std::vector< double > > generateZeroMeanTermsFromReference(
        const std::map< int, std::vector< std::complex< double > > >& loveNumbersReference )
{
    std::map< int, std::vector< double > > zeroMap;

    for( const auto& it : loveNumbersReference )
    {
        const int degree = it.first;
        const std::size_t numberOfOrders = it.second.size( );
        zeroMap[ degree ] = std::vector< double >( numberOfOrders, 0.0 );
    }

    return zeroMap;
}

// This function takes a map with key=degree and value=a vector of coefficients,
// and maps it onto an n x n matrix.
Eigen::MatrixXd convertSphericalHarmonicCoefficientMapToMatrix( const std::map< int, std::vector< double > >& inputMap, int maximumDegree )
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumDegree + 1 );  // Initialize with zeros

    for( const auto& it : inputMap )
    {
        const int degree = it.first;
        const std::vector< double >& vec = it.second;

        if( degree >= ( maximumDegree + 1 ) )
        {
            continue;  // Skip out-of-bounds rows
        }

        for( std::size_t order = 0; order < vec.size( ) && order < static_cast< std::size_t >( maximumDegree + 1 ); ++order )
        {
            matrix( degree, order ) = vec[ order ];
        }
    }

    return matrix;
}

//! Function to calculate solid body tide gravity field variations due to single body at a set of degrees and orders
//! from perturbing body's Cartesian state.
// (is this overload even used currently?)
//
// meanForcingCosineTerms argument do not have defaults, assuming function is called from within a BasicSolidBodyTideGravityFieldVariations
// object which defaults to zero mean forcing and populates these fields meanForcingSineTerms argument do not have defaults, assuming
// function is called from within a BasicSSolidBodyTideGravityFieldVariations object which defaults to zero mean forcing and populates these
// fields
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
        const std::map< int, std::vector< std::complex< double > > > loveNumbers,
        const double massRatio,
        const double referenceRadius,
        const Eigen::Vector3d& relativeBodyFixedPosition,
        const int maximumDegree,
        const int maximumOrder,
        std::map< int, std::vector< double > > meanCosineForcing,
        std::map< int, std::vector< double > > meanSineForcing )

{
    if( meanCosineForcing.empty( ) )
    {
        meanCosineForcing = generateZeroMeanTermsFromReference( loveNumbers );
    }
    if( meanSineForcing.empty( ) )
    {
        meanSineForcing = generateZeroMeanTermsFromReference( loveNumbers );
    }

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
            currentCorrections = calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude( loveNumbers.at( n ).at( m ),
                                                                                                    massRatio,
                                                                                                    referenceRadius,
                                                                                                    relativeBodyFixedPosition,
                                                                                                    n,
                                                                                                    m,
                                                                                                    meanCosineForcing.at( n ).at( m ),
                                                                                                    meanSineForcing.at( n ).at( m ) );
            cosineCorrections( n, m ) = currentCorrections.real( );
            if( m > 0 )
            {
                sineCorrections( n, m ) = -currentCorrections.imag( );
            }
        }
    }

    return std::make_pair( cosineCorrections, sineCorrections );
}

void SolidBodyTideGravityFieldVariations::initializeTimeDerivative( const int maximumForcingDegree, const int maximumForcingOrder )
{
    maximumForcingDegree_ = maximumForcingDegree;
    maximumForcingOrder_ = maximumForcingOrder;
    const int rows = maximumForcingDegree_ + 1;
    const int columns = maximumForcingOrder_ + 1;
    tidalForcing_ = Eigen::MatrixXcd::Zero( rows, columns );
    tidalForcingRates_ = Eigen::MatrixXcd::Zero( rows, columns );
    firstRecurrenceFactors_ = Eigen::MatrixXd::Zero( rows, columns );
    secondRecurrenceFactors_ = Eigen::MatrixXd::Zero( rows, columns );
    Eigen::MatrixXd normalization = Eigen::MatrixXd::Zero( rows, columns );

    for( int n = 0; n < rows; ++n )
    {
        for( int m = 0; m <= std::min( n, maximumForcingOrder_ ); ++m )
        {
            normalization( n, m ) = basic_mathematics::calculateLegendreGeodesyNormalizationFactor( n, m ) / ( 2.0 * n + 1.0 );
            if( n == m && n > 0 )
            {
                firstRecurrenceFactors_( n, m ) = ( 2.0 * n - 1.0 ) * normalization( n, m ) / normalization( n - 1, m - 1 );
            }
            else if( n > m )
            {
                firstRecurrenceFactors_( n, m ) = ( 2.0 * n - 1.0 ) * normalization( n, m ) / ( ( n - m ) * normalization( n - 1, m ) );
                if( n > m + 1 )
                {
                    secondRecurrenceFactors_( n, m ) = ( n + m - 1.0 ) * normalization( n, m ) / ( ( n - m ) * normalization( n - 2, m ) );
                }
            }
        }
    }
}

template< bool computeTimeDerivative >
void SolidBodyTideGravityFieldVariations::updateTidalForcing( const Eigen::Vector3d& position, const Eigen::Vector3d& velocity )
{
    const double inverseDistanceSquared = 1.0 / position.squaredNorm( );
    const double relativeRadialRate = position.dot( velocity ) * inverseDistanceSquared;
    const double radiusRatio = deformedBodyReferenceRadius_ * std::sqrt( inverseDistanceSquared );
    const double radiusRatioSquared = radiusRatio * radiusRatio;
    const double radiusRatioSquaredRate = -2.0 * relativeRadialRate * radiusRatioSquared;
    const Eigen::Vector3d scaledPosition = deformedBodyReferenceRadius_ * inverseDistanceSquared * position;
    const Eigen::Vector3d scaledVelocity =
            deformedBodyReferenceRadius_ * inverseDistanceSquared * ( velocity - 2.0 * relativeRadialRate * position );
    const std::complex< double > xy( scaledPosition.x( ), -scaledPosition.y( ) );
    const std::complex< double > xyRate( scaledVelocity.x( ), -scaledVelocity.y( ) );
    const double z = scaledPosition.z( );
    const double zRate = scaledVelocity.z( );

    // These are the normalized cosine/sine forcing terms, including mass ratio and 1 / (2n + 1).
    // The Cartesian recurrence and its product-rule derivative also remain regular at the poles.
    tidalForcing_( 0, 0 ) = massRatio * radiusRatio;
    if constexpr( computeTimeDerivative )
    {
        tidalForcingRates_( 0, 0 ) = -relativeRadialRate * tidalForcing_( 0, 0 );
    }
    for( int m = 1; m <= maximumForcingOrder_; ++m )
    {
        const double factor = firstRecurrenceFactors_( m, m );
        tidalForcing_( m, m ) = factor * xy * tidalForcing_( m - 1, m - 1 );
        if constexpr( computeTimeDerivative )
        {
            tidalForcingRates_( m, m ) = factor * ( xyRate * tidalForcing_( m - 1, m - 1 ) + xy * tidalForcingRates_( m - 1, m - 1 ) );
        }
    }
    for( int m = 0; m <= maximumForcingOrder_ && m < maximumForcingDegree_; ++m )
    {
        const double factor = firstRecurrenceFactors_( m + 1, m );
        tidalForcing_( m + 1, m ) = factor * z * tidalForcing_( m, m );
        if constexpr( computeTimeDerivative )
        {
            tidalForcingRates_( m + 1, m ) = factor * ( zRate * tidalForcing_( m, m ) + z * tidalForcingRates_( m, m ) );
        }
        for( int n = m + 2; n <= maximumForcingDegree_; ++n )
        {
            const double a = firstRecurrenceFactors_( n, m );
            const double b = secondRecurrenceFactors_( n, m );
            tidalForcing_( n, m ) = a * z * tidalForcing_( n - 1, m ) - b * radiusRatioSquared * tidalForcing_( n - 2, m );
            if constexpr( computeTimeDerivative )
            {
                tidalForcingRates_( n, m ) = a * ( zRate * tidalForcing_( n - 1, m ) + z * tidalForcingRates_( n - 1, m ) ) -
                        b * ( radiusRatioSquaredRate * tidalForcing_( n - 2, m ) + radiusRatioSquared * tidalForcingRates_( n - 2, m ) );
            }
        }
    }
}

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > SolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivative(
        const double time )
{
    if( !canComputeTimeDerivative_ )
    {
        return GravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivative( time );
    }

    Eigen::MatrixXd cosineRates = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineRates = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    const Eigen::Matrix3d rotationRate = rotationDerivativeFunction_( time );

    for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); ++i )
    {
        setBodyGeometryParameters( i, time );
        const Eigen::Vector3d relativeBodyFixedVelocity = toDeformedBodyFrameRotation * relativeDeformingBodyState_.tail< 3 >( ) +
                rotationRate * relativeDeformingBodyState_.head< 3 >( );
        updateTidalForcing< true >( relativeDeformingBodyFixedPosition_, relativeBodyFixedVelocity );
        addTidalCorrectionTimeDerivatives( cosineRates, sineRates );
    }
    return std::make_pair( cosineRates, sineRates );
}

void BasicSolidBodyTideGravityFieldVariations::addTidalCorrectionTimeDerivatives( Eigen::MatrixXd& cosineRates, Eigen::MatrixXd& sineRates )
{
    // Constant mean-forcing offsets have zero derivative.
    for( const auto& degree : loveNumbers_ )
    {
        const int n = degree.first;
        for( unsigned int m = 0; m < degree.second.size( ) && m <= static_cast< unsigned int >( n ); ++m )
        {
            const std::complex< double > rate = degree.second[ m ] * tidalForcingRates_( n, m );
            cosineRates( n - minimumDegree_, m ) += rate.real( );
            if( m != 0 )
            {
                sineRates( n - minimumDegree_, m ) -= rate.imag( );
            }
        }
    }
}

void ModeCoupledSolidBodyTideGravityFieldVariations::addTidalCorrectionTimeDerivatives( Eigen::MatrixXd& cosineRates,
                                                                                        Eigen::MatrixXd& sineRates )
{
    for( const auto& forcing : loveNumbers_ )
    {
        const std::complex< double > forcingRate = tidalForcingRates_( forcing.first.first, forcing.first.second );
        for( const auto& response : forcing.second )
        {
            const int row = response.first.first - minimumDegree_;
            const int order = response.first.second;
            const std::complex< double > rate = response.second * forcingRate;
            cosineRates( row, order ) += rate.real( );
            if( order != 0 )
            {
                sineRates( row, order ) -= rate.imag( );
            }
        }
    }
}

//! Sets current properties (mass state) of body involved in tidal deformation.
void SolidBodyTideGravityFieldVariations::setBodyGeometryParameters( const int bodyIndex, const double evaluationTime )
{
    // Calculate current state and orientation of deformed body.
    if( bodyIndex == 0 )
    {
        const Eigen::Vector6d deformedBodyState = deformedBodyStateFunction_( evaluationTime );
        deformedBodyPosition = deformedBodyState.head< 3 >( );
        deformedBodyVelocity_ = deformedBodyState.tail< 3 >( );
        toDeformedBodyFrameRotation = deformedBodyOrientationFunction_( evaluationTime );
        inverseDeformedBodyMass_ = 1.0 / deformedBodyMass_( );
    }

    // Calculate current state of body causing deformation.
    relativeDeformingBodyState_ = deformingBodyStateFunctions_[ bodyIndex ]( evaluationTime );
    relativeDeformingBodyState_.head< 3 >( ) -= deformedBodyPosition;
    relativeDeformingBodyState_.tail< 3 >( ) -= deformedBodyVelocity_;
    relativeDeformingBodyFixedPosition_ = toDeformedBodyFrameRotation * relativeDeformingBodyState_.head< 3 >( );
    Eigen::Vector3d relativeDeformingBodySphericalPosition =
            coordinate_conversions::convertCartesianToSpherical( relativeDeformingBodyFixedPosition_ );
    massRatio = deformingBodyMasses_[ bodyIndex ]( ) * inverseDeformedBodyMass_;

    // Set geometric parameters of body causing deformation.
    radiusRatio = deformedBodyReferenceRadius_ / relativeDeformingBodySphericalPosition.x( );
    iLongitude = mathematical_constants::COMPLEX_I * relativeDeformingBodySphericalPosition.z( );
    sineOfLatitude = std::sin( mathematical_constants::PI / 2.0 - relativeDeformingBodySphericalPosition.y( ) );
}

//! Function for calculating spherical harmonic coefficient corrections.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > SolidBodyTideGravityFieldVariations::calculateBasicSphericalHarmonicsCorrections(
        const double time )
{
    // Initialize corrections to zero.
    Eigen::MatrixXd cTermCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sTermCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );

    // Iterate over all bodies causing deformation and calculate and add associated corrections
    for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
    {
        setBodyGeometryParameters( i, time );

        updateTidalForcing< false >( relativeDeformingBodyFixedPosition_, Eigen::Vector3d::Zero( ) );

        // Calculate all correction functions.
        for( unsigned int j = 0; j < correctionFunctions.size( ); j++ )
        {
            correctionFunctions[ j ]( cTermCorrections, sTermCorrections );
        }
    }

    return std::make_pair( cTermCorrections, sTermCorrections );
}

//! Calculates basic solid body gravity field corrections due to single body.
void BasicSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections( Eigen::MatrixXd& cTermCorrections,
                                                                                 Eigen::MatrixXd& sTermCorrections )
{
    //    // Initialize power of radiusRatio^(N+1) (calculation starts at N=2)

    currentCosineCorrections_.setZero( );
    currentSineCorrections_.setZero( );

    //    // Iterate over all love
    std::complex< double > stokesCoefficientCorrection( 0.0, 0.0 );

    for( const auto& loveNumberIt : loveNumbers_ )
    {
        unsigned int n = static_cast< unsigned int >( loveNumberIt.first );
        const auto& meanCosine = meanForcingCosineTerms_.at( n );
        const auto& meanSine = meanForcingSineTerms_.at( n );

        for( unsigned int m = 0; ( m <= n && m < loveNumberIt.second.size( ) ); m++ )
        {
            // Values and rates use the same normalized forcing recurrence.
            stokesCoefficientCorrection =
                    loveNumberIt.second[ m ] * ( tidalForcing_( n, m ) - std::complex< double >( meanCosine[ m ], -meanSine[ m ] ) );

            currentCosineCorrections_( n - 2, m ) += stokesCoefficientCorrection.real( );
            if( m != 0 )
            {
                currentSineCorrections_( n - 2, m ) -= stokesCoefficientCorrection.imag( );
            }
        }
    }

    cTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_ - minimumOrder_ + 1 ) += currentCosineCorrections_;
    sTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_ - minimumOrder_ + 1 ) += currentSineCorrections_;
}

int getBasicTideMaximumDegree( const std::map< int, std::vector< std::complex< double > > >& loveNumbers )
{
    if( loveNumbers.empty( ) )
    {
        throw std::runtime_error( "Error creating basic tidal variation: Love numbers are empty." );
    }
    return loveNumbers.rbegin( )->first;
}

int getModeCoupledMaximumResponseDegree( const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > >& loveNumbers )
{
    if( loveNumbers.empty( ) )
    {
        throw std::runtime_error( "Error creating mode-coupled tidal variation: Love numbers are empty." );
    }
    int maximumDegree = 0;
    for( auto it : loveNumbers )
    {
        for( auto it2 : it.second )
        {
            if( it2.first.first > maximumDegree )
            {
                maximumDegree = it2.first.first;
            }
        }
    }
    return maximumDegree;
}

int getModeCoupledMaximumResponseOrder( const std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > >& loveNumbers )
{
    int maximumOrder = 0;
    for( auto it : loveNumbers )
    {
        for( auto it2 : it.second )
        {
            if( it2.first.second > maximumOrder )
            {
                maximumOrder = it2.first.second;
            }
        }
    }
    return maximumOrder;
}

//! Calculates basic solid body gravity field corrections due to single body.
void ModeCoupledSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections( Eigen::MatrixXd& cTermCorrections,
                                                                                       Eigen::MatrixXd& sTermCorrections )
{
    //    // Initialize power of radiusRatio^(N+1) (calculation starts at N=2)

    currentCosineCorrections_.setZero( );
    currentSineCorrections_.setZero( );

    //    // Iterate over all love
    std::complex< double > stokesCoefficientCorrection( 0.0, 0.0 );

    for( const auto& loveNumberIt : loveNumbers_ )
    {
        // Retrieve forcing degree/order
        std::pair< int, int > forcingDegreeOrder = loveNumberIt.first;
        int n = forcingDegreeOrder.first;
        int m = forcingDegreeOrder.second;

        const std::complex< double > unityLoveNumberStokesCoefficientCorrection = tidalForcing_( n, m );

        // Compute response at required degrees/orders with required love numbers
        for( const auto& responseIt : loveNumberIt.second )
        {
            int nResponse = responseIt.first.first;
            int mResponse = responseIt.first.second;
            stokesCoefficientCorrection = std::complex< double >( responseIt.second, 0.0 ) * unityLoveNumberStokesCoefficientCorrection;
            currentCosineCorrections_( nResponse - 2, mResponse ) += stokesCoefficientCorrection.real( );
            if( mResponse != 0 )
            {
                currentSineCorrections_( nResponse - 2, mResponse ) -= stokesCoefficientCorrection.imag( );
            }
        }
    }

    cTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_ - minimumOrder_ + 1 ) += currentCosineCorrections_;
    sTermCorrections.block( 0, 0, maximumDegree_ - minimumDegree_ + 1, maximumOrder_ - minimumOrder_ + 1 ) += currentSineCorrections_;
}

}  // namespace gravitation

}  // namespace tudat
