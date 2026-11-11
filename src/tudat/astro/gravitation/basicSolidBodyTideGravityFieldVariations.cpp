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
#include "tudat/astro/ephemerides/rotationalEphemeris.h"

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

//! Function to compute the time derivative of the tidal forcing of a single degree and order, for unit Love number
std::complex< double > SolidBodyTideGravityFieldVariations::calculateTidalForcingTimeDerivative( const int degree,
                                                                                                 const int order,
                                                                                                 const double deformingBodyMassRatio,
                                                                                                 const Eigen::Vector3d& sphericalPosition,
                                                                                                 const Eigen::Vector3d& sphericalVelocity )
{
    // The forcing (R/r)^(n+1) * Pnm * exp( -i * m * longitude ) is the potential of a single spherical harmonic term (with
    // mu = R), so that its time derivative is the gradient of this term, multiplied by the body-fixed velocity (here, spherical
    // gradient multiplied by sphericalVelocity). The spherical harmonics cache must be updated to the current position.
    // Near the poles, the gradient is computed without the (singular) Legendre polynomial derivative, as for the spherical harmonic
    // acceleration (same colatitude cutoff)
    Eigen::Vector3d cosineTermSphericalGradient, sineTermSphericalGradient;
    if( std::fabs( mathematical_constants::PI / 2.0 - std::fabs( sphericalPosition( 1 ) ) ) < 1.0E-6 )
    {
        cosineTermSphericalGradient = basic_mathematics::computePotentialGradientNearPole(
                sphericalPosition, 1.0, degree, order, 1.0, 0.0, sphericalHarmonicsCache_ );
        sineTermSphericalGradient = basic_mathematics::computePotentialGradientNearPole(
                sphericalPosition, 1.0, degree, order, 0.0, 1.0, sphericalHarmonicsCache_ );
    }
    else
    {
        const double legendrePolynomial = sphericalHarmonicsCache_.getLegendreCache( ).getLegendrePolynomial( degree, order );
        const double legendrePolynomialDerivative =
                sphericalHarmonicsCache_.getLegendreCache( ).getLegendrePolynomialDerivative( degree, order );
        cosineTermSphericalGradient = basic_mathematics::computePotentialGradient( sphericalPosition,
                                                                                   1.0,
                                                                                   degree,
                                                                                   order,
                                                                                   1.0,
                                                                                   0.0,
                                                                                   legendrePolynomial,
                                                                                   legendrePolynomialDerivative,
                                                                                   sphericalHarmonicsCache_ );
        sineTermSphericalGradient = basic_mathematics::computePotentialGradient( sphericalPosition,
                                                                                 1.0,
                                                                                 degree,
                                                                                 order,
                                                                                 0.0,
                                                                                 1.0,
                                                                                 legendrePolynomial,
                                                                                 legendrePolynomialDerivative,
                                                                                 sphericalHarmonicsCache_ );
    }
    return deformingBodyMassRatio / ( 2.0 * static_cast< double >( degree ) + 1.0 ) *
            std::complex< double >( cosineTermSphericalGradient.dot( sphericalVelocity ),
                                    -sineTermSphericalGradient.dot( sphericalVelocity ) );
}

//! Function to compute the state of a body causing deformation w.r.t. the deformed body, in the deformed body-fixed frame
Eigen::Vector6d SolidBodyTideGravityFieldVariations::getBodyFixedRelativeState( const int bodyIndex, const double time )
{
    return ephemerides::transformStateToFrameFromRotations< double >(
            deformingBodyStateFunctions_[ bodyIndex ]( time ) - deformedBodyStateFunction_( time ),
            deformedBodyOrientationFunction_( time ),
            rotationDerivativeFunction_( time ) );
}

//! Function to update the spherical harmonics cache to the body-fixed position of a body causing deformation
Eigen::Vector3d SolidBodyTideGravityFieldVariations::updateSphericalHarmonicsCache( const Eigen::Vector3d& bodyFixedPosition )
{
    // Convert to spherical position, with latitude instead of colatitude
    Eigen::Vector3d sphericalPosition = coordinate_conversions::convertCartesianToSpherical( bodyFixedPosition );
    sphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - sphericalPosition( 1 );
    sphericalHarmonicsCache_.update(
            sphericalPosition( 0 ), std::sin( sphericalPosition( 1 ) ), sphericalPosition( 2 ), deformedBodyReferenceRadius_ );
    return sphericalPosition;
}

//! Function for calculating time derivatives of spherical harmonic coefficient corrections due to selected bodies.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > SolidBodyTideGravityFieldVariations::calculateTimeDerivativeOfCorrections(
        const double time,
        const std::vector< int >& bodyIndices,
        const int maximumDegree,
        const int angularVelocityComponent )
{
    Eigen::MatrixXd cosineRates = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineRates = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    if( bodyIndices.size( ) == 0 )
    {
        return std::make_pair( cosineRates, sineRates );
    }

    // Retrieve state, orientation, rotation derivative and mass of deformed body
    const Eigen::Vector6d deformedBodyState = deformedBodyStateFunction_( time );
    const Eigen::Quaterniond deformedBodyOrientation = deformedBodyOrientationFunction_( time );
    const Eigen::Matrix3d deformedBodyRotationDerivative = rotationDerivativeFunction_( time );
    const double deformedBodyMass = deformedBodyMass_( );

    for( const int bodyIndex : bodyIndices )
    {
        // Compute state of body causing deformation w.r.t. deformed body, in deformed body-fixed frame
        const Eigen::Vector6d bodyFixedRelativeState = ephemerides::transformStateToFrameFromRotations< double >(
                deformingBodyStateFunctions_[ bodyIndex ]( time ) - deformedBodyState,
                deformedBodyOrientation,
                deformedBodyRotationDerivative );
        const Eigen::Vector3d bodyFixedPosition = bodyFixedRelativeState.segment( 0, 3 );
        const double currentMassRatio = deformingBodyMasses_[ bodyIndex ]( ) / deformedBodyMass;

        // The time derivatives are linear in the body-fixed velocity, of which the part due to the rotation of the deformed body
        // is -omega x r. The partial w.r.t. an angular velocity component is obtained by replacing the velocity with its partial.
        const Eigen::Vector3d bodyFixedVelocity = ( angularVelocityComponent >= 0 )
                ? Eigen::Vector3d( bodyFixedPosition.cross( Eigen::Vector3d::Unit( angularVelocityComponent ) ) )
                : Eigen::Vector3d( bodyFixedRelativeState.segment( 3, 3 ) );

        // Update spherical harmonics cache, and premultiply velocity by transpose of spherical-to-Cartesian gradient matrix, once for
        // all degrees and orders
        const Eigen::Vector3d sphericalPosition = updateSphericalHarmonicsCache( bodyFixedPosition );
        const Eigen::Vector3d sphericalVelocity =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( bodyFixedPosition ).transpose( ) * bodyFixedVelocity;

        addTidalCorrectionTimeDerivatives(
                [ & ]( const int degree, const int order ) {
                    return calculateTidalForcingTimeDerivative( degree, order, currentMassRatio, sphericalPosition, sphericalVelocity );
                },
                cosineRates,
                sineRates,
                maximumDegree );
    }
    return std::make_pair( cosineRates, sineRates );
}

//! Function for calculating time derivatives of spherical harmonic coefficient corrections.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > SolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivative(
        const double time )
{
    if( rotationDerivativeFunction_ == nullptr )
    {
        return GravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivative( time );
    }

    std::vector< int > bodyIndices;
    for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
    {
        bodyIndices.push_back( i );
    }
    return calculateTimeDerivativeOfCorrections( time, bodyIndices, maximumDegree_ );
}

//! Function for calculating time derivatives of spherical harmonic coefficient corrections due to a single body causing deformation
std::pair< Eigen::MatrixXd, Eigen::MatrixXd >
SolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivativeDueToBody( const double time,
                                                                                                    const int bodyIndex,
                                                                                                    const int maximumDegree )
{
    if( rotationDerivativeFunction_ == nullptr )
    {
        return GravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivative( time );
    }
    return calculateTimeDerivativeOfCorrections( time, { bodyIndex }, maximumDegree );
}

//! Function for calculating partials of time derivatives of coefficient corrections w.r.t. angular velocity of deformed body
std::pair< Eigen::MatrixXd, Eigen::MatrixXd >
SolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivativePartialWrtAngularVelocity(
        const double time,
        const int angularVelocityComponent,
        const int maximumDegree )
{
    std::vector< int > bodyIndices;
    if( rotationDerivativeFunction_ != nullptr )
    {
        for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
        {
            bodyIndices.push_back( i );
        }
    }
    return calculateTimeDerivativeOfCorrections( time, bodyIndices, maximumDegree, angularVelocityComponent );
}

//! Function to compute the partial of the time derivative of the tidal forcing w.r.t. the body-fixed position of the deforming body
std::complex< double > SolidBodyTideGravityFieldVariations::calculateTidalForcingTimeDerivativePartialWrtPosition(
        const int degree,
        const int order,
        const double deformingBodyMassRatio,
        const Eigen::Vector6d& bodyFixedRelativeState,
        const Eigen::Vector3d& bodyFixedPositionVariation,
        const std::function< Eigen::Matrix3d( const int, const int, const bool, const Eigen::Vector3d& ) >& forcingHessianFunction )
{
    // The time derivative is the gradient of the single spherical harmonic forcing term, multiplied by the body-fixed velocity. Its
    // partial w.r.t. position is obtained from the Hessian of this term.
    const Eigen::Vector3d bodyFixedPosition = bodyFixedRelativeState.segment( 0, 3 );
    const Eigen::Matrix3d cosineTermHessian = forcingHessianFunction( degree, order, true, bodyFixedPosition );
    const Eigen::Matrix3d sineTermHessian = forcingHessianFunction( degree, order, false, bodyFixedPosition );
    return deformingBodyMassRatio / ( 2.0 * static_cast< double >( degree ) + 1.0 ) *
            std::complex< double >( ( cosineTermHessian * bodyFixedPositionVariation ).dot( bodyFixedRelativeState.segment( 3, 3 ) ),
                                    -( sineTermHessian * bodyFixedPositionVariation ).dot( bodyFixedRelativeState.segment( 3, 3 ) ) );
}

//! Function for calculating partials of coefficient corrections w.r.t. orientation of deformed body
std::vector< std::pair< Eigen::MatrixXd, Eigen::MatrixXd > >
SolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrectionsPartialWrtOrientation(
        const double time,
        const std::vector< Eigen::Matrix3d >& rotationToBodyFixedFramePartials,
        const int maximumDegree )
{
    std::vector< std::pair< Eigen::MatrixXd, Eigen::MatrixXd > > correctionPartials(
            rotationToBodyFixedFramePartials.size( ),
            std::make_pair( Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ ),
                            Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ ) ) );
    if( rotationDerivativeFunction_ == nullptr )
    {
        return correctionPartials;
    }

    for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
    {
        const Eigen::Vector3d bodyFixedPosition = getBodyFixedRelativeState( i, time ).segment( 0, 3 );
        const Eigen::Vector3d inertialRelativePosition =
                ( deformingBodyStateFunctions_[ i ]( time ) - deformedBodyStateFunction_( time ) ).segment( 0, 3 );
        const double currentMassRatio = deformingBodyMasses_[ i ]( ) / deformedBodyMass_( );
        const Eigen::Vector3d sphericalPosition = updateSphericalHarmonicsCache( bodyFixedPosition );
        const Eigen::Matrix3d gradientMatrixTranspose =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( bodyFixedPosition ).transpose( );

        // The corrections are linear in the forcing, of which the partial w.r.t. body-fixed position is given by its gradient. The
        // partial w.r.t. orientation is therefore obtained from the time derivative, with the velocity replaced by the position partial
        for( unsigned int j = 0; j < rotationToBodyFixedFramePartials.size( ); j++ )
        {
            const Eigen::Vector3d sphericalPositionPartial =
                    gradientMatrixTranspose * ( rotationToBodyFixedFramePartials.at( j ) * inertialRelativePosition );
            addTidalCorrectionTimeDerivatives(
                    [ & ]( const int degree, const int order ) {
                        return calculateTidalForcingTimeDerivative(
                                degree, order, currentMassRatio, sphericalPosition, sphericalPositionPartial );
                    },
                    correctionPartials.at( j ).first,
                    correctionPartials.at( j ).second,
                    maximumDegree );
        }
    }
    return correctionPartials;
}

//! Function for calculating partials of time derivatives of coefficient corrections w.r.t. orientation of deformed body
std::vector< std::pair< Eigen::MatrixXd, Eigen::MatrixXd > >
SolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrectionsTimeDerivativePartialWrtOrientation(
        const double time,
        const std::vector< Eigen::Matrix3d >& rotationToBodyFixedFramePartials,
        const std::function< Eigen::Matrix3d( const int, const int, const bool, const Eigen::Vector3d& ) >& forcingHessianFunction,
        const int maximumDegree )
{
    std::vector< std::pair< Eigen::MatrixXd, Eigen::MatrixXd > > correctionRatePartials(
            rotationToBodyFixedFramePartials.size( ),
            std::make_pair( Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ ),
                            Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ ) ) );
    if( rotationDerivativeFunction_ == nullptr )
    {
        return correctionRatePartials;
    }

    // Partial of rotation matrix derivative w.r.t. orientation, for constant angular velocity: dR/dt = -[omega x] R = ( dR/dt R^T ) R
    const Eigen::Matrix3d rotationToBodyFixedFrame = deformedBodyOrientationFunction_( time ).toRotationMatrix( );
    const Eigen::Matrix3d angularVelocityMatrix = rotationDerivativeFunction_( time ) * rotationToBodyFixedFrame.transpose( );

    for( unsigned int i = 0; i < deformingBodyStateFunctions_.size( ); i++ )
    {
        const Eigen::Vector6d bodyFixedRelativeState = getBodyFixedRelativeState( i, time );
        const Eigen::Vector3d bodyFixedPosition = bodyFixedRelativeState.segment( 0, 3 );
        const Eigen::Vector6d inertialRelativeState = deformingBodyStateFunctions_[ i ]( time ) - deformedBodyStateFunction_( time );
        const double currentMassRatio = deformingBodyMasses_[ i ]( ) / deformedBodyMass_( );
        const Eigen::Vector3d sphericalPosition = updateSphericalHarmonicsCache( bodyFixedPosition );
        const Eigen::Matrix3d gradientMatrixTranspose =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( bodyFixedPosition ).transpose( );

        for( unsigned int j = 0; j < rotationToBodyFixedFramePartials.size( ); j++ )
        {
            // Partials of body-fixed position and velocity of body causing deformation w.r.t. orientation
            const Eigen::Vector3d positionPartial = rotationToBodyFixedFramePartials.at( j ) * inertialRelativeState.segment( 0, 3 );
            const Eigen::Vector3d sphericalVelocityPartial = gradientMatrixTranspose *
                    ( rotationToBodyFixedFramePartials.at( j ) * inertialRelativeState.segment( 3, 3 ) +
                      angularVelocityMatrix * positionPartial );

            // Time derivative of forcing is linear in velocity, and depends on position through the gradient of the forcing
            addTidalCorrectionTimeDerivatives(
                    [ & ]( const int degree, const int order ) {
                        return calculateTidalForcingTimeDerivative(
                                       degree, order, currentMassRatio, sphericalPosition, sphericalVelocityPartial ) +
                                calculateTidalForcingTimeDerivativePartialWrtPosition(
                                        degree, order, currentMassRatio, bodyFixedRelativeState, positionPartial, forcingHessianFunction );
                    },
                    correctionRatePartials.at( j ).first,
                    correctionRatePartials.at( j ).second,
                    maximumDegree );
        }
    }
    return correctionRatePartials;
}

//! Function for calculating partials of coefficient corrections w.r.t. gravitational parameter of a body causing deformation
std::pair< Eigen::MatrixXd, Eigen::MatrixXd >
SolidBodyTideGravityFieldVariations::calculateSphericalHarmonicsCorrectionsPartialWrtDeformingBodyGravitationalParameter(
        const double time,
        const int bodyIndex )
{
    Eigen::MatrixXd cosineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineCorrections = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd cosineCorrectionsWithoutForcing = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    Eigen::MatrixXd sineCorrectionsWithoutForcing = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );

    // Set geometry of deformed body (index 0) and body causing deformation
    setBodyGeometryParameters( 0, time );
    setBodyGeometryParameters( bodyIndex, time );

    // The corrections are linear in the mass ratio, apart from (constant) mean forcing terms, which are removed by evaluating the
    // corrections for zero mass ratio
    const double gravitationalParameter = deformingBodyMasses_[ bodyIndex ]( );
    massRatio = gravitationalParameter / deformedBodyMass_( );
    for( unsigned int j = 0; j < correctionFunctions.size( ); j++ )
    {
        correctionFunctions[ j ]( cosineCorrections, sineCorrections );
    }
    massRatio = 0.0;
    for( unsigned int j = 0; j < correctionFunctions.size( ); j++ )
    {
        correctionFunctions[ j ]( cosineCorrectionsWithoutForcing, sineCorrectionsWithoutForcing );
    }

    return std::make_pair( ( cosineCorrections - cosineCorrectionsWithoutForcing ) / gravitationalParameter,
                           ( sineCorrections - sineCorrectionsWithoutForcing ) / gravitationalParameter );
}

//! Function to compute the time derivative of the tidal forcing of a single degree, order and body causing deformation
std::complex< double > SolidBodyTideGravityFieldVariations::calculateUnitLoveNumberTidalForcingTimeDerivative( const int degree,
                                                                                                               const int order,
                                                                                                               const int bodyIndex,
                                                                                                               const double time )
{
    if( rotationDerivativeFunction_ == nullptr )
    {
        throw std::runtime_error(
                "Error when computing time derivative of tidal forcing, no rotation matrix derivative function is defined" );
    }
    const Eigen::Vector6d bodyFixedRelativeState = getBodyFixedRelativeState( bodyIndex, time );
    const Eigen::Vector3d bodyFixedPosition = bodyFixedRelativeState.segment( 0, 3 );
    const Eigen::Vector3d sphericalPosition = updateSphericalHarmonicsCache( bodyFixedPosition );
    return calculateTidalForcingTimeDerivative(
            degree,
            order,
            deformingBodyMasses_[ bodyIndex ]( ) / deformedBodyMass_( ),
            sphericalPosition,
            coordinate_conversions::getSphericalToCartesianGradientMatrix( bodyFixedPosition ).transpose( ) *
                    bodyFixedRelativeState.segment( 3, 3 ) );
}

//! Function to add time derivatives of corrections due to the current body causing deformation.
void BasicSolidBodyTideGravityFieldVariations::addTidalCorrectionTimeDerivatives(
        const std::function< std::complex< double >( const int, const int ) >& unitLoveNumberForcingRateFunction,
        Eigen::MatrixXd& cosineRates,
        Eigen::MatrixXd& sineRates,
        const int maximumDegree )
{
    // Mean forcing terms are constant, and do not contribute
    for( const auto& loveNumberIt : loveNumbers_ )
    {
        int n = loveNumberIt.first;
        if( n > maximumDegree )
        {
            continue;
        }
        for( int m = 0; ( m <= n && m < static_cast< int >( loveNumberIt.second.size( ) ) ); m++ )
        {
            std::complex< double > stokesCoefficientCorrectionRate = loveNumberIt.second[ m ] * unitLoveNumberForcingRateFunction( n, m );
            cosineRates( n - minimumDegree_, m ) += stokesCoefficientCorrectionRate.real( );
            if( m != 0 )
            {
                sineRates( n - minimumDegree_, m ) -= stokesCoefficientCorrectionRate.imag( );
            }
        }
    }
}

//! Function to add time derivatives of corrections due to the current body causing deformation.
void ModeCoupledSolidBodyTideGravityFieldVariations::addTidalCorrectionTimeDerivatives(
        const std::function< std::complex< double >( const int, const int ) >& unitLoveNumberForcingRateFunction,
        Eigen::MatrixXd& cosineRates,
        Eigen::MatrixXd& sineRates,
        const int maximumDegree )
{
    for( const auto& loveNumberIt : loveNumbers_ )
    {
        // Only compute forcing if it contributes to a response at or below the maximum degree
        std::complex< double > unityLoveNumberForcingRate;
        bool isForcingRateComputed = false;
        for( const auto& responseIt : loveNumberIt.second )
        {
            int nResponse = responseIt.first.first;
            int mResponse = responseIt.first.second;
            if( nResponse > maximumDegree )
            {
                continue;
            }
            if( !isForcingRateComputed )
            {
                unityLoveNumberForcingRate = unitLoveNumberForcingRateFunction( loveNumberIt.first.first, loveNumberIt.first.second );
                isForcingRateComputed = true;
            }
            std::complex< double > stokesCoefficientCorrectionRate = responseIt.second * unityLoveNumberForcingRate;
            cosineRates( nResponse - minimumDegree_, mResponse ) += stokesCoefficientCorrectionRate.real( );
            if( mResponse != 0 )
            {
                sineRates( nResponse - minimumDegree_, mResponse ) -= stokesCoefficientCorrectionRate.imag( );
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
        deformedBodyPosition = std::move( deformedBodyStateFunction_( evaluationTime ) ).segment( 0, 3 );
        toDeformedBodyFrameRotation = deformedBodyOrientationFunction_( evaluationTime );
    }

    // Calculate current state of body causing deformation.
    Eigen::Vector3d relativeDeformingBodyPosition = toDeformedBodyFrameRotation *
            ( std::move( deformingBodyStateFunctions_[ bodyIndex ]( evaluationTime ) ).segment( 0, 3 ) - deformedBodyPosition );
    Eigen::Vector3d relativeDeformingBodySphericalPosition =
            coordinate_conversions::convertCartesianToSpherical( relativeDeformingBodyPosition );

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
void BasicSolidBodyTideGravityFieldVariations::addBasicSolidBodyTideCorrections( Eigen::MatrixXd& cTermCorrections,
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
                    calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude( loveNumbers_[ n ][ m ],
                                                                                       massRatio,
                                                                                       radiusRatioPower,
                                                                                       tideAmplitude,
                                                                                       tideArgument,
                                                                                       n,
                                                                                       m,
                                                                                       meanForcingCosineTerms_[ n ][ m ],
                                                                                       meanForcingSineTerms_[ n ][ m ] );

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
                        std::complex< double >( 1.0, 0.0 ), massRatio, radiusRatioPower, tideAmplitude, tideArgument, n, m );

        // Compute response at required degrees/orders with required love numbers
        for( auto responseIt : loveNumberIt.second )
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
