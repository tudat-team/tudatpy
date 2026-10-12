/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/math/basic/ellipsoidalHarmonics.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace basic_mathematics
{

namespace
{

//! Evaluates a polynomial whose coefficients are given in ascending powers.
double evaluateAscending( const Eigen::VectorXd& coefficients, const double x )
{
    double value = 0.0;
    for( int i = coefficients.size( ) - 1; i >= 0; i-- )
    {
        value = value * x + coefficients( i );
    }
    return value;
}

//! Evaluates the same coefficients read in descending powers, which gives the reciprocal
//! polynomial. The second-kind Lame integrand needs this form.
double evaluateDescending( const Eigen::VectorXd& coefficients, const double x )
{
    double value = 0.0;
    for( int i = 0; i < coefficients.size( ); i++ )
    {
        value = value * x + coefficients( i );
    }
    return value;
}

}  // namespace

int getNumberOfTermsInGroup( const int degree, const EllipsoidalHarmonicTermGroup group )
{
    switch( group )
    {
        case ellipsoidal_harmonic_k_terms:
            return getNumberOfKTerms( degree );
        case ellipsoidal_harmonic_l_terms:
            return getNumberOfLTerms( degree );
        case ellipsoidal_harmonic_m_terms:
            return getNumberOfMTerms( degree );
        case ellipsoidal_harmonic_n_terms:
            return getNumberOfNTerms( degree );
        default:
            throw std::runtime_error( "Error when counting ellipsoidal harmonic terms, group not recognised." );
    }
}

void getTermGroupAndIndex( const int degree, const int fileOrderIndex, EllipsoidalHarmonicTermGroup& group, int& indexInGroup )
{
    if( fileOrderIndex < 1 || fileOrderIndex > 2 * degree + 1 )
    {
        throw std::runtime_error( "Error when identifying ellipsoidal harmonic term, order index " + std::to_string( fileOrderIndex ) +
                                  " is out of range for degree " + std::to_string( degree ) + "." );
    }

    const int halfDegree = degree / 2;

    if( fileOrderIndex <= halfDegree + 1 )
    {
        group = ellipsoidal_harmonic_k_terms;
        indexInGroup = fileOrderIndex - 1;
    }
    else if( fileOrderIndex <= degree + 1 )
    {
        group = ellipsoidal_harmonic_l_terms;
        indexInGroup = fileOrderIndex - halfDegree - 2;
    }
    else if( fileOrderIndex <= 2 * degree + 1 - halfDegree )
    {
        group = ellipsoidal_harmonic_m_terms;
        indexInGroup = fileOrderIndex - degree - 2;
    }
    else
    {
        group = ellipsoidal_harmonic_n_terms;
        indexInGroup = fileOrderIndex - ( 2 * degree + 1 - halfDegree ) - 1;
    }
}

EllipsoidalHarmonicCoefficients createEllipsoidalHarmonicCoefficients( const int maximumDegree )
{
    EllipsoidalHarmonicCoefficients coefficients;
    coefficients.kTerms.resize( maximumDegree + 1 );
    coefficients.lTerms.resize( maximumDegree + 1 );
    coefficients.mTerms.resize( maximumDegree + 1 );
    coefficients.nTerms.resize( maximumDegree + 1 );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        coefficients.kTerms[ degree ] = Eigen::VectorXd::Zero( getNumberOfKTerms( degree ) );
        coefficients.lTerms[ degree ] = Eigen::VectorXd::Zero( getNumberOfLTerms( degree ) );
        coefficients.mTerms[ degree ] = Eigen::VectorXd::Zero( getNumberOfMTerms( degree ) );
        coefficients.nTerms[ degree ] = Eigen::VectorXd::Zero( getNumberOfNTerms( degree ) );
    }

    return coefficients;
}

Eigen::Vector3d convertCartesianToEllipsoidalCoordinates( const Eigen::Vector3d& bodyFixedPosition, const double h, const double k )
{
    const double xSquared = bodyFixedPosition.x( ) * bodyFixedPosition.x( );
    const double ySquared = bodyFixedPosition.y( ) * bodyFixedPosition.y( );
    const double zSquared = bodyFixedPosition.z( ) * bodyFixedPosition.z( );
    const double hSquared = h * h;
    const double kSquared = k * k;

    // Coefficients of the cubic in rho^2 whose largest root is the ellipsoidal radius. It comes
    // from requiring x^2/rho^2 + y^2/(rho^2-h^2) + z^2/(rho^2-k^2) = 1.
    const double alpha = hSquared + kSquared;
    const double beta = hSquared * kSquared;
    const double secondCoefficient = -( xSquared + ySquared + zSquared + alpha );
    const double firstCoefficient = xSquared * alpha + ySquared * kSquared + zSquared * hSquared + beta;
    const double zerothCoefficient = -xSquared * beta;

    // Trigonometric solution, which is the stable branch when all three roots are real.
    double p = ( secondCoefficient * secondCoefficient - 3.0 * firstCoefficient ) / 9.0;
    const double q = ( 9.0 * secondCoefficient * firstCoefficient - 27.0 * zerothCoefficient -
                       2.0 * secondCoefficient * secondCoefficient * secondCoefficient ) /
            54.0;

    double cosineArgument = q / std::sqrt( p * p * p );
    cosineArgument = std::max( -1.0, std::min( 1.0, cosineArgument ) );
    const double theta = std::acos( cosineArgument ) / 3.0;

    const double shift = secondCoefficient / 3.0;
    p = 2.0 * std::sqrt( p );

    const double rho = std::sqrt( p * std::cos( theta ) - shift );

    const double gamma = std::atan2( bodyFixedPosition.z( ) / std::sqrt( rho * rho - kSquared ),
                                     bodyFixedPosition.y( ) / std::sqrt( rho * rho - hSquared ) );

    double sineOfPhi = bodyFixedPosition.x( ) / rho;
    sineOfPhi = std::max( -1.0, std::min( 1.0, sineOfPhi ) );
    const double phi = std::asin( sineOfPhi );

    return ( Eigen::Vector3d( ) << rho, phi, gamma ).finished( );
}

Eigen::Vector3d convertEllipsoidalToCartesianCoordinates( const Eigen::Vector3d& ellipsoidalCoordinates, const double h, const double k )
{
    const double rho = ellipsoidalCoordinates( 0 );
    const double phi = ellipsoidalCoordinates( 1 );
    const double gamma = ellipsoidalCoordinates( 2 );

    return ( Eigen::Vector3d( ) << rho * std::sin( phi ),
             std::sqrt( rho * rho - h * h ) * std::cos( phi ) * std::cos( gamma ),
             std::sqrt( rho * rho - k * k ) * std::cos( phi ) * std::sin( gamma ) )
            .finished( );
}

void computeGaussLegendreNodesAndWeights( const int order, Eigen::VectorXd& nodes, Eigen::VectorXd& weights )
{
    if( order < 1 )
    {
        throw std::runtime_error( "Error when creating Gauss-Legendre rule, order must be at least 1." );
    }

    nodes = Eigen::VectorXd::Zero( order );
    weights = Eigen::VectorXd::Zero( order );

    // The nodes are symmetric about zero, so only half of them are solved for.
    for( int i = 0; i < ( order + 1 ) / 2; i++ )
    {
        double root = std::cos( mathematical_constants::PI * ( i + 0.75 ) / ( order + 0.5 ) );
        double derivative = 1.0;

        for( int iteration = 0; iteration < 100; iteration++ )
        {
            // Legendre polynomial of the given order, by its three-term recurrence.
            double previous = 1.0;
            double current = root;
            for( int degree = 2; degree <= order; degree++ )
            {
                const double next = ( ( 2 * degree - 1 ) * root * current - ( degree - 1 ) * previous ) / degree;
                previous = current;
                current = next;
            }

            if( order == 1 )
            {
                current = root;
                previous = 1.0;
            }

            derivative = order * ( root * current - previous ) / ( root * root - 1.0 );
            const double step = -current / derivative;
            root += step;

            if( std::fabs( step ) < 1.0E-15 )
            {
                break;
            }
        }

        const double weight = 2.0 / ( ( 1.0 - root * root ) * derivative * derivative );

        nodes( i ) = -root;
        nodes( order - 1 - i ) = root;
        weights( i ) = weight;
        weights( order - 1 - i ) = weight;
    }
}

double evaluateLameFunctionFirstKind( const Eigen::VectorXd& polynomialCoefficients,
                                      const int degree,
                                      const EllipsoidalHarmonicTermGroup group,
                                      const double rho,
                                      const double h,
                                      const double k )
{
    const double scaledRadiusSquared = ( rho / h ) * ( rho / h );
    const double focalRatioSquared = ( k * k ) / ( h * h );
    const double polynomial = evaluateAscending( polynomialCoefficients, scaledRadiusSquared );
    const bool degreeIsOdd = ( degree % 2 != 0 );

    double value;
    bool multiplyByScaledRadius;

    switch( group )
    {
        case ellipsoidal_harmonic_k_terms:
            value = polynomial;
            multiplyByScaledRadius = degreeIsOdd;
            break;
        case ellipsoidal_harmonic_l_terms:
            value = std::sqrt( scaledRadiusSquared - 1.0 ) * polynomial;
            multiplyByScaledRadius = !degreeIsOdd;
            break;
        case ellipsoidal_harmonic_m_terms:
            value = std::sqrt( scaledRadiusSquared - focalRatioSquared ) * polynomial;
            multiplyByScaledRadius = !degreeIsOdd;
            break;
        default:
            value = std::sqrt( scaledRadiusSquared - focalRatioSquared ) * std::sqrt( scaledRadiusSquared - 1.0 ) * polynomial;
            multiplyByScaledRadius = degreeIsOdd;
            break;
    }

    if( multiplyByScaledRadius )
    {
        value *= rho / h;
    }

    return value;
}

double evaluateLameFunctionSecondKind( const Eigen::VectorXd& polynomialCoefficients,
                                       const int degree,
                                       const EllipsoidalHarmonicTermGroup group,
                                       const double rho,
                                       const double h,
                                       const double k,
                                       const Eigen::VectorXd& quadratureNodes,
                                       const Eigen::VectorXd& quadratureWeights )
{
    if( rho <= k )
    {
        throw std::runtime_error( "Error when evaluating ellipsoidal harmonics, the ellipsoidal radius " + std::to_string( rho ) +
                                  " is not larger than the focal parameter " + std::to_string( k ) +
                                  "; the point lies inside the reference ellipsoid." );
    }

    const double focalRatioSquared = ( k * k ) / ( h * h );
    const double upperLimit = h / rho;

    // Integral over s from 0 to h/rho. The integrand is smooth here because rho > k keeps both
    // square roots away from zero, so a fixed-order rule converges quickly.
    double integral = 0.0;
    for( int i = 0; i < quadratureNodes.size( ); i++ )
    {
        const double s = 0.5 * upperLimit * ( quadratureNodes( i ) + 1.0 );
        const double sSquared = s * s;
        const double scaledSquared = sSquared * focalRatioSquared;
        const double polynomial = evaluateDescending( polynomialCoefficients, sSquared );

        double lameFactor;
        switch( group )
        {
            case ellipsoidal_harmonic_k_terms:
                lameFactor = polynomial;
                break;
            case ellipsoidal_harmonic_l_terms:
                lameFactor = std::sqrt( 1.0 - sSquared ) * polynomial;
                break;
            case ellipsoidal_harmonic_m_terms:
                lameFactor = std::sqrt( 1.0 - scaledSquared ) * polynomial;
                break;
            default:
                lameFactor = std::sqrt( 1.0 - sSquared ) * std::sqrt( 1.0 - scaledSquared ) * polynomial;
                break;
        }

        const double integrand =
                std::pow( sSquared, degree ) / ( h * lameFactor * lameFactor * std::sqrt( ( 1.0 - scaledSquared ) * ( 1.0 - sSquared ) ) );

        integral += 0.5 * upperLimit * quadratureWeights( i ) * integrand;
    }

    return ( 2 * degree + 1 ) * evaluateLameFunctionFirstKind( polynomialCoefficients, degree, group, rho, h, k ) * integral;
}

void computeRadialFunctions( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                             const double rho,
                             const double h,
                             const double k,
                             const Eigen::VectorXd& quadratureNodes,
                             const Eigen::VectorXd& quadratureWeights,
                             EllipsoidalHarmonicCoefficients& radialFunctions )
{
    const int maximumDegree = lameFunctionCoefficients.getMaximumDegree( );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const Eigen::MatrixXd& polynomials = lameFunctionCoefficients.getGroup( degree, group );
            Eigen::VectorXd& values = radialFunctions.getGroup( degree, group );

            for( int term = 0; term < getNumberOfTermsInGroup( degree, group ); term++ )
            {
                values( term ) = evaluateLameFunctionSecondKind(
                        polynomials.col( term ), degree, group, rho, h, k, quadratureNodes, quadratureWeights );
            }
        }
    }
}

EllipsoidalHarmonicCoefficients normaliseCoefficientsToReferenceSurface(
        const EllipsoidalHarmonicCoefficients& fileCoefficients,
        const EllipsoidalHarmonicCoefficients& radialFunctionsAtReferenceSurface )
{
    if( fileCoefficients.getMaximumDegree( ) != radialFunctionsAtReferenceSurface.getMaximumDegree( ) )
    {
        throw std::runtime_error(
                "Error when normalising ellipsoidal harmonic coefficients, the coefficients and "
                "the reference radial functions have different maximum degrees." );
    }

    EllipsoidalHarmonicCoefficients normalised = fileCoefficients;

    for( int degree = 0; degree <= normalised.getMaximumDegree( ); degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            normalised.getGroup( degree, group ) =
                    normalised.getGroup( degree, group ).cwiseQuotient( radialFunctionsAtReferenceSurface.getGroup( degree, group ) );
        }
    }

    return normalised;
}

double sumEllipsoidalHarmonicTerms( const EllipsoidalHarmonicCoefficients& coefficients,
                                    const EllipsoidalHarmonicCoefficients& radialFunctions,
                                    const EllipsoidalHarmonicMatrices& transformFactors,
                                    const LegendreCache& legendreCache,
                                    const double gamma )
{
    const int maximumDegree = coefficients.getMaximumDegree( );
    double potential = 0.0;

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        // Each group pairs its terms with a different sequence of Legendre orders and with
        // either a cosine or a sine of a multiple of gamma.
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const int numberOfTerms = getNumberOfTermsInGroup( degree, group );

            if( numberOfTerms == 0 )
            {
                continue;
            }

            Eigen::VectorXd surfaceTerms = Eigen::VectorXd::Zero( numberOfTerms );

            for( int term = 0; term < numberOfTerms; term++ )
            {
                int legendreOrder;
                double angleFactor;

                switch( group )
                {
                    case ellipsoidal_harmonic_k_terms:
                        legendreOrder = 2 * term;
                        angleFactor = std::cos( 2.0 * term * gamma );
                        break;
                    case ellipsoidal_harmonic_l_terms:
                        legendreOrder = 2 * ( term + 1 ) - 1;
                        angleFactor = std::cos( ( 2.0 * ( term + 1 ) - 1.0 ) * gamma );
                        break;
                    case ellipsoidal_harmonic_m_terms:
                        legendreOrder = 2 * ( term + 1 ) - 1;
                        angleFactor = std::sin( ( 2.0 * ( term + 1 ) - 1.0 ) * gamma );
                        break;
                    default:
                        legendreOrder = 2 * ( term + 1 );
                        angleFactor = std::sin( 2.0 * ( term + 1 ) * gamma );
                        break;
                }

                surfaceTerms( term ) = legendreCache.getLegendrePolynomial( degree, legendreOrder ) * angleFactor;
            }

            const Eigen::VectorXd weightedCoefficients =
                    coefficients.getGroup( degree, group ).cwiseProduct( radialFunctions.getGroup( degree, group ) );

            potential += weightedCoefficients.dot( transformFactors.getGroup( degree, group ) * surfaceTerms );
        }
    }

    return potential;
}

double computeEllipsoidalHarmonicGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition,
                                                         const double h,
                                                         const double k,
                                                         const EllipsoidalHarmonicCoefficients& coefficients,
                                                         const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                                         const EllipsoidalHarmonicMatrices& transformFactors,
                                                         const int quadratureOrder )
{
    const int maximumDegree = coefficients.getMaximumDegree( );

    const Eigen::Vector3d ellipsoidalCoordinates = convertCartesianToEllipsoidalCoordinates( bodyFixedPosition, h, k );

    Eigen::VectorXd quadratureNodes;
    Eigen::VectorXd quadratureWeights;
    computeGaussLegendreNodesAndWeights( quadratureOrder, quadratureNodes, quadratureWeights );

    EllipsoidalHarmonicCoefficients radialFunctions = createEllipsoidalHarmonicCoefficients( maximumDegree );
    computeRadialFunctions(
            lameFunctionCoefficients, ellipsoidalCoordinates( 0 ), h, k, quadratureNodes, quadratureWeights, radialFunctions );

    LegendreCache legendreCache( maximumDegree, maximumDegree, true );

    // Only the polynomial values are needed. Their derivatives are singular where the long
    // axis of the body crosses the surface of constant ellipsoidal radius, and asking for them
    // there is an error, so they are switched off.
    legendreCache.setComputeFirstDerivatives( false );
    legendreCache.setComputeSecondDerivatives( false );
    legendreCache.update( std::sin( ellipsoidalCoordinates( 1 ) ) );

    return sumEllipsoidalHarmonicTerms( coefficients, radialFunctions, transformFactors, legendreCache, ellipsoidalCoordinates( 2 ) );
}

}  // namespace basic_mathematics

}  // namespace tudat
