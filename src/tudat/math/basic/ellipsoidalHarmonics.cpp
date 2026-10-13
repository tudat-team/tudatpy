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
#include <vector>

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

//! Derivative of a polynomial whose coefficients are given in ascending order of power.
double evaluateAscendingDerivative( const Eigen::VectorXd& coefficients, const double x )
{
    double value = 0.0;
    for( int i = coefficients.size( ) - 1; i >= 1; i-- )
    {
        value = value * x + i * coefficients( i );
    }
    return value;
}

/*! Everything about the second-kind integral that does not depend on which term is evaluated.
 *
 * The integral runs over s from 0 to h/rho, and all of its terms share the same nodes. Most of
 * the integrand is the same for every term too: only the Lame polynomial and the degree differ.
 * Those shared parts are gathered here once per ellipsoidal radius, which takes the square
 * roots and the powers of s out of the per-term loop. For a degree 15 field that is 256 terms
 * reading one table instead of rebuilding it 256 times.
 *
 * The last entry is the upper limit of the integral rather than a quadrature node. Its weight is
 * zero so it does not contribute to the sum, and the integrand there is what the derivative with
 * respect to the ellipsoidal radius needs.
 */
struct SecondKindIntegralNodes {
    //! Square of the integration variable at each node.
    Eigen::VectorXd sSquared;

    //! One minus that, which is the squared Lame factor of the L group.
    Eigen::VectorXd oneMinusSSquared;

    //! One minus the focal-ratio-scaled square, the squared Lame factor of the M group.
    Eigen::VectorXd oneMinusScaledSquared;

    //! Reciprocal of the square root in the denominator of the integrand.
    Eigen::VectorXd inverseRootOfProduct;

    //! Quadrature weight already scaled by the width of the interval.
    Eigen::VectorXd scaledWeight;

    //! Entry ( i, n ) is sSquared( i ) to the power n, built by repeated multiplication.
    Eigen::MatrixXd powersOfSSquared;

    //! Index of the entry that holds the upper limit rather than a quadrature node.
    int upperLimitIndex;
};

/*! Gathers the shared part of the second-kind integral at one ellipsoidal radius.
 *
 * @param rho Ellipsoidal radius, which must be larger than k.
 * @param h First focal parameter of the reference ellipsoid.
 * @param k Second focal parameter of the reference ellipsoid.
 * @param quadratureNodes Gauss-Legendre nodes on [-1, 1].
 * @param quadratureWeights Matching Gauss-Legendre weights.
 * @param maximumDegree Highest degree that will be asked for.
 * @param nodes Output.
 */
void buildSecondKindIntegralNodes( const double rho,
                                   const double h,
                                   const double k,
                                   const Eigen::VectorXd& quadratureNodes,
                                   const Eigen::VectorXd& quadratureWeights,
                                   const int maximumDegree,
                                   SecondKindIntegralNodes& nodes )
{
    if( rho <= k )
    {
        throw std::runtime_error( "Error when evaluating ellipsoidal harmonics, the ellipsoidal radius " + std::to_string( rho ) +
                                  " is not larger than the focal parameter " + std::to_string( k ) +
                                  "; the point lies inside the reference ellipsoid." );
    }

    const double focalRatioSquared = ( k * k ) / ( h * h );
    const double upperLimit = h / rho;
    const int numberOfQuadratureNodes = static_cast< int >( quadratureNodes.size( ) );
    const int numberOfEntries = numberOfQuadratureNodes + 1;

    nodes.upperLimitIndex = numberOfQuadratureNodes;
    nodes.sSquared.resize( numberOfEntries );
    nodes.oneMinusSSquared.resize( numberOfEntries );
    nodes.oneMinusScaledSquared.resize( numberOfEntries );
    nodes.inverseRootOfProduct.resize( numberOfEntries );
    nodes.scaledWeight.resize( numberOfEntries );
    nodes.powersOfSSquared.resize( numberOfEntries, maximumDegree + 1 );

    for( int i = 0; i < numberOfEntries; i++ )
    {
        // The last entry is the upper limit, and carries no weight.
        const bool isUpperLimit = ( i == nodes.upperLimitIndex );
        const double s = isUpperLimit ? upperLimit : 0.5 * upperLimit * ( quadratureNodes( i ) + 1.0 );

        const double sSquared = s * s;
        nodes.sSquared( i ) = sSquared;
        nodes.oneMinusSSquared( i ) = 1.0 - sSquared;
        nodes.oneMinusScaledSquared( i ) = 1.0 - sSquared * focalRatioSquared;
        nodes.inverseRootOfProduct( i ) = 1.0 / std::sqrt( nodes.oneMinusScaledSquared( i ) * nodes.oneMinusSSquared( i ) );
        nodes.scaledWeight( i ) = isUpperLimit ? 0.0 : 0.5 * upperLimit * quadratureWeights( i );

        // Powers of s^2 by repeated multiplication, in place of a library power call per term.
        nodes.powersOfSSquared( i, 0 ) = 1.0;
        for( int degree = 1; degree <= maximumDegree; degree++ )
        {
            nodes.powersOfSSquared( i, degree ) = nodes.powersOfSSquared( i, degree - 1 ) * sSquared;
        }
    }
}

/*! The part of the squared Lame factor that the group contributes, from the shared table.
 *
 * The integrand divides by the square of the Lame factor, so the square roots that the groups
 * put in front of their polynomials can be taken squared and never evaluated as roots at all.
 */
double getSquaredGroupFactor( const EllipsoidalHarmonicTermGroup group, const SecondKindIntegralNodes& nodes, const int entry )
{
    switch( group )
    {
        case ellipsoidal_harmonic_k_terms:
            return 1.0;
        case ellipsoidal_harmonic_l_terms:
            return nodes.oneMinusSSquared( entry );
        case ellipsoidal_harmonic_m_terms:
            return nodes.oneMinusScaledSquared( entry );
        default:
            return nodes.oneMinusSSquared( entry ) * nodes.oneMinusScaledSquared( entry );
    }
}

/*! The second-kind integral of one term, and the integrand at the upper limit.
 *
 * @param polynomialCoefficients Lame polynomial of the term, in ascending powers.
 * @param degree Degree of the term.
 * @param group Group the term belongs to.
 * @param h First focal parameter of the reference ellipsoid.
 * @param nodes Shared table from buildSecondKindIntegralNodes.
 * @param integral Output.
 * @param integrandAtUpperLimit Output, needed only for the derivative.
 */
void integrateSecondKindTerm( const Eigen::VectorXd& polynomialCoefficients,
                              const int degree,
                              const EllipsoidalHarmonicTermGroup group,
                              const double h,
                              const SecondKindIntegralNodes& nodes,
                              double& integral,
                              double& integrandAtUpperLimit )
{
    integral = 0.0;

    for( int i = 0; i <= nodes.upperLimitIndex; i++ )
    {
        const double polynomial = evaluateDescending( polynomialCoefficients, nodes.sSquared( i ) );
        const double squaredLameFactor = getSquaredGroupFactor( group, nodes, i ) * polynomial * polynomial;

        const double integrand = nodes.powersOfSSquared( i, degree ) * nodes.inverseRootOfProduct( i ) / ( h * squaredLameFactor );

        if( i == nodes.upperLimitIndex )
        {
            integrandAtUpperLimit = integrand;
        }
        else
        {
            integral += nodes.scaledWeight( i ) * integrand;
        }
    }
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
    // Builds the shared table for a single term, so this is the convenient route rather than the
    // fast one. Anything evaluating many terms at one radius should use computeRadialFunctions.
    SecondKindIntegralNodes nodes;
    buildSecondKindIntegralNodes( rho, h, k, quadratureNodes, quadratureWeights, degree, nodes );

    double integral;
    double integrandAtUpperLimit;
    integrateSecondKindTerm( polynomialCoefficients, degree, group, h, nodes, integral, integrandAtUpperLimit );

    return ( 2 * degree + 1 ) * evaluateLameFunctionFirstKind( polynomialCoefficients, degree, group, rho, h, k ) * integral;
}

double evaluateLameFunctionFirstKindDerivative( const Eigen::VectorXd& polynomialCoefficients,
                                                const int degree,
                                                const EllipsoidalHarmonicTermGroup group,
                                                const double rho,
                                                const double h,
                                                const double k )
{
    const double scaledRadiusSquared = ( rho / h ) * ( rho / h );
    const double focalRatioSquared = ( k * k ) / ( h * h );

    const double polynomial = evaluateAscending( polynomialCoefficients, scaledRadiusSquared );
    // The polynomial is in (rho/h)^2, so differentiating it brings out a factor 2*rho/h^2.
    const double polynomialDerivative = evaluateAscendingDerivative( polynomialCoefficients, scaledRadiusSquared ) * 2.0 * rho / ( h * h );

    // Square-root factor that the group puts in front of the polynomial, and its derivative.
    double rootFactor;
    double rootFactorDerivative;
    switch( group )
    {
        case ellipsoidal_harmonic_k_terms:
            rootFactor = 1.0;
            rootFactorDerivative = 0.0;
            break;
        case ellipsoidal_harmonic_l_terms:
            rootFactor = std::sqrt( scaledRadiusSquared - 1.0 );
            rootFactorDerivative = rho / ( h * h * rootFactor );
            break;
        case ellipsoidal_harmonic_m_terms:
            rootFactor = std::sqrt( scaledRadiusSquared - focalRatioSquared );
            rootFactorDerivative = rho / ( h * h * rootFactor );
            break;
        default: {
            const double firstRoot = std::sqrt( scaledRadiusSquared - focalRatioSquared );
            const double secondRoot = std::sqrt( scaledRadiusSquared - 1.0 );
            rootFactor = firstRoot * secondRoot;
            rootFactorDerivative = rho / ( h * h ) * ( secondRoot / firstRoot + firstRoot / secondRoot );
            break;
        }
    }

    double derivative = rootFactorDerivative * polynomial + rootFactor * polynomialDerivative;

    // One parity class of each group carries an extra factor rho/h.
    const bool degreeIsOdd = ( degree % 2 != 0 );
    const bool hasExtraRadiusFactor =
            ( group == ellipsoidal_harmonic_k_terms || group == ellipsoidal_harmonic_n_terms ) ? degreeIsOdd : !degreeIsOdd;
    if( hasExtraRadiusFactor )
    {
        derivative = derivative * ( rho / h ) + rootFactor * polynomial / h;
    }

    return derivative;
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

    // The nodes and everything derived from them are shared by every term at this radius.
    SecondKindIntegralNodes nodes;
    buildSecondKindIntegralNodes( rho, h, k, quadratureNodes, quadratureWeights, maximumDegree, nodes );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const Eigen::MatrixXd& polynomials = lameFunctionCoefficients.getGroup( degree, group );
            Eigen::VectorXd& values = radialFunctions.getGroup( degree, group );

            for( int term = 0; term < getNumberOfTermsInGroup( degree, group ); term++ )
            {
                double integral;
                double integrandAtUpperLimit;
                integrateSecondKindTerm( polynomials.col( term ), degree, group, h, nodes, integral, integrandAtUpperLimit );

                values( term ) =
                        ( 2 * degree + 1 ) * evaluateLameFunctionFirstKind( polynomials.col( term ), degree, group, rho, h, k ) * integral;
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

int getTermLegendreOrder( const EllipsoidalHarmonicTermGroup group, const int indexInGroup )
{
    switch( group )
    {
        case ellipsoidal_harmonic_k_terms:
            return 2 * indexInGroup;
        case ellipsoidal_harmonic_l_terms:
            return 2 * indexInGroup + 1;
        case ellipsoidal_harmonic_m_terms:
            return 2 * indexInGroup + 1;
        default:
            return 2 * indexInGroup + 2;
    }
}

bool termUsesSineOfAngle( const EllipsoidalHarmonicTermGroup group )
{
    return group == ellipsoidal_harmonic_m_terms || group == ellipsoidal_harmonic_n_terms;
}

void computeRadialFunctionsAndDerivatives( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                           const double rho,
                                           const double h,
                                           const double k,
                                           const Eigen::VectorXd& quadratureNodes,
                                           const Eigen::VectorXd& quadratureWeights,
                                           EllipsoidalHarmonicCoefficients& radialFunctions,
                                           EllipsoidalHarmonicCoefficients& radialFunctionDerivatives )
{
    const int maximumDegree = lameFunctionCoefficients.getMaximumDegree( );

    // The nodes and everything derived from them are shared by every term at this radius.
    SecondKindIntegralNodes nodes;
    buildSecondKindIntegralNodes( rho, h, k, quadratureNodes, quadratureWeights, maximumDegree, nodes );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const Eigen::MatrixXd& polynomials = lameFunctionCoefficients.getGroup( degree, group );
            Eigen::VectorXd& values = radialFunctions.getGroup( degree, group );
            Eigen::VectorXd& derivatives = radialFunctionDerivatives.getGroup( degree, group );

            for( int term = 0; term < getNumberOfTermsInGroup( degree, group ); term++ )
            {
                double integral;
                double integrandAtUpperLimit;
                integrateSecondKindTerm( polynomials.col( term ), degree, group, h, nodes, integral, integrandAtUpperLimit );

                const double firstKind = evaluateLameFunctionFirstKind( polynomials.col( term ), degree, group, rho, h, k );
                const double firstKindDerivative =
                        evaluateLameFunctionFirstKindDerivative( polynomials.col( term ), degree, group, rho, h, k );

                values( term ) = ( 2 * degree + 1 ) * firstKind * integral;

                // The integral runs up to h/rho, so its derivative is the integrand there times
                // the derivative of that limit. No second quadrature is needed.
                derivatives( term ) =
                        ( 2 * degree + 1 ) * ( firstKindDerivative * integral - firstKind * ( h / ( rho * rho ) ) * integrandAtUpperLimit );
            }
        }
    }
}

Eigen::Matrix3d computeEllipsoidalCoordinateJacobian( const Eigen::Vector3d& ellipsoidalCoordinates, const double h, const double k )
{
    const double rho = ellipsoidalCoordinates( 0 );
    const double sineOfPhi = std::sin( ellipsoidalCoordinates( 1 ) );
    const double cosineOfPhi = std::cos( ellipsoidalCoordinates( 1 ) );
    const double sineOfGamma = std::sin( ellipsoidalCoordinates( 2 ) );
    const double cosineOfGamma = std::cos( ellipsoidalCoordinates( 2 ) );

    const double firstRadius = std::sqrt( rho * rho - h * h );
    const double secondRadius = std::sqrt( rho * rho - k * k );

    Eigen::Matrix3d jacobian;
    jacobian << sineOfPhi, rho * cosineOfPhi, 0.0, ( rho / firstRadius ) * cosineOfPhi * cosineOfGamma,
            -firstRadius * sineOfPhi * cosineOfGamma, -firstRadius * cosineOfPhi * sineOfGamma,
            ( rho / secondRadius ) * cosineOfPhi * sineOfGamma, -secondRadius * sineOfPhi * sineOfGamma,
            secondRadius * cosineOfPhi * cosineOfGamma;

    return jacobian;
}

Eigen::Vector3d computeEllipsoidalHarmonicPotentialGradient( const EllipsoidalHarmonicCoefficients& coefficients,
                                                             const EllipsoidalHarmonicCoefficients& radialFunctions,
                                                             const EllipsoidalHarmonicCoefficients& radialFunctionDerivatives,
                                                             const EllipsoidalHarmonicMatrices& transformFactors,
                                                             const LegendreCache& legendreCache,
                                                             const Eigen::Vector3d& ellipsoidalCoordinates,
                                                             const double h,
                                                             const double k )
{
    const int maximumDegree = coefficients.getMaximumDegree( );
    const double cosineOfPhi = std::cos( ellipsoidalCoordinates( 1 ) );
    const double gamma = ellipsoidalCoordinates( 2 );

    // Derivatives of the potential with respect to the three ellipsoidal coordinates.
    Eigen::Vector3d gradientInEllipsoidalCoordinates = Eigen::Vector3d::Zero( );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const int numberOfTerms = getNumberOfTermsInGroup( degree, group );

            if( numberOfTerms == 0 )
            {
                continue;
            }

            Eigen::VectorXd surfaceTerms = Eigen::VectorXd::Zero( numberOfTerms );
            Eigen::VectorXd surfaceTermsWrtPhi = Eigen::VectorXd::Zero( numberOfTerms );
            Eigen::VectorXd surfaceTermsWrtGamma = Eigen::VectorXd::Zero( numberOfTerms );

            for( int term = 0; term < numberOfTerms; term++ )
            {
                const int legendreOrder = getTermLegendreOrder( group, term );
                const double angle = legendreOrder * gamma;

                const double legendreValue = legendreCache.getLegendrePolynomial( degree, legendreOrder );
                // The cache differentiates with respect to sin( phi ), so the chain rule adds a
                // factor cos( phi ) to get the derivative with respect to phi itself.
                const double legendreDerivative = legendreCache.getLegendrePolynomialDerivative( degree, legendreOrder ) * cosineOfPhi;

                if( termUsesSineOfAngle( group ) )
                {
                    surfaceTerms( term ) = legendreValue * std::sin( angle );
                    surfaceTermsWrtPhi( term ) = legendreDerivative * std::sin( angle );
                    surfaceTermsWrtGamma( term ) = legendreOrder * legendreValue * std::cos( angle );
                }
                else
                {
                    surfaceTerms( term ) = legendreValue * std::cos( angle );
                    surfaceTermsWrtPhi( term ) = legendreDerivative * std::cos( angle );
                    surfaceTermsWrtGamma( term ) = -legendreOrder * legendreValue * std::sin( angle );
                }
            }

            const Eigen::VectorXd& groupCoefficients = coefficients.getGroup( degree, group );
            const Eigen::MatrixXd& groupTransform = transformFactors.getGroup( degree, group );

            const Eigen::VectorXd weightedByValue = groupCoefficients.cwiseProduct( radialFunctions.getGroup( degree, group ) );
            const Eigen::VectorXd weightedByDerivative =
                    groupCoefficients.cwiseProduct( radialFunctionDerivatives.getGroup( degree, group ) );

            gradientInEllipsoidalCoordinates( 0 ) += weightedByDerivative.dot( groupTransform * surfaceTerms );
            gradientInEllipsoidalCoordinates( 1 ) += weightedByValue.dot( groupTransform * surfaceTermsWrtPhi );
            gradientInEllipsoidalCoordinates( 2 ) += weightedByValue.dot( groupTransform * surfaceTermsWrtGamma );
        }
    }

    // Convert to Cartesian. With J the derivatives of position with respect to the ellipsoidal
    // coordinates, the Cartesian gradient solves J^T g = gradientInEllipsoidalCoordinates.
    const Eigen::Matrix3d jacobian = computeEllipsoidalCoordinateJacobian( ellipsoidalCoordinates, h, k );

    return jacobian.transpose( ).colPivHouseholderQr( ).solve( gradientInEllipsoidalCoordinates );
}

void computeRadialFunctionsAndSecondDerivatives( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                                 const double rho,
                                                 const double h,
                                                 const double k,
                                                 const Eigen::VectorXd& quadratureNodes,
                                                 const Eigen::VectorXd& quadratureWeights,
                                                 EllipsoidalHarmonicCoefficients& radialFunctions,
                                                 EllipsoidalHarmonicCoefficients& radialFunctionDerivatives,
                                                 EllipsoidalHarmonicCoefficients& radialFunctionSecondDerivatives )
{
    const int maximumDegree = lameFunctionCoefficients.getMaximumDegree( );

    computeRadialFunctionsAndDerivatives(
            lameFunctionCoefficients, rho, h, k, quadratureNodes, quadratureWeights, radialFunctions, radialFunctionDerivatives );

    // The first derivatives are analytic and depend on the ellipsoidal radius alone, so the
    // second ones are taken by differencing them in that one variable.
    const double step = rho * 1.0E-6;

    EllipsoidalHarmonicCoefficients valuesForward = createEllipsoidalHarmonicCoefficients( maximumDegree );
    EllipsoidalHarmonicCoefficients derivativesForward = createEllipsoidalHarmonicCoefficients( maximumDegree );
    EllipsoidalHarmonicCoefficients valuesBackward = createEllipsoidalHarmonicCoefficients( maximumDegree );
    EllipsoidalHarmonicCoefficients derivativesBackward = createEllipsoidalHarmonicCoefficients( maximumDegree );

    computeRadialFunctionsAndDerivatives(
            lameFunctionCoefficients, rho + step, h, k, quadratureNodes, quadratureWeights, valuesForward, derivativesForward );
    computeRadialFunctionsAndDerivatives(
            lameFunctionCoefficients, rho - step, h, k, quadratureNodes, quadratureWeights, valuesBackward, derivativesBackward );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            radialFunctionSecondDerivatives.getGroup( degree, group ) =
                    ( derivativesForward.getGroup( degree, group ) - derivativesBackward.getGroup( degree, group ) ) / ( 2.0 * step );
        }
    }
}

void computeEllipsoidalCoordinateSecondDerivatives( const Eigen::Vector3d& ellipsoidalCoordinates,
                                                    const double h,
                                                    const double k,
                                                    std::vector< Eigen::Matrix3d >& secondDerivatives )
{
    const double rho = ellipsoidalCoordinates( 0 );
    const double sineOfPhi = std::sin( ellipsoidalCoordinates( 1 ) );
    const double cosineOfPhi = std::cos( ellipsoidalCoordinates( 1 ) );
    const double sineOfGamma = std::sin( ellipsoidalCoordinates( 2 ) );
    const double cosineOfGamma = std::cos( ellipsoidalCoordinates( 2 ) );

    const double firstRadius = std::sqrt( rho * rho - h * h );
    const double secondRadius = std::sqrt( rho * rho - k * k );

    // First and second derivatives of the two radii with respect to the ellipsoidal radius.
    const double firstRadiusDerivative = rho / firstRadius;
    const double firstRadiusSecondDerivative = -( h * h ) / ( firstRadius * firstRadius * firstRadius );
    const double secondRadiusDerivative = rho / secondRadius;
    const double secondRadiusSecondDerivative = -( k * k ) / ( secondRadius * secondRadius * secondRadius );

    secondDerivatives.resize( 3 );

    // x = rho * sin( phi )
    Eigen::Matrix3d& forX = secondDerivatives[ 0 ];
    forX.setZero( );
    forX( 0, 1 ) = cosineOfPhi;
    forX( 1, 0 ) = cosineOfPhi;
    forX( 1, 1 ) = -rho * sineOfPhi;

    // y = sqrt( rho^2 - h^2 ) * cos( phi ) * cos( gamma )
    Eigen::Matrix3d& forY = secondDerivatives[ 1 ];
    forY( 0, 0 ) = firstRadiusSecondDerivative * cosineOfPhi * cosineOfGamma;
    forY( 0, 1 ) = -firstRadiusDerivative * sineOfPhi * cosineOfGamma;
    forY( 0, 2 ) = -firstRadiusDerivative * cosineOfPhi * sineOfGamma;
    forY( 1, 1 ) = -firstRadius * cosineOfPhi * cosineOfGamma;
    forY( 1, 2 ) = firstRadius * sineOfPhi * sineOfGamma;
    forY( 2, 2 ) = -firstRadius * cosineOfPhi * cosineOfGamma;
    forY( 1, 0 ) = forY( 0, 1 );
    forY( 2, 0 ) = forY( 0, 2 );
    forY( 2, 1 ) = forY( 1, 2 );

    // z = sqrt( rho^2 - k^2 ) * cos( phi ) * sin( gamma )
    Eigen::Matrix3d& forZ = secondDerivatives[ 2 ];
    forZ( 0, 0 ) = secondRadiusSecondDerivative * cosineOfPhi * sineOfGamma;
    forZ( 0, 1 ) = -secondRadiusDerivative * sineOfPhi * sineOfGamma;
    forZ( 0, 2 ) = secondRadiusDerivative * cosineOfPhi * cosineOfGamma;
    forZ( 1, 1 ) = -secondRadius * cosineOfPhi * sineOfGamma;
    forZ( 1, 2 ) = -secondRadius * sineOfPhi * cosineOfGamma;
    forZ( 2, 2 ) = -secondRadius * cosineOfPhi * sineOfGamma;
    forZ( 1, 0 ) = forZ( 0, 1 );
    forZ( 2, 0 ) = forZ( 0, 2 );
    forZ( 2, 1 ) = forZ( 1, 2 );
}

Eigen::Matrix3d computeEllipsoidalHarmonicPotentialHessian( const EllipsoidalHarmonicCoefficients& coefficients,
                                                            const EllipsoidalHarmonicCoefficients& radialFunctions,
                                                            const EllipsoidalHarmonicCoefficients& radialFunctionDerivatives,
                                                            const EllipsoidalHarmonicCoefficients& radialFunctionSecondDerivatives,
                                                            const EllipsoidalHarmonicMatrices& transformFactors,
                                                            const LegendreCache& legendreCache,
                                                            const Eigen::Vector3d& ellipsoidalCoordinates,
                                                            const double h,
                                                            const double k )
{
    const int maximumDegree = coefficients.getMaximumDegree( );
    const double sineOfPhi = std::sin( ellipsoidalCoordinates( 1 ) );
    const double cosineOfPhi = std::cos( ellipsoidalCoordinates( 1 ) );
    const double gamma = ellipsoidalCoordinates( 2 );

    Eigen::Vector3d gradientInEllipsoidalCoordinates = Eigen::Vector3d::Zero( );
    Eigen::Matrix3d hessianInEllipsoidalCoordinates = Eigen::Matrix3d::Zero( );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const int numberOfTerms = getNumberOfTermsInGroup( degree, group );

            if( numberOfTerms == 0 )
            {
                continue;
            }

            // The angular part of each term and its derivatives up to second order.
            Eigen::VectorXd surfaceTerms = Eigen::VectorXd::Zero( numberOfTerms );
            Eigen::VectorXd wrtPhi = Eigen::VectorXd::Zero( numberOfTerms );
            Eigen::VectorXd wrtGamma = Eigen::VectorXd::Zero( numberOfTerms );
            Eigen::VectorXd wrtPhiPhi = Eigen::VectorXd::Zero( numberOfTerms );
            Eigen::VectorXd wrtPhiGamma = Eigen::VectorXd::Zero( numberOfTerms );
            Eigen::VectorXd wrtGammaGamma = Eigen::VectorXd::Zero( numberOfTerms );

            for( int term = 0; term < numberOfTerms; term++ )
            {
                const int legendreOrder = getTermLegendreOrder( group, term );
                const double angle = legendreOrder * gamma;

                const double value = legendreCache.getLegendrePolynomial( degree, legendreOrder );
                // The cache differentiates with respect to sin( phi ), so the chain rule turns
                // those into derivatives with respect to phi.
                const double wrtSine = legendreCache.getLegendrePolynomialDerivative( degree, legendreOrder );
                const double wrtSineSine = legendreCache.getLegendrePolynomialSecondDerivative( degree, legendreOrder );

                const double legendreWrtPhi = wrtSine * cosineOfPhi;
                const double legendreWrtPhiPhi = wrtSineSine * cosineOfPhi * cosineOfPhi - wrtSine * sineOfPhi;

                const double angleFactor = termUsesSineOfAngle( group ) ? std::sin( angle ) : std::cos( angle );
                // Derivative of the sine or cosine of the angle with respect to gamma.
                const double angleFactorWrtGamma =
                        legendreOrder * ( termUsesSineOfAngle( group ) ? std::cos( angle ) : -std::sin( angle ) );
                const double angleFactorWrtGammaGamma = -legendreOrder * legendreOrder * angleFactor;

                surfaceTerms( term ) = value * angleFactor;
                wrtPhi( term ) = legendreWrtPhi * angleFactor;
                wrtGamma( term ) = value * angleFactorWrtGamma;
                wrtPhiPhi( term ) = legendreWrtPhiPhi * angleFactor;
                wrtPhiGamma( term ) = legendreWrtPhi * angleFactorWrtGamma;
                wrtGammaGamma( term ) = value * angleFactorWrtGammaGamma;
            }

            const Eigen::VectorXd& groupCoefficients = coefficients.getGroup( degree, group );
            const Eigen::MatrixXd& groupTransform = transformFactors.getGroup( degree, group );

            const Eigen::VectorXd weightedByValue = groupCoefficients.cwiseProduct( radialFunctions.getGroup( degree, group ) );
            const Eigen::VectorXd weightedByDerivative =
                    groupCoefficients.cwiseProduct( radialFunctionDerivatives.getGroup( degree, group ) );
            const Eigen::VectorXd weightedBySecondDerivative =
                    groupCoefficients.cwiseProduct( radialFunctionSecondDerivatives.getGroup( degree, group ) );

            gradientInEllipsoidalCoordinates( 0 ) += weightedByDerivative.dot( groupTransform * surfaceTerms );
            gradientInEllipsoidalCoordinates( 1 ) += weightedByValue.dot( groupTransform * wrtPhi );
            gradientInEllipsoidalCoordinates( 2 ) += weightedByValue.dot( groupTransform * wrtGamma );

            hessianInEllipsoidalCoordinates( 0, 0 ) += weightedBySecondDerivative.dot( groupTransform * surfaceTerms );
            hessianInEllipsoidalCoordinates( 0, 1 ) += weightedByDerivative.dot( groupTransform * wrtPhi );
            hessianInEllipsoidalCoordinates( 0, 2 ) += weightedByDerivative.dot( groupTransform * wrtGamma );
            hessianInEllipsoidalCoordinates( 1, 1 ) += weightedByValue.dot( groupTransform * wrtPhiPhi );
            hessianInEllipsoidalCoordinates( 1, 2 ) += weightedByValue.dot( groupTransform * wrtPhiGamma );
            hessianInEllipsoidalCoordinates( 2, 2 ) += weightedByValue.dot( groupTransform * wrtGammaGamma );
        }
    }

    hessianInEllipsoidalCoordinates( 1, 0 ) = hessianInEllipsoidalCoordinates( 0, 1 );
    hessianInEllipsoidalCoordinates( 2, 0 ) = hessianInEllipsoidalCoordinates( 0, 2 );
    hessianInEllipsoidalCoordinates( 2, 1 ) = hessianInEllipsoidalCoordinates( 1, 2 );

    // Change variables to Cartesian. With J the derivatives of position with respect to the
    // ellipsoidal coordinates and K its inverse, the Cartesian gradient is K^T times the
    // ellipsoidal one, and the Hessian picks up a term from the curvature of the coordinates.
    const Eigen::Matrix3d jacobian = computeEllipsoidalCoordinateJacobian( ellipsoidalCoordinates, h, k );
    const Eigen::Matrix3d inverseJacobian = jacobian.inverse( );

    const Eigen::Vector3d cartesianGradient = inverseJacobian.transpose( ) * gradientInEllipsoidalCoordinates;

    Eigen::Matrix3d hessian = inverseJacobian.transpose( ) * hessianInEllipsoidalCoordinates * inverseJacobian;

    std::vector< Eigen::Matrix3d > coordinateSecondDerivatives;
    computeEllipsoidalCoordinateSecondDerivatives( ellipsoidalCoordinates, h, k, coordinateSecondDerivatives );

    for( int component = 0; component < 3; component++ )
    {
        hessian -= cartesianGradient( component ) *
                ( inverseJacobian.transpose( ) * coordinateSecondDerivatives[ component ] * inverseJacobian );
    }

    return hessian;
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
