#define BOOST_TEST_MAIN

#include <cmath>
#include <functional>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/gravitation/mutualForcePotential.h"

namespace tudat
{

namespace unit_tests
{

namespace
{

double computeFactorial( const int number )
{
    double value = 1.0;
    for( int i = 2; i <= number; i++ )
    {
        value *= static_cast< double >( i );
    }
    return value;
}

double computeKroneckerDeltaZero( const int order )
{
    return ( order == 0 ) ? 1.0 : 0.0;
}

double computeSigmaSignFunctionFromLiterature( const int order )
{
    return ( order >= 0 ) ? 1.0 : ( ( std::abs( order ) % 2 == 0 ) ? 1.0 : -1.0 );
}

double computeSFunctionFromLiterature( const int order )
{
    return ( order > 0 ) ? 1.0 : ( ( order < 0 ) ? -1.0 : 0.0 );
}

double computeGammaFromLiterature( const int degree1, const int order1, const int degree2, const int order2 )
{
    // Require valid body-1 degree/order input before evaluating the literature factorial expression.
    BOOST_REQUIRE_LE( std::abs( order1 ), degree1 );
    // Require valid body-2 degree/order input before evaluating the literature factorial expression.
    BOOST_REQUIRE_LE( std::abs( order2 ), degree2 );

    if( ( degree1 == 0 && order1 == 0 ) || ( degree2 == 0 && order2 == 0 ) )
    {
        return 1.0 / std::sqrt( 4.0 * mathematical_constants::PI );
    }

    const double numerator = ( 2.0 * static_cast< double >( degree1 ) + 1.0 ) * ( 2.0 * static_cast< double >( degree2 ) + 1.0 ) *
            computeFactorial( degree1 + degree2 - order1 - order2 ) * computeFactorial( degree1 + degree2 + order1 + order2 );

    const double denominator = computeFactorial( degree1 + order1 ) * computeFactorial( degree1 - order1 ) *
            computeFactorial( degree2 + order2 ) * computeFactorial( degree2 - order2 ) * 4.0 * mathematical_constants::PI *
            static_cast< double >( 2 * degree1 + 2 * degree2 + 1 );

    return std::sqrt( numerator / denominator );
}

double computeNormalizedEffectiveMultiplierFromLiterature( const int degree1, const int order1, const int degree2, const int order2 )
{
    // Require valid body-1 degree/order input before evaluating the normalized multiplier.
    BOOST_REQUIRE_LE( std::abs( order1 ), degree1 );
    // Require valid body-2 degree/order input before evaluating the normalized multiplier.
    BOOST_REQUIRE_LE( std::abs( order2 ), degree2 );

    // Literature sources used for this factor:
    // - Compere & Lemaître (2014): cross-body gamma normalization.
    // - Dirkx et al. (2019): normalized gamma-bar from Eq. (21), followed by the
    //   real coefficient scaling from Eqs. (47)-(48).
    const double gammaTilde = computeGammaFromLiterature( degree1, order1, degree2, order2 );
    const double twoMinusDeltaOrder1 = 2.0 - computeKroneckerDeltaZero( order1 );
    const double twoMinusDeltaOrder2 = 2.0 - computeKroneckerDeltaZero( order2 );
    const double twoMinusDeltaCombinedOrder = 2.0 - computeKroneckerDeltaZero( order1 + order2 );
    const double onePlusDeltaOrder1 = 1.0 + computeKroneckerDeltaZero( order1 );
    const double onePlusDeltaOrder2 = 1.0 + computeKroneckerDeltaZero( order2 );
    const double parityFactor = ( degree1 % 2 == 0 ) ? 1.0 : -1.0;

    const double gammaBar = parityFactor * gammaTilde *
            std::sqrt( 4.0 * mathematical_constants::PI * twoMinusDeltaOrder1 * twoMinusDeltaOrder2 / twoMinusDeltaCombinedOrder ) *
            computeSigmaSignFunctionFromLiterature( order1 + order2 );

    return computeSigmaSignFunctionFromLiterature( order1 ) * computeSigmaSignFunctionFromLiterature( order2 ) * gammaBar *
            onePlusDeltaOrder1 * onePlusDeltaOrder2 / 4.0;
}

std::pair< double, double > computeExpectedEffectiveCoefficientsFromLiterature( const Eigen::MatrixXd& cosineBody1Coefficients,
                                                                                const Eigen::MatrixXd& sineBody1Coefficients,
                                                                                const Eigen::MatrixXd& transformedCosineBody2Coefficients,
                                                                                const Eigen::MatrixXd& transformedSineBody2Coefficients,
                                                                                const int degree1,
                                                                                const int order1,
                                                                                const int degree2,
                                                                                const int order2 )
{
    // Require valid body-1 degree/order input before reading the coefficient arrays.
    BOOST_REQUIRE_LE( std::abs( order1 ), degree1 );
    // Require valid body-2 degree/order input before reading the coefficient arrays.
    BOOST_REQUIRE_LE( std::abs( order2 ), degree2 );

    const double cosineBody1 = cosineBody1Coefficients( degree1, std::abs( order1 ) );
    const double sineBody1 = sineBody1Coefficients( degree1, std::abs( order1 ) );
    const double cosineBody2 = transformedCosineBody2Coefficients( degree2, std::abs( order2 ) );
    const double sineBody2 = transformedSineBody2Coefficients( degree2, std::abs( order2 ) );

    // Literature form for real effective coefficients:
    // C_eff = nu * ( C1*C2 - s_m1*s_m2*S1*S2 )
    // S_eff = s_(m1+m2) * nu * ( s_m2*C1*S2 + s_m1*S1*C2 )
    const double multiplier = computeNormalizedEffectiveMultiplierFromLiterature( degree1, order1, degree2, order2 );
    const double sOrder1 = computeSFunctionFromLiterature( order1 );
    const double sOrder2 = computeSFunctionFromLiterature( order2 );
    const double sCombinedOrder = computeSFunctionFromLiterature( order1 + order2 );

    const double expectedCosine = ( cosineBody1 * cosineBody2 - sOrder1 * sOrder2 * sineBody1 * sineBody2 ) * multiplier;
    const double expectedSine = sCombinedOrder * ( sOrder2 * cosineBody1 * sineBody2 + sOrder1 * sineBody1 * cosineBody2 ) * multiplier;

    return std::make_pair( expectedCosine, expectedSine );
}

void validateSingleCombinationAgainstLiterature( gravitation::EffectiveMutualSphericalHarmonicsField& mutualField,
                                                 const Eigen::MatrixXd& cosineBody1Coefficients,
                                                 const Eigen::MatrixXd& sineBody1Coefficients,
                                                 const Eigen::MatrixXd& transformedCosineBody2Coefficients,
                                                 const Eigen::MatrixXd& transformedSineBody2Coefficients,
                                                 const int degree1,
                                                 const int order1,
                                                 const int degree2,
                                                 const int order2,
                                                 const double tolerance )
{
    const int effectiveIndex = mutualField.getEffectiveIndex( degree1, order1, degree2, order2 );

    double computedCosine = TUDAT_NAN;
    double computedSine = TUDAT_NAN;
    mutualField.getCurrentEffectiveCoefficients( degree1, order1, degree2, order2, effectiveIndex, computedCosine, computedSine );

    const std::pair< double, double > expectedCoefficients =
            computeExpectedEffectiveCoefficientsFromLiterature( cosineBody1Coefficients,
                                                                sineBody1Coefficients,
                                                                transformedCosineBody2Coefficients,
                                                                transformedSineBody2Coefficients,
                                                                degree1,
                                                                order1,
                                                                degree2,
                                                                order2 );

    // Verify the effective cosine coefficient against the independent literature expression.
    BOOST_CHECK_SMALL( std::fabs( computedCosine - expectedCoefficients.first ), tolerance );
    // The m=0 sine basis function is identically zero in the potential. The
    // physical coefficient check for signed pairs with m1 + m2 = 0 is therefore
    // carried by the cosine coefficient.
    if( order1 + order2 != 0 )
    {
        // Verify the effective sine coefficient when the associated sine basis is not identically zero.
        BOOST_CHECK_SMALL( std::fabs( computedSine - expectedCoefficients.second ), tolerance );
    }
}

void validateRequestedCombinationAgainstLiterature( gravitation::EffectiveMutualSphericalHarmonicsField& mutualField,
                                                    const Eigen::MatrixXd& cosineBody1Coefficients,
                                                    const Eigen::MatrixXd& sineBody1Coefficients,
                                                    const Eigen::MatrixXd& transformedCosineBody2Coefficients,
                                                    const Eigen::MatrixXd& transformedSineBody2Coefficients,
                                                    const int degree1,
                                                    const int order1,
                                                    const int degree2,
                                                    const int order2,
                                                    const double tolerance )
{
    for( int combinationCase = 0; combinationCase < 4; combinationCase++ )
    {
        int signedOrder1 = 0;
        int signedOrder2 = 0;
        const bool useCombination =
                gravitation::getSignedOrdersForCombinationCase( combinationCase, order1, order2, signedOrder1, signedOrder2 );
        if( useCombination )
        {
            validateSingleCombinationAgainstLiterature( mutualField,
                                                        cosineBody1Coefficients,
                                                        sineBody1Coefficients,
                                                        transformedCosineBody2Coefficients,
                                                        transformedSineBody2Coefficients,
                                                        degree1,
                                                        signedOrder1,
                                                        degree2,
                                                        signedOrder2,
                                                        tolerance );
        }
    }
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_MutualForcePotential )

BOOST_AUTO_TEST_CASE( testEffectiveMutualCoefficientsForRequestedDegreeOrderCombinations )
{
    // Test rationale:
    // Validate the effective mutual coefficients C_eff and S_eff produced by
    // EffectiveMutualSphericalHarmonicsField against closed-form literature expressions:
    // 1) cross-body normalization term gamma (Compere & Lemaître, 2014),
    // 2) effective mutual coefficient combination/scaling in normalized form (Dirkx et al., 2019).
    //
    // What is tested:
    // - (0,0)x(2,m), (2,m)x(0,0), and (2,m1)x(2,m2) combinations.
    // - identity and arbitrary relative orientations (through transformed body-2 coefficients).
    //
    // Why this matters:
    // the full-two-body acceleration/torque models depend directly on these effective coefficients;
    // a mismatch here propagates into all figure-figure interaction terms.
    Eigen::MatrixXd cosineBody1Coefficients = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineBody1Coefficients = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd cosineBody2Coefficients = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineBody2Coefficients = Eigen::MatrixXd::Zero( 3, 3 );

    cosineBody1Coefficients( 0, 0 ) = 1.0;
    cosineBody2Coefficients( 0, 0 ) = 1.0;

    cosineBody1Coefficients( 2, 0 ) = 0.21;
    cosineBody1Coefficients( 2, 1 ) = -0.11;
    sineBody1Coefficients( 2, 1 ) = 0.14;
    cosineBody1Coefficients( 2, 2 ) = 0.08;
    sineBody1Coefficients( 2, 2 ) = -0.09;

    cosineBody2Coefficients( 2, 0 ) = -0.17;
    cosineBody2Coefficients( 2, 1 ) = 0.16;
    sineBody2Coefficients( 2, 1 ) = -0.12;
    cosineBody2Coefficients( 2, 2 ) = 0.07;
    sineBody2Coefficients( 2, 2 ) = 0.19;

    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations;

    for( int m = 0; m <= 2; m++ )
    {
        coefficientCombinations.push_back( std::make_tuple( 0, 0, 2, m ) );
    }
    for( int m = 0; m <= 2; m++ )
    {
        coefficientCombinations.push_back( std::make_tuple( 2, m, 0, 0 ) );
    }
    for( int m1 = 0; m1 <= 2; m1++ )
    {
        for( int m2 = 0; m2 <= 2; m2++ )
        {
            coefficientCombinations.push_back( std::make_tuple( 2, m1, 2, m2 ) );
        }
    }

    gravitation::EffectiveMutualSphericalHarmonicsField mutualField(
            coefficientCombinations,
            [ cosineBody1Coefficients ]( ) { return cosineBody1Coefficients; },
            [ sineBody1Coefficients ]( ) { return sineBody1Coefficients; },
            [ cosineBody2Coefficients ]( ) { return cosineBody2Coefficients; },
            [ sineBody2Coefficients ]( ) { return sineBody2Coefficients; },
            []( ) { return 1.0; },
            1.0,
            1.0,
            true );

    const double tolerance = 25.0 * std::numeric_limits< double >::epsilon( );

    std::vector< Eigen::Matrix3d > coefficientRotations;
    coefficientRotations.push_back( Eigen::Matrix3d::Identity( ) );
    coefficientRotations.push_back(
            ( Eigen::AngleAxisd( 0.63, Eigen::Vector3d::UnitX( ) ) * Eigen::AngleAxisd( -0.37, Eigen::Vector3d::UnitY( ) ) *
              Eigen::AngleAxisd( 1.11, Eigen::Vector3d::UnitZ( ) ) )
                    .toRotationMatrix( ) );

    for( const Eigen::Matrix3d& currentRotationMatrix : coefficientRotations )
    {
        mutualField.computeCurrentEffectiveCoefficients( Eigen::Quaterniond( currentRotationMatrix ) );

        const Eigen::MatrixXd transformedCosineBody2Coefficients = mutualField.getTransformedCosineCoefficientsOfBody2( );
        const Eigen::MatrixXd transformedSineBody2Coefficients = mutualField.getTransformedSineCoefficientsOfBody2( );

        // Case 1: l=m=0 for body 1, l=2 and m=0..2 for body 2
        for( int m2 = 0; m2 <= 2; m2++ )
        {
            validateRequestedCombinationAgainstLiterature( mutualField,
                                                           cosineBody1Coefficients,
                                                           sineBody1Coefficients,
                                                           transformedCosineBody2Coefficients,
                                                           transformedSineBody2Coefficients,
                                                           0,
                                                           0,
                                                           2,
                                                           m2,
                                                           tolerance );
        }

        // Case 2: l=m=0 for body 2, l=2 and m=0..2 for body 1
        for( int m1 = 0; m1 <= 2; m1++ )
        {
            validateRequestedCombinationAgainstLiterature( mutualField,
                                                           cosineBody1Coefficients,
                                                           sineBody1Coefficients,
                                                           transformedCosineBody2Coefficients,
                                                           transformedSineBody2Coefficients,
                                                           2,
                                                           m1,
                                                           0,
                                                           0,
                                                           tolerance );
        }

        // Case 3: l=2 for both bodies and m1,m2=0..2
        for( int m1 = 0; m1 <= 2; m1++ )
        {
            for( int m2 = 0; m2 <= 2; m2++ )
            {
                validateRequestedCombinationAgainstLiterature( mutualField,
                                                               cosineBody1Coefficients,
                                                               sineBody1Coefficients,
                                                               transformedCosineBody2Coefficients,
                                                               transformedSineBody2Coefficients,
                                                               2,
                                                               m1,
                                                               2,
                                                               m2,
                                                               tolerance );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
