#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <algorithm>
#include <memory>
#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/simulation/simulation.h"
#include "tudat/math/statistics/randomVariableGenerator.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::statistics;
using namespace tudat::gravitation;

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > generateCosineSineCoefficients(
        std::shared_ptr< RandomVariableGenerator< double > > randomNumberGenerator,
        const int maximumDegree,
        const int maximumOrder )
{
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );

    for( int i = 0; i <= maximumDegree; i++ )
    {
        for( int j = 0; ( ( j <= i ) && ( j <= maximumOrder ) ); j++ )
        {
            cosineCoefficients( i, j ) = randomNumberGenerator->getRandomVariableValue( );
            if( j != 0 )
            {
                sineCoefficients( i, j ) = randomNumberGenerator->getRandomVariableValue( );
            }
        }
    }

    // cosineCoefficients( 1, 0 ) = randomNumberGenerator->getRandomVariableValue( );
    //     sineCoefficients( 2, 1 ) = randomNumberGenerator->getRandomVariableValue( );

    cosineCoefficients( 0, 0 ) = 1.0;

    return std::make_pair( cosineCoefficients, sineCoefficients );
}

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getDegreeTwoDegreeTwoInteractions( )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse;
    for( unsigned int m1 = 0; m1 <= 2; m1++ )
    {
        for( unsigned int m2 = 0; m2 <= 2; m2++ )
        {
            coefficientCombinationsToUse.push_back( std::make_tuple( 2, m1, 2, m2 ) );
        }
    }
    return coefficientCombinationsToUse;
}

std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > createMutualExtendedBodyAccelerationModel(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const Eigen::Vector3d& positionOfBody1,
        const Eigen::Vector3d& positionOfBody2,
        const double gravitationalParameter,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const Eigen::MatrixXd& cosineCoefficientsOfBody1,
        const Eigen::MatrixXd& sineCoefficientsOfBody1,
        const Eigen::MatrixXd& cosineCoefficientsOfBody2,
        const Eigen::MatrixXd& sineCoefficientsOfBody2,
        const Eigen::Quaterniond& rotationToBody1,
        const Eigen::Quaterniond& rotationToBody2 )
{
    return std::make_shared< FullTwoBodySphericalHarmonicAcceleration >( [ = ]( ) { return positionOfBody1; },
                                                                         [ = ]( ) { return positionOfBody2; },
                                                                         [ = ]( ) { return gravitationalParameter; },
                                                                         equatorialRadiusOfBody1,
                                                                         equatorialRadiusOfBody2,
                                                                         [ = ]( ) { return cosineCoefficientsOfBody1; },
                                                                         [ = ]( ) { return sineCoefficientsOfBody1; },
                                                                         [ = ]( ) { return cosineCoefficientsOfBody2; },
                                                                         [ = ]( ) { return sineCoefficientsOfBody2; },
                                                                         coefficientCombinationsToUse,
                                                                         [ = ]( ) { return rotationToBody1; },
                                                                         [ = ]( ) { return rotationToBody2; },
                                                                         false,
                                                                         true );
}

void runDegreeTwoCrossTermValidationCase( const Eigen::Vector3d& positionOfBody1,
                                          const Eigen::Quaterniond& rotationToBody1,
                                          const Eigen::Quaterniond& rotationToBody2,
                                          const Eigen::MatrixXd& cosineCoefficientsOfBody1,
                                          const Eigen::MatrixXd& sineCoefficientsOfBody1,
                                          const Eigen::MatrixXd& cosineCoefficientsOfBody2,
                                          const Eigen::MatrixXd& sineCoefficientsOfBody2,
                                          const double currentTime,
                                          const double gravitationalParameter,
                                          const double equatorialRadiusOfBody1,
                                          const double equatorialRadiusOfBody2 )
{
    // Helper validation for a single geometry/orientation case:
    // Isolate and verify l1=2,l2=2 figure-figure acceleration terms from the full interaction summation
    // (Dirkx et al., 2019 mutual-potential interaction decomposition).
    const Eigen::Vector3d positionOfBody2 = Eigen::Vector3d::Zero( );

    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > fullDegreeTwoInteractions =
            FullTwoBodySphericalHarmonicAccelerationSettings( 2, 2, 2, 2 ).coefficientCombinationsToUse_;
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body1OnlyInteractions =
            FullTwoBodySphericalHarmonicAccelerationSettings( 2, 2, 0, 0 ).coefficientCombinationsToUse_;
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2OnlyInteractions =
            FullTwoBodySphericalHarmonicAccelerationSettings( 0, 0, 2, 2 ).coefficientCombinationsToUse_;
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > centralInteractions =
            FullTwoBodySphericalHarmonicAccelerationSettings( 0, 0, 0, 0 ).coefficientCombinationsToUse_;
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > degreeTwoCrossInteractions =
            getDegreeTwoDegreeTwoInteractions( );

    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > fullDegreeTwoAcceleration =
            createMutualExtendedBodyAccelerationModel( fullDegreeTwoInteractions,
                                                       positionOfBody1,
                                                       positionOfBody2,
                                                       gravitationalParameter,
                                                       equatorialRadiusOfBody1,
                                                       equatorialRadiusOfBody2,
                                                       cosineCoefficientsOfBody1,
                                                       sineCoefficientsOfBody1,
                                                       cosineCoefficientsOfBody2,
                                                       sineCoefficientsOfBody2,
                                                       rotationToBody1,
                                                       rotationToBody2 );

    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > body1OnlyAcceleration =
            createMutualExtendedBodyAccelerationModel( body1OnlyInteractions,
                                                       positionOfBody1,
                                                       positionOfBody2,
                                                       gravitationalParameter,
                                                       equatorialRadiusOfBody1,
                                                       equatorialRadiusOfBody2,
                                                       cosineCoefficientsOfBody1,
                                                       sineCoefficientsOfBody1,
                                                       cosineCoefficientsOfBody2,
                                                       sineCoefficientsOfBody2,
                                                       rotationToBody1,
                                                       rotationToBody2 );

    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > body2OnlyAcceleration =
            createMutualExtendedBodyAccelerationModel( body2OnlyInteractions,
                                                       positionOfBody1,
                                                       positionOfBody2,
                                                       gravitationalParameter,
                                                       equatorialRadiusOfBody1,
                                                       equatorialRadiusOfBody2,
                                                       cosineCoefficientsOfBody1,
                                                       sineCoefficientsOfBody1,
                                                       cosineCoefficientsOfBody2,
                                                       sineCoefficientsOfBody2,
                                                       rotationToBody1,
                                                       rotationToBody2 );

    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > centralAcceleration =
            createMutualExtendedBodyAccelerationModel( centralInteractions,
                                                       positionOfBody1,
                                                       positionOfBody2,
                                                       gravitationalParameter,
                                                       equatorialRadiusOfBody1,
                                                       equatorialRadiusOfBody2,
                                                       cosineCoefficientsOfBody1,
                                                       sineCoefficientsOfBody1,
                                                       cosineCoefficientsOfBody2,
                                                       sineCoefficientsOfBody2,
                                                       rotationToBody1,
                                                       rotationToBody2 );

    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > degreeTwoCrossAcceleration =
            createMutualExtendedBodyAccelerationModel( degreeTwoCrossInteractions,
                                                       positionOfBody1,
                                                       positionOfBody2,
                                                       gravitationalParameter,
                                                       equatorialRadiusOfBody1,
                                                       equatorialRadiusOfBody2,
                                                       cosineCoefficientsOfBody1,
                                                       sineCoefficientsOfBody1,
                                                       cosineCoefficientsOfBody2,
                                                       sineCoefficientsOfBody2,
                                                       rotationToBody1,
                                                       rotationToBody2 );

    fullDegreeTwoAcceleration->updateMembers( currentTime );
    body1OnlyAcceleration->updateMembers( currentTime );
    body2OnlyAcceleration->updateMembers( currentTime );
    centralAcceleration->updateMembers( currentTime );
    degreeTwoCrossAcceleration->updateMembers( currentTime );

    // Check 1: isolate l1=2,l2=2 terms from the full degree-2 model.
    const Eigen::Vector3d isolatedCrossTermAcceleration = fullDegreeTwoAcceleration->getAcceleration( ) -
            body1OnlyAcceleration->getAcceleration( ) - body2OnlyAcceleration->getAcceleration( ) + centralAcceleration->getAcceleration( );
    const Eigen::Vector3d directCrossTermAcceleration = degreeTwoCrossAcceleration->getAcceleration( );
    const Eigen::Vector3d isolationDifference = isolatedCrossTermAcceleration - directCrossTermAcceleration;
    const double isolationScale = std::max( 1.0, isolatedCrossTermAcceleration.norm( ) );
    BOOST_CHECK_SMALL( isolationDifference.norm( ) / isolationScale, 5.0E-14 );

    // Check 2: cross-term acceleration scales linearly with degree-2 coefficients of body 2.
    const double scalingFactor = -1.35;
    Eigen::MatrixXd scaledCosineCoefficientsOfBody2 = cosineCoefficientsOfBody2;
    Eigen::MatrixXd scaledSineCoefficientsOfBody2 = sineCoefficientsOfBody2;
    for( int m = 0; m <= 2; m++ )
    {
        scaledCosineCoefficientsOfBody2( 2, m ) *= scalingFactor;
        if( m > 0 )
        {
            scaledSineCoefficientsOfBody2( 2, m ) *= scalingFactor;
        }
    }

    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > scaledCrossAcceleration =
            createMutualExtendedBodyAccelerationModel( degreeTwoCrossInteractions,
                                                       positionOfBody1,
                                                       positionOfBody2,
                                                       gravitationalParameter,
                                                       equatorialRadiusOfBody1,
                                                       equatorialRadiusOfBody2,
                                                       cosineCoefficientsOfBody1,
                                                       sineCoefficientsOfBody1,
                                                       scaledCosineCoefficientsOfBody2,
                                                       scaledSineCoefficientsOfBody2,
                                                       rotationToBody1,
                                                       rotationToBody2 );
    scaledCrossAcceleration->updateMembers( currentTime );

    const Eigen::Vector3d scaledDifference = scaledCrossAcceleration->getAcceleration( ) - scalingFactor * directCrossTermAcceleration;
    const double scalingReference = std::max( 1.0, directCrossTermAcceleration.norm( ) );
    BOOST_CHECK_SMALL( scaledDifference.norm( ) / scalingReference, 5.0E-14 );

    // Check 3: zero degree-2 coefficients of body 2 must null the degree-2 cross-term acceleration.
    Eigen::MatrixXd zeroedCosineCoefficientsOfBody2 = cosineCoefficientsOfBody2;
    Eigen::MatrixXd zeroedSineCoefficientsOfBody2 = sineCoefficientsOfBody2;
    for( int m = 0; m <= 2; m++ )
    {
        zeroedCosineCoefficientsOfBody2( 2, m ) = 0.0;
        if( m > 0 )
        {
            zeroedSineCoefficientsOfBody2( 2, m ) = 0.0;
        }
    }

    const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > zeroedCrossAcceleration =
            createMutualExtendedBodyAccelerationModel( degreeTwoCrossInteractions,
                                                       positionOfBody1,
                                                       positionOfBody2,
                                                       gravitationalParameter,
                                                       equatorialRadiusOfBody1,
                                                       equatorialRadiusOfBody2,
                                                       cosineCoefficientsOfBody1,
                                                       sineCoefficientsOfBody1,
                                                       zeroedCosineCoefficientsOfBody2,
                                                       zeroedSineCoefficientsOfBody2,
                                                       rotationToBody1,
                                                       rotationToBody2 );
    zeroedCrossAcceleration->updateMembers( currentTime );

    BOOST_CHECK_SMALL( zeroedCrossAcceleration->getAcceleration( ).norm( ), 1.0E-22 );
}

std::shared_ptr< tudat::simulation_setup::GravityFieldSettings > getDummyJovianSystemGravityField( const std::string& bodyName,
                                                                                                   const int isNormalized )
{
    std::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    std::vector< double > randomNumberSettings;
    randomNumberSettings.push_back( 0.0 );
    randomNumberSettings.push_back( 1.0E-4 );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;

    if( bodyName == "Jupiter" )
    {
        std::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator =
                createBoostContinuousRandomVariableGenerator( normal_boost_distribution, randomNumberSettings, 0.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >( getBodyGravitationalParameter( "Jupiter" ),
                                                                                           getAverageRadius( "Jupiter" ),
                                                                                           coefficients.first,
                                                                                           coefficients.second,
                                                                                           "IAU_Jupiter_SIMPLIFIED" );
    }
    else if( bodyName == "Io" )
    {
        std::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator =
                createBoostContinuousRandomVariableGenerator( normal_boost_distribution, randomNumberSettings, 1.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                5.959916033410404E012, getAverageRadius( "Io" ), coefficients.first, coefficients.second, "IAU_Io_SIMPLIFIED" );
    }
    else if( bodyName == "Europa" )
    {
        std::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator =
                createBoostContinuousRandomVariableGenerator( normal_boost_distribution, randomNumberSettings, 2.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                3.202738774922892E12, getAverageRadius( "Europa" ), coefficients.first, coefficients.second, "IAU_Europa_SIMPLIFIED" );
    }

    return gravityFieldSettings;
}

// void getBodyCoefficientEulerAnglePartials(
//         const Eigen::Quaterniond nominalQuaternion,
//         std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualShField,
//         std::vector< Eigen::MatrixXd >& numericalTransformedCosineCoefficientsOfBody2Partials,
//         std::vector< Eigen::MatrixXd >& numericalTransformedSineCoefficientsOfBody2Partials,
//         const double perturbation )
//{
//     numericalTransformedCosineCoefficientsOfBody2Partials.resize( 3 );
//     numericalTransformedSineCoefficientsOfBody2Partials.resize( 3 );

//    Eigen::Vector3d nominalEulerAngles = basic_mathematics::get313EulerAnglesFromQuaternion( nominalQuaternion );
//    Eigen::Vector3d perturbedEulerAngles;
//    Eigen::MatrixXd upPerturbedCosineCoefficients, upPerturbedSineCoefficients;
//    Eigen::MatrixXd downPerturbedCosineCoefficients, downPerturbedSineCoefficients;

//    for( unsigned int i = 0; i < 3; i ++ )
//    {
//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) += perturbation;
//        mutualShField->computeCurrentEffectiveCoefficients( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0
//        ) ); upPerturbedCosineCoefficients = mutualShField->getTransformedCosineCoefficientsOfBody2( ); upPerturbedSineCoefficients =
//        mutualShField->getTransformedSineCoefficientsOfBody2( );

//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) -= perturbation;
//        mutualShField->computeCurrentEffectiveCoefficients( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0
//        ) ); downPerturbedCosineCoefficients = mutualShField->getTransformedCosineCoefficientsOfBody2( ); downPerturbedSineCoefficients =
//        mutualShField->getTransformedSineCoefficientsOfBody2( );

//        numericalTransformedCosineCoefficientsOfBody2Partials[ i ] =
//                ( upPerturbedCosineCoefficients - downPerturbedCosineCoefficients ) / ( 2.0 * perturbation );
//        numericalTransformedSineCoefficientsOfBody2Partials[ i ] =
//                ( upPerturbedSineCoefficients - downPerturbedSineCoefficients ) / ( 2.0 * perturbation );

//    }
//}

Eigen::Matrix2d getEffectiveMutualPotentialCoefficientWrtBody2Coefficient(
        const int degreeOfBody1,
        const int orderOfBody1,
        const int degreeOfBody2,
        const int orderOfBody2,
        const Eigen::MatrixXd& nominalTransformedCosineCoefficients,
        const Eigen::MatrixXd& nominalTransformedSineCoefficients,
        const std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualPotentialField,
        const double coefficientPerturbation )
{
    Eigen::MatrixXd perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    Eigen::MatrixXd perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedCosineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) += coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients( perturbedTransformedCosineCoefficients,
                                                                                                perturbedTransformedSineCoefficients );
    double upPerturbedCosineCoefficientFromCosinePerturbation =
            mutualPotentialField->getEffectiveCosineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double upPerturbedSineCoefficientFromCosinePerturbation =
            mutualPotentialField->getEffectiveSineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedCosineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) -= coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients( perturbedTransformedCosineCoefficients,
                                                                                                perturbedTransformedSineCoefficients );
    double downPerturbedCosineCoefficientFromCosinePerturbation =
            mutualPotentialField->getEffectiveCosineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double downPerturbedSineCoefficientFromCosinePerturbation =
            mutualPotentialField->getEffectiveSineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedSineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) += coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients( perturbedTransformedCosineCoefficients,
                                                                                                perturbedTransformedSineCoefficients );
    double upPerturbedCosineCoefficientFromSinePerturbation =
            mutualPotentialField->getEffectiveCosineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double upPerturbedSineCoefficientFromSinePerturbation =
            mutualPotentialField->getEffectiveSineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedSineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) -= coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients( perturbedTransformedCosineCoefficients,
                                                                                                perturbedTransformedSineCoefficients );
    double downPerturbedCosineCoefficientFromSinePerturbation =
            mutualPotentialField->getEffectiveCosineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double downPerturbedSineCoefficientFromSinePerturbation =
            mutualPotentialField->getEffectiveSineCoefficient( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

    return ( ( Eigen::Matrix2d( ) << ( upPerturbedCosineCoefficientFromCosinePerturbation -
                                       downPerturbedCosineCoefficientFromCosinePerturbation ),
               ( upPerturbedCosineCoefficientFromSinePerturbation - downPerturbedCosineCoefficientFromSinePerturbation ),
               ( upPerturbedSineCoefficientFromCosinePerturbation - downPerturbedSineCoefficientFromCosinePerturbation ),
               ( upPerturbedSineCoefficientFromSinePerturbation - downPerturbedSineCoefficientFromSinePerturbation ) )
                     .finished( ) ) /
            ( 2.0 * coefficientPerturbation );
}

BOOST_AUTO_TEST_SUITE( test_mutual_spherical_harmonic_gravity )

BOOST_AUTO_TEST_CASE( testMutualSphericalHarmonicGravity )
{
    // Test rationale:
    // Cross-validate the full-two-body spherical-harmonic acceleration model against independent model families:
    // 1) one-way spherical harmonics,
    // 2) legacy mutual spherical harmonics (single-point interactions),
    // 3) third-body mutual-extended model decomposition.
    //
    // Why this matters:
    // these equivalence checks confirm that specific interaction subsets selected in
    // FullTwoBodySphericalHarmonicAccelerationSettings collapse to known reference formulations.
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialTime = 1.0E7;
    double finalTime = 1.2E7;

    int expansionDegree = 5;

    // for( int isNormalized = 1; isNormalized < 1; isNormalized++ )
    int isNormalized = 1;
    {
        // Get body settings.
        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialTime, finalTime );
        bodySettings.at( "Jupiter" )->gravityFieldSettings = getDummyJovianSystemGravityField( "Jupiter", isNormalized );
        bodySettings.at( "Io" )->gravityFieldSettings = getDummyJovianSystemGravityField( "Io", isNormalized );
        bodySettings.at( "Europa" )->gravityFieldSettings = getDummyJovianSystemGravityField( "Europa", isNormalized );

        bodySettings.at( "Jupiter" )->rotationModelSettings =
                std::make_shared< SimpleRotationModelSettings >( "ECLIPJ2000",
                                                                 "IAU_Jupiter_SIMPLIFIED",
                                                                 Eigen::Quaterniond( Eigen::AngleAxisd( -1.4, Eigen::Vector3d::UnitZ( ) ) *
                                                                                     Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitX( ) ) *
                                                                                     Eigen::AngleAxisd( 2.6, Eigen::Vector3d::UnitZ( ) ) ),
                                                                 0.0,
                                                                 0.0 );

        bodySettings.at( "Io" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000",
                "IAU_Io_SIMPLIFIED",
                Eigen::Quaterniond( Eigen::AngleAxisd( mathematical_constants::PI * 3.0 / 8.0, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( mathematical_constants::PI * 5.0 / 8.0, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( mathematical_constants::PI * -2.5 / 8.0, Eigen::Vector3d::UnitZ( ) ) ),
                0.0,
                0.0 );
        bodySettings.at( "Europa" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000",
                "IAU_Europa_SIMPLIFIED",
                Eigen::Quaterniond( Eigen::AngleAxisd( mathematical_constants::PI * 3.0 / 8.0, Eigen::Vector3d::UnitZ( ) ) ),
                0.0,
                0.0 );

        // Create bodies needed in simulation
        SystemOfBodies bodyMap = createSystemOfBodies( bodySettings );

        // Set current state and rotation of bodies.
        double currentTime = 1.1E7;
        bodyMap.at( "Jupiter" )->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap.at( "Jupiter" )->setStateFromEphemeris( currentTime );
        bodyMap.at( "Io" )->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap.at( "Io" )->setStateFromEphemeris( currentTime );
        bodyMap.at( "Europa" )->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap.at( "Europa" )->setStateFromEphemeris( currentTime );

        // Retrieve gravity fields.
        std::shared_ptr< SphericalHarmonicsGravityField > jupiterGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( ( bodyMap.at( "Jupiter" ) )->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > ioGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( ( bodyMap.at( "Io" ) )->getGravityFieldModel( ) );
        std::shared_ptr< SphericalHarmonicsGravityField > europaGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( ( bodyMap.at( "Europa" ) )->getGravityFieldModel( ) );

        // Create central gravity acceleration (mu = Io + Jupiter)
        std::shared_ptr< AccelerationSettings > centralGravitySettings =
                std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity );
        std::shared_ptr< CentralGravitationalAccelerationModel3d > centralGravity =
                std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >( createAccelerationModel( bodyMap.at( "Io" ),
                                                                                                               bodyMap.at( "Jupiter" ),
                                                                                                               centralGravitySettings,
                                                                                                               "Io",
                                                                                                               "Jupiter",
                                                                                                               bodyMap.at( "Jupiter" ),
                                                                                                               "Jupiter" ) );

        // Calculate central gravity acceleration.
        centralGravity->updateMembers( );

        // Create spherical harmonic gravity of Jupiter on Io, Jupiter-fixed (mu = Io + Jupiter)
        std::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnIoFromJupiterSettings =
                std::make_shared< SphericalHarmonicAccelerationSettings >( expansionDegree, expansionDegree );
        std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnIoFromJupiter =
                std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >( createAccelerationModel(
                        bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), sphericalHarmonicGravityOnIoFromJupiterSettings, "Io", "Jupiter" ) );

        // Calculate spherical harmonic gravity of Jupiter on Io.
        sphericalHarmonicGravityOnIoFromJupiter->updateMembers( );
        Eigen::Vector3d sphericalHarmonicGravityOnIoFromJupiterAcceleration = sphericalHarmonicGravityOnIoFromJupiter->getAcceleration( );
        double sphericalHarmonicPotentialAtIoFromJupiter = jupiterGravityField->getGravitationalPotential(
                sphericalHarmonicGravityOnIoFromJupiter->getCurrentRelativePosition( ), expansionDegree, expansionDegree );

        // Create spherical harmonic gravity of Io on Jupiter, Io-fixed (mu = Io + Jupiter)
        std::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnJupiterFromIoSettings =
                std::make_shared< SphericalHarmonicAccelerationSettings >( expansionDegree, expansionDegree );
        std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnJupiterFromIo =
                std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >( createAccelerationModel(
                        bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), sphericalHarmonicGravityOnJupiterFromIoSettings, "Jupiter", "Io" ) );

        // Calculate spherical harmonic gravity of Io on Jupiter.
        sphericalHarmonicGravityOnJupiterFromIo->updateMembers( );
        Eigen::Vector3d sphericalHarmonicGravityOnJupiterFromIoAcceleration = sphericalHarmonicGravityOnJupiterFromIo->getAcceleration( );
        double sphericalHarmonicPotentialAtJupiterFromIo = ioGravityField->getGravitationalPotential(
                sphericalHarmonicGravityOnJupiterFromIo->getCurrentRelativePosition( ), expansionDegree, expansionDegree );

        // Create mutual spherical harmonic gravity between Io and Jupiter on Io, Jupiter fixed (mu = Io + Jupiter)
        std::shared_ptr< AccelerationSettings > mutualDirectJupiterIoShGravitySettings =
                std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                        expansionDegree, expansionDegree, expansionDegree, expansionDegree );
        std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualDirectJupiterIoShGravity =
                std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >( createAccelerationModel(
                        bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), mutualDirectJupiterIoShGravitySettings, "Io", "Jupiter" ) );

        mutualDirectJupiterIoShGravity->updateMembers( );
        Eigen::Vector3d mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration = mutualDirectJupiterIoShGravity->getAcceleration( );

        {
            // Full-two-body settings [l1<=N, l2=0] should match the one-way Io-on-Jupiter SH model after
            // converting accelerated-body convention with the gravitational-parameter scale factor.
            std::shared_ptr< AccelerationSettings > extendedBodySettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( expansionDegree, expansionDegree, 0, 0 );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d precomputedAcceleration = ( ( -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) ) /
                                                        5.959916033410404E012 * sphericalHarmonicGravityOnJupiterFromIoAcceleration );

            Eigen::Vector3d accelerationDifference = precomputedAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Full-two-body settings [l1=0, l2<=N] should match Jupiter-on-Io one-way SH acceleration directly.
            std::shared_ptr< AccelerationSettings > extendedBodySettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 0, 0, expansionDegree, expansionDegree );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );
            ;

            Eigen::Vector3d accelerationDifference = sphericalHarmonicGravityOnIoFromJupiterAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Same [l1=0, l2<=N] check, now swapping accelerated/accelerating bodies.
            std::shared_ptr< AccelerationSettings > extendedBodySettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 0, 0, expansionDegree, expansionDegree );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d accelerationDifference = sphericalHarmonicGravityOnJupiterFromIoAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Same [l1<=N, l2=0] check with swapped body roles and corresponding scale/sign conversion.
            std::shared_ptr< AccelerationSettings > extendedBodySettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( expansionDegree, expansionDegree, 0, 0 );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d precomputedAcceleration = -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) /
                    5.959916033410404E012 * sh1ExtendedGravityAcceleration;

            Eigen::Vector3d accelerationDifference = sphericalHarmonicGravityOnIoFromJupiterAcceleration - precomputedAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Full-two-body settings restricted to extended-single-point interactions should reproduce
            // the legacy mutual SH acceleration model.
            std::shared_ptr< AccelerationSettings > extendedBodySettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >(
                            getExtendedSinglePointMassInteractions( expansionDegree, expansionDegree, expansionDegree, expansionDegree ) );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d accelerationDifference =
                    mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Same extended-single-point equivalence check for swapped body roles.
            std::shared_ptr< AccelerationSettings > extendedBodySettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >(
                            getExtendedSinglePointMassInteractions( expansionDegree, expansionDegree, expansionDegree, expansionDegree ) );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d precomputedAcceleration = -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) /
                    5.959916033410404E012 * sh1ExtendedGravityAcceleration;

            Eigen::Vector3d accelerationDifference = mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration - precomputedAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Third-body mutual-extended consistency:
            // a_third-body = a(undergoing wrt exerting) - a(central wrt exerting),
            // where each direct term comes from a matching full-two-body model.
            const int maximumDegreeOfBodyUndergoingAcceleration = 4;
            const int maximumDegreeOfBodyExertingAcceleration = 3;
            const int maximumDegreeOfCentralBody = 6;

            std::shared_ptr< AccelerationSettings > mutualThirdBodyExtendedSettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( maximumDegreeOfBodyUndergoingAcceleration,
                                                                                          maximumDegreeOfBodyUndergoingAcceleration,
                                                                                          maximumDegreeOfBodyExertingAcceleration,
                                                                                          maximumDegreeOfBodyExertingAcceleration,
                                                                                          maximumDegreeOfCentralBody,
                                                                                          maximumDegreeOfCentralBody );

            std::shared_ptr< ThirdBodyFullTwoBodySphericalHarmonicsGravitationalAccelerationModel > thirdBodyMutualExtendedModel =
                    std::dynamic_pointer_cast< ThirdBodyFullTwoBodySphericalHarmonicsGravitationalAccelerationModel >(
                            createAccelerationModel( bodyMap.at( "Europa" ),
                                                     bodyMap.at( "Io" ),
                                                     mutualThirdBodyExtendedSettings,
                                                     "Europa",
                                                     "Io",
                                                     bodyMap.at( "Jupiter" ),
                                                     "Jupiter" ) );
            BOOST_REQUIRE( thirdBodyMutualExtendedModel != nullptr );
            thirdBodyMutualExtendedModel->updateMembers( );
            const Eigen::Vector3d thirdBodyMutualExtendedAcceleration = thirdBodyMutualExtendedModel->getAcceleration( );

            std::shared_ptr< AccelerationSettings > directUndergoingSettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( maximumDegreeOfBodyUndergoingAcceleration,
                                                                                          maximumDegreeOfBodyUndergoingAcceleration,
                                                                                          maximumDegreeOfBodyExertingAcceleration,
                                                                                          maximumDegreeOfBodyExertingAcceleration );

            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > directAccelerationOnUndergoingBody =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >( createAccelerationModel(
                            bodyMap.at( "Europa" ), bodyMap.at( "Io" ), directUndergoingSettings, "Europa", "Io" ) );
            BOOST_REQUIRE( directAccelerationOnUndergoingBody != nullptr );
            directAccelerationOnUndergoingBody->updateMembers( );
            const Eigen::Vector3d directUndergoingAcceleration = directAccelerationOnUndergoingBody->getAcceleration( );

            std::shared_ptr< AccelerationSettings > directCentralSettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( maximumDegreeOfCentralBody,
                                                                                          maximumDegreeOfCentralBody,
                                                                                          maximumDegreeOfBodyExertingAcceleration,
                                                                                          maximumDegreeOfBodyExertingAcceleration );

            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > directAccelerationOnCentralBody =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >( createAccelerationModel(
                            bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), directCentralSettings, "Jupiter", "Io" ) );
            BOOST_REQUIRE( directAccelerationOnCentralBody != nullptr );
            directAccelerationOnCentralBody->updateMembers( );
            const Eigen::Vector3d directCentralAcceleration = directAccelerationOnCentralBody->getAcceleration( );

            const Eigen::Vector3d directAccelerationFromThirdBodyModel =
                    thirdBodyMutualExtendedModel->getAccelerationModelForBodyUndergoingAcceleration( )->getAcceleration( );
            const Eigen::Vector3d centralAccelerationFromThirdBodyModel =
                    thirdBodyMutualExtendedModel->getAccelerationModelForCentralBody( )->getAcceleration( );

            BOOST_CHECK_GT( directAccelerationFromThirdBodyModel.norm( ), 1.0E-16 );
            BOOST_CHECK_GT( centralAccelerationFromThirdBodyModel.norm( ), 1.0E-16 );
            BOOST_CHECK_GT( thirdBodyMutualExtendedAcceleration.norm( ), 1.0E-16 );

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL(
                        directAccelerationFromThirdBodyModel( i ) - directUndergoingAcceleration( i ),
                        15.0 * std::numeric_limits< double >::epsilon( ) * std::max( directAccelerationFromThirdBodyModel.norm( ), 1.0 ) );

                BOOST_CHECK_SMALL(
                        centralAccelerationFromThirdBodyModel( i ) - directCentralAcceleration( i ),
                        15.0 * std::numeric_limits< double >::epsilon( ) * std::max( centralAccelerationFromThirdBodyModel.norm( ), 1.0 ) );

                BOOST_CHECK_SMALL(
                        thirdBodyMutualExtendedAcceleration( i ) -
                                ( directAccelerationFromThirdBodyModel( i ) - centralAccelerationFromThirdBodyModel( i ) ),
                        15.0 * std::numeric_limits< double >::epsilon( ) * std::max( thirdBodyMutualExtendedAcceleration.norm( ), 1.0 ) );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testDegreeTwoCrossTermMutualAccelerationIsolation )
{
    // Test rationale:
    // Validate explicit isolation of degree-2/degree-2 figure-figure acceleration terms by two independent paths:
    // 1) subtraction from the full degree-2 interaction model,
    // 2) direct model containing only (l1=2,l2=2) combinations.
    //
    // Additional linearity/zero-response checks verify expected dependence on body-2 degree-2 coefficients.
    const double gravitationalParameter = 5.0E5;
    const double equatorialRadiusOfBody1 = 1300.0;
    const double equatorialRadiusOfBody2 = 1100.0;

    Eigen::MatrixXd cosineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd cosineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );

    cosineCoefficientsOfBody1( 0, 0 ) = 1.0;
    cosineCoefficientsOfBody2( 0, 0 ) = 1.0;

    cosineCoefficientsOfBody1( 2, 0 ) = 0.23;
    cosineCoefficientsOfBody1( 2, 1 ) = -0.11;
    sineCoefficientsOfBody1( 2, 1 ) = 0.14;
    cosineCoefficientsOfBody1( 2, 2 ) = 0.09;
    sineCoefficientsOfBody1( 2, 2 ) = -0.06;

    cosineCoefficientsOfBody2( 2, 0 ) = -0.19;
    cosineCoefficientsOfBody2( 2, 1 ) = 0.16;
    sineCoefficientsOfBody2( 2, 1 ) = -0.08;
    cosineCoefficientsOfBody2( 2, 2 ) = 0.13;
    sineCoefficientsOfBody2( 2, 2 ) = 0.12;

    runDegreeTwoCrossTermValidationCase(
            Eigen::Vector3d( 5100.0, -2200.0, 3600.0 ),
            Eigen::Quaterniond( Eigen::AngleAxisd( 0.7, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( -0.4, Eigen::Vector3d::UnitX( ) ) *
                                Eigen::AngleAxisd( 0.2, Eigen::Vector3d::UnitY( ) ) ),
            Eigen::Quaterniond( Eigen::AngleAxisd( -0.3, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 0.5, Eigen::Vector3d::UnitY( ) ) *
                                Eigen::AngleAxisd( 0.25, Eigen::Vector3d::UnitX( ) ) ),
            cosineCoefficientsOfBody1,
            sineCoefficientsOfBody1,
            cosineCoefficientsOfBody2,
            sineCoefficientsOfBody2,
            1234.0,
            gravitationalParameter,
            equatorialRadiusOfBody1,
            equatorialRadiusOfBody2 );

    cosineCoefficientsOfBody1( 2, 0 ) = -0.17;
    cosineCoefficientsOfBody1( 2, 1 ) = 0.13;
    sineCoefficientsOfBody1( 2, 1 ) = -0.07;
    cosineCoefficientsOfBody1( 2, 2 ) = 0.04;
    sineCoefficientsOfBody1( 2, 2 ) = 0.15;

    cosineCoefficientsOfBody2( 2, 0 ) = 0.12;
    cosineCoefficientsOfBody2( 2, 1 ) = -0.18;
    sineCoefficientsOfBody2( 2, 1 ) = 0.09;
    cosineCoefficientsOfBody2( 2, 2 ) = -0.11;
    sineCoefficientsOfBody2( 2, 2 ) = 0.05;

    runDegreeTwoCrossTermValidationCase( Eigen::Vector3d( -4300.0, 3100.0, 2800.0 ),
                                         Eigen::Quaterniond( Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitX( ) ) *
                                                             Eigen::AngleAxisd( 0.45, Eigen::Vector3d::UnitY( ) ) *
                                                             Eigen::AngleAxisd( 1.0, Eigen::Vector3d::UnitZ( ) ) ),
                                         Eigen::Quaterniond( Eigen::AngleAxisd( 0.35, Eigen::Vector3d::UnitX( ) ) *
                                                             Eigen::AngleAxisd( -0.6, Eigen::Vector3d::UnitZ( ) ) *
                                                             Eigen::AngleAxisd( 0.15, Eigen::Vector3d::UnitY( ) ) ),
                                         cosineCoefficientsOfBody1,
                                         sineCoefficientsOfBody1,
                                         cosineCoefficientsOfBody2,
                                         sineCoefficientsOfBody2,
                                         5678.0,
                                         gravitationalParameter,
                                         equatorialRadiusOfBody1,
                                         equatorialRadiusOfBody2 );
}

BOOST_AUTO_TEST_CASE( testFullTwoBodyAccelerationAsMutualPotentialGradient )
{
    const auto coefficientCombinations = FullTwoBodySphericalHarmonicAccelerationSettings( 4, 2, 3, 2 ).coefficientCombinationsToUse_;
    const double gravitationalParameter = 5.0E5, radius1 = 1300.0, radius2 = 900.0, currentTime = 42.0;
    const Eigen::Vector3d nominalPosition( 5100.0, -2300.0, 3700.0 );
    Eigen::MatrixXd cosine1 = Eigen::MatrixXd::Zero( 5, 3 ), sine1 = Eigen::MatrixXd::Zero( 5, 3 );
    Eigen::MatrixXd cosine2 = Eigen::MatrixXd::Zero( 4, 3 ), sine2 = Eigen::MatrixXd::Zero( 4, 3 );
    auto fillCoefficients = []( Eigen::MatrixXd& cosineCoefficients, Eigen::MatrixXd& sineCoefficients, const double scale ) {
        cosineCoefficients( 0, 0 ) = 1.0;
        for( int degree = 1; degree < cosineCoefficients.rows( ); degree++ )
        {
            for( int order = 0; order <= std::min( degree, static_cast< int >( cosineCoefficients.cols( ) ) - 1 ); order++ )
            {
                cosineCoefficients( degree, order ) = scale * ( 0.3 + 0.7 * degree - 0.2 * order );
                if( order > 0 )
                {
                    sineCoefficients( degree, order ) = scale * ( -0.4 * degree + 0.6 * order );
                }
            }
        }
    };
    fillCoefficients( cosine1, sine1, 2.0E-2 );
    fillCoefficients( cosine2, sine2, -1.5E-2 );

    const Eigen::Quaterniond rotationToBody1( Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitZ( ) ) *
                                              Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitX( ) ) );
    const Eigen::Quaterniond rotationToBody2( Eigen::AngleAxisd( -0.3, Eigen::Vector3d::UnitY( ) ) *
                                              Eigen::AngleAxisd( 0.5, Eigen::Vector3d::UnitZ( ) ) );
    auto model = createMutualExtendedBodyAccelerationModel( coefficientCombinations,
                                                            nominalPosition,
                                                            Eigen::Vector3d::Zero( ),
                                                            gravitationalParameter,
                                                            radius1,
                                                            radius2,
                                                            cosine1,
                                                            sine1,
                                                            cosine2,
                                                            sine2,
                                                            rotationToBody1,
                                                            rotationToBody2 );
    model->updateMembers( currentTime );
    const Eigen::Vector3d bodyFixedPosition = model->getCurrentBodyFixedRelativePosition( );
    auto getPotential = [ & ]( const Eigen::Vector3d& position ) {
        return model->getEffectiveMutualPotentialField( )->getGravitationalPotential( position, model->getSphericalHarmonicsCache( ) );
    };
    const double perturbation = 1.0E-1;
    Eigen::Vector3d finiteDifferenceGradientInBody1Frame;
    for( int i = 0; i < 3; i++ )
    {
        finiteDifferenceGradientInBody1Frame( i ) = ( getPotential( bodyFixedPosition + perturbation * Eigen::Vector3d::Unit( i ) ) -
                                                      getPotential( bodyFixedPosition - perturbation * Eigen::Vector3d::Unit( i ) ) ) /
                ( 2.0 * perturbation );
    }
    const Eigen::Vector3d finiteDifferenceGradient =
            model->getCurrentRotationFromInertialToBody1( ).inverse( ) * finiteDifferenceGradientInBody1Frame;
    BOOST_CHECK_SMALL( ( model->getAcceleration( ) - finiteDifferenceGradient ).norm( ) / model->getAcceleration( ).norm( ), 1.0E-8 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
