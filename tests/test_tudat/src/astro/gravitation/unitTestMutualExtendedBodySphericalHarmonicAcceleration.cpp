#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"


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
        boost::shared_ptr< RandomVariableGenerator< double > > randomNumberGenerator,
        const int maximumDegree, const int maximumOrder )
{
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );

//    for( int i = 0; i <= maximumDegree; i++ )
//    {
//        for( int j = 0; ( ( j <= i ) && ( j <= maximumOrder ) ); j++ )
//        {
//          cosineCoefficients( i, j ) = randomNumberGenerator->getRandomVariableValue( );
//            if( j != 0 )
//            {
////                sineCoefficients( i, j ) = randomNumberGenerator->getRandomVariableValue( );
//            }
//        }
//    }

//    sineCoefficients( 2, 1 ) = randomNumberGenerator->getRandomVariableValue( );
    sineCoefficients( 2, 1 ) = randomNumberGenerator->getRandomVariableValue( );

    cosineCoefficients( 0, 0 ) = 1.0;


    return std::make_pair( cosineCoefficients, sineCoefficients );
}

boost::shared_ptr< tudat::simulation_setup::GravityFieldSettings > getDummyJovianSystemGravityField(
        const std::string& bodyName,
        const int isNormalized )
{
    boost::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    std::vector< double > randomNumberSettings;
    randomNumberSettings.push_back( 0.0 );
    randomNumberSettings.push_back( 1.0E-4 );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;

    if( bodyName == "Jupiter" )
    {
        boost::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator = createBoostContinuousRandomVariableGenerator(
                    normal_boost_distribution, randomNumberSettings , 0.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ),
                  coefficients.first, coefficients.second, "IAU_Jupiter_SIMPLIFIED" );
    }
    else if( bodyName == "Io" )
    {
        boost::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator = createBoostContinuousRandomVariableGenerator(
                    normal_boost_distribution, randomNumberSettings , 1.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( 5.959916033410404E012, getAverageRadius( "Io" ),
                  coefficients.first, coefficients.second, "IAU_Io_SIMPLIFIED" );
    }
    else if( bodyName == "Europa" )
    {
        boost::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator = createBoostContinuousRandomVariableGenerator(
                    normal_boost_distribution, randomNumberSettings , 2.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( 3.202738774922892E12, getAverageRadius( "Europa" ),
                  coefficients.first, coefficients.second, "IAU_Europa_SIMPLIFIED" );
    }

    return gravityFieldSettings;

}



//void getBodyCoefficientEulerAnglePartials(
//        const Eigen::Quaterniond nominalQuaternion,
//        boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualShField,
//        std::vector< Eigen::MatrixXd >& numericalTransformedCosineCoefficientsOfBody2Partials,
//        std::vector< Eigen::MatrixXd >& numericalTransformedSineCoefficientsOfBody2Partials,
//        const double perturbation )
//{
//    numericalTransformedCosineCoefficientsOfBody2Partials.resize( 3 );
//    numericalTransformedSineCoefficientsOfBody2Partials.resize( 3 );


//    Eigen::Vector3d nominalEulerAngles = basic_mathematics::get313EulerAnglesFromQuaternion( nominalQuaternion );
//    Eigen::Vector3d perturbedEulerAngles;
//    Eigen::MatrixXd upPerturbedCosineCoefficients, upPerturbedSineCoefficients;
//    Eigen::MatrixXd downPerturbedCosineCoefficients, downPerturbedSineCoefficients;

//    for( unsigned int i = 0; i < 3; i ++ )
//    {
//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) += perturbation;
//        mutualShField->computeCurrentEffectiveCoefficients( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0 ) );
//        upPerturbedCosineCoefficients = mutualShField->getTransformedCosineCoefficientsOfBody2( );
//        upPerturbedSineCoefficients = mutualShField->getTransformedSineCoefficientsOfBody2( );


//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) -= perturbation;
//        mutualShField->computeCurrentEffectiveCoefficients( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0 ) );
//        downPerturbedCosineCoefficients = mutualShField->getTransformedCosineCoefficientsOfBody2( );
//        downPerturbedSineCoefficients = mutualShField->getTransformedSineCoefficientsOfBody2( );

//        numericalTransformedCosineCoefficientsOfBody2Partials[ i ] =
//                ( upPerturbedCosineCoefficients - downPerturbedCosineCoefficients ) / ( 2.0 * perturbation );
//        numericalTransformedSineCoefficientsOfBody2Partials[ i ] =
//                ( upPerturbedSineCoefficients - downPerturbedSineCoefficients ) / ( 2.0 * perturbation );

//    }
//}

Eigen::Matrix2d getEffectiveMutualPotentialCoefficientWrtBody2Coefficient(
        const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2,
        const Eigen::MatrixXd& nominalTransformedCosineCoefficients,
        const Eigen::MatrixXd& nominalTransformedSineCoefficients,
        const boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualPotentialField,
        const double coefficientPerturbation )
{
    Eigen::MatrixXd perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    Eigen::MatrixXd perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;


    perturbedTransformedCosineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) += coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double upPerturbedCosineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double upPerturbedSineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );



    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedCosineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) -= coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double downPerturbedCosineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double downPerturbedSineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );



    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedSineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) += coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double upPerturbedCosineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double upPerturbedSineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );



    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedSineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) -= coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double downPerturbedCosineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double downPerturbedSineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

    return ( ( Eigen::Matrix2d( )<<
               ( upPerturbedCosineCoefficientFromCosinePerturbation  - downPerturbedCosineCoefficientFromCosinePerturbation ),
               ( upPerturbedCosineCoefficientFromSinePerturbation  - downPerturbedCosineCoefficientFromSinePerturbation ),
               ( upPerturbedSineCoefficientFromCosinePerturbation  - downPerturbedSineCoefficientFromCosinePerturbation ),
               ( upPerturbedSineCoefficientFromSinePerturbation  - downPerturbedSineCoefficientFromSinePerturbation ) ).finished( ) ) /
            ( 2.0 * coefficientPerturbation );


}

BOOST_AUTO_TEST_SUITE( test_mutual_spherical_harmonic_gravity )

BOOST_AUTO_TEST_CASE( testMutualSphericalHarmonicGravity )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_jup310_small.bsp" } );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialTime = 1.0E7;
    double finalTime = 1.2E7;

    int expanstionDegree = 3;

    //for( int isNormalized = 1; isNormalized < 1; isNormalized++ )
    int isNormalized = 1;
    {
        // Get body settings.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodyNames, initialTime, finalTime );
        bodySettings[ "Jupiter" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Jupiter", isNormalized );
        bodySettings[ "Io" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Io", isNormalized );
        bodySettings[ "Europa" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Europa", isNormalized );

        bodySettings[ "Jupiter" ]->rotationModelSettings = boost::make_shared<
                SimpleRotationModelSettings >( "ECLIPJ2000", "IAU_Jupiter_SIMPLIFIED", Eigen::Quaterniond(
                                                   Eigen::AngleAxisd( -1.4, Eigen::Vector3d::UnitZ( ) ) *
                                                   Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitX( ) ) *
                                                   Eigen::AngleAxisd( 2.6, Eigen::Vector3d::UnitZ( ) ) ), 0.0, 0.0 );

        bodySettings[ "Io" ]->rotationModelSettings = boost::make_shared<
                SimpleRotationModelSettings >( "ECLIPJ2000", "IAU_Io_SIMPLIFIED", Eigen::Quaterniond(
                                                   Eigen::AngleAxisd( mathematical_constants::PI * 3.0 / 8.0, Eigen::Vector3d::UnitZ( ) ) *
                                                   Eigen::AngleAxisd( mathematical_constants::PI * 5.0 / 8.0, Eigen::Vector3d::UnitX( ) ) *
                                                   Eigen::AngleAxisd( mathematical_constants::PI * -2.5 / 8.0, Eigen::Vector3d::UnitZ( ) ) ), 0.0, 0.0 );
        bodySettings[ "Europa" ]->rotationModelSettings = boost::make_shared<
                SimpleRotationModelSettings >( "ECLIPJ2000", "IAU_Europa_SIMPLIFIED", Eigen::Quaterniond(
                                                   Eigen::AngleAxisd( mathematical_constants::PI * 3.0 / 8.0, Eigen::Vector3d::UnitZ( ) ) ), 0.0, 0.0 );

        // Create bodies needed in simulation
        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Set current state and rotation of bodies.
        double currentTime = 1.1E7;
        bodyMap[ "Jupiter" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap[ "Jupiter" ]->setStateFromEphemeris( currentTime );
        bodyMap[ "Io" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap[ "Io" ]->setStateFromEphemeris( currentTime );
        bodyMap[ "Europa" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap[ "Europa" ]->setStateFromEphemeris( currentTime );

        // Retrieve gravity fields.
        boost::shared_ptr< SphericalHarmonicsGravityField > jupiterGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    ( bodyMap.at( "Jupiter" ) )->getGravityFieldModel( ) );
        boost::shared_ptr< SphericalHarmonicsGravityField > ioGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    ( bodyMap.at( "Io" ) )->getGravityFieldModel( ) );
        boost::shared_ptr< SphericalHarmonicsGravityField > europaGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    ( bodyMap.at( "Europa" ) )->getGravityFieldModel( ) );

        // Create central gravity acceleration (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > centralGravitySettings = boost::make_shared< AccelerationSettings >(
                    central_gravity );
        boost::shared_ptr< CentralGravitationalAccelerationModel3d > centralGravity =
                boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), centralGravitySettings, "Io", "Jupiter",
                                             bodyMap.at( "Jupiter" ), "Jupiter" ) );

        // Calculate central gravity acceleration.
        centralGravity->updateMembers( );
        Eigen::Vector3d centralGravityAcceleration = centralGravity->getAcceleration( );

        // Create spherical harmonic gravity of Jupiter on Io, Jupiter-fixed (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnIoFromJupiterSettings =
                boost::make_shared< SphericalHarmonicAccelerationSettings >( expanstionDegree, expanstionDegree );
        boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnIoFromJupiter =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel(  bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), sphericalHarmonicGravityOnIoFromJupiterSettings,
                                              "Io", "Jupiter" ) );

        // Calculate spherical harmonic gravity of Jupiter on Io.
        sphericalHarmonicGravityOnIoFromJupiter->updateMembers( );
        Eigen::Vector3d sphericalHarmonicGravityOnIoFromJupiterAcceleration = sphericalHarmonicGravityOnIoFromJupiter->getAcceleration( );
        double sphericalHarmonicPotentialAtIoFromJupiter = jupiterGravityField->getGravitationalPotential(
                    sphericalHarmonicGravityOnIoFromJupiter->getCurrentRelativePosition( ), expanstionDegree, expanstionDegree );

        // Create spherical harmonic gravity of Io on Jupiter, Io-fixed (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnJupiterFromIoSettings =
                boost::make_shared< SphericalHarmonicAccelerationSettings >( expanstionDegree, expanstionDegree );
        boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnJupiterFromIo =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), sphericalHarmonicGravityOnJupiterFromIoSettings,
                                             "Jupiter", "Io" ) );

        // Calculate spherical harmonic gravity of Io on Jupiter.
        sphericalHarmonicGravityOnJupiterFromIo->updateMembers( );
        Eigen::Vector3d sphericalHarmonicGravityOnJupiterFromIoAcceleration = sphericalHarmonicGravityOnJupiterFromIo->getAcceleration( );
        double sphericalHarmonicPotentialAtJupiterFromIo = ioGravityField->getGravitationalPotential(
                    sphericalHarmonicGravityOnJupiterFromIo->getCurrentRelativePosition( ), expanstionDegree, expanstionDegree );

        // Create mutual spherical harmonic gravity between Io and Jupiter on Io, Jupiter fixed (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > mutualDirectJupiterIoShGravitySettings =
                boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( expanstionDegree, expanstionDegree, expanstionDegree, expanstionDegree );
        boost::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualDirectJupiterIoShGravity =
                boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), mutualDirectJupiterIoShGravitySettings,
                                             "Io", "Jupiter" ) );

        mutualDirectJupiterIoShGravity->updateMembers( );
        Eigen::Vector3d mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration = mutualDirectJupiterIoShGravity->getAcceleration( );


        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        expanstionDegree, expanstionDegree, 0, 0 );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d relativePosition = sphericalHarmonicGravityOnJupiterFromIo->getCurrentInertialRelativePosition( );
            Eigen::Vector3d pointMassAcceleration =
                    spice_interface::getBodyGravitationalParameter( "Jupiter" ) * relativePosition /
                    ( std::pow( relativePosition.norm( ), 3 ) );

            std::cout<<"Extended body:  "<<( sh1ExtendedGravityAcceleration - pointMassAcceleration ).transpose( )<<std::endl;

            Eigen::Vector3d precomputedAcceleration = ( ( -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) ) /
                                                        5.959916033410404E012 * sphericalHarmonicGravityOnJupiterFromIoAcceleration );

            std::cout<<"Original:  "<<( precomputedAcceleration - pointMassAcceleration ).transpose( )<<std::endl;

            std::cout<<"Difference: "<<( ( precomputedAcceleration - sh1ExtendedGravityAcceleration ) ).transpose( )<<std::endl;

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( precomputedAcceleration, sh1ExtendedGravityAcceleration, 5.0E-14 );

        }

        sleep( 10000.0 );

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        0, 0, expanstionDegree, expanstionDegree );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sphericalHarmonicGravityOnIoFromJupiterAcceleration, sh1ExtendedGravityAcceleration, 1.0E-14 );

        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        0, 0, expanstionDegree, expanstionDegree );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sphericalHarmonicGravityOnJupiterFromIoAcceleration, sh1ExtendedGravityAcceleration, 5.0E-14 );
        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        expanstionDegree, expanstionDegree, 0, 0 );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( sphericalHarmonicGravityOnIoFromJupiterAcceleration ),
                        ( -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) / 5.959916033410404E012 * sh1ExtendedGravityAcceleration ), 1.0E-14 );
        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        getExtendedSinglePointMassInteractions( expanstionDegree, expanstionDegree, expanstionDegree, expanstionDegree ) );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ),  extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration, sh1ExtendedGravityAcceleration, 5.0E-14 );
        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        getExtendedSinglePointMassInteractions( expanstionDegree, expanstionDegree, expanstionDegree, expanstionDegree ) );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ),  extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration,
                        ( -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) / 5.959916033410404E012 * sh1ExtendedGravityAcceleration ), 5.0E-14 );
        }

        {
            std::cout<<"Test start"<<std::endl;
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        expanstionDegree, expanstionDegree, expanstionDegree, expanstionDegree );
            //getExtendedSinglePointMassInteractions( expanstionDegree, expanstionDegree, expanstionDegree, expanstionDegree ) );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ),  extendedBodySettings, "Io", "Jupiter" ) );

//            boost::shared_ptr< MutualExtendedBodySphericalHarmonicTorque > sh1ExtendedTorque =
//                    boost::make_shared< MutualExtendedBodySphericalHarmonicTorque >(
//                        sh1ExtendedGravity );

            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

//            sh1ExtendedTorque->updateMembers( );

//            {
//                std::vector< Eigen::MatrixXd > transformedCosineCoefficientsOfBody2Partials;
//                std::vector< Eigen::MatrixXd > transformedSineCoefficientsOfBody2Partials;

////                sh1ExtendedTorque->getBodyCoefficientPartialsWrtEulerAngles(
////                            transformedCosineCoefficientsOfBody2Partials,
////                            transformedSineCoefficientsOfBody2Partials );

//                std::vector< Eigen::MatrixXd > numericalTransformedCosineCoefficientsOfBody2Partials;
//                std::vector< Eigen::MatrixXd > numericalTransformedSineCoefficientsOfBody2Partials;


//                getBodyCoefficientEulerAnglePartials(
//                            sh1ExtendedGravity->getCurrentRotationFromBody2ToBody1( ), sh1ExtendedGravity->getEffectiveMutualPotentialField( ),
//                            numericalTransformedCosineCoefficientsOfBody2Partials,
//                            numericalTransformedSineCoefficientsOfBody2Partials,
//                            std::pow( 10, -5 ) );

//                for( unsigned int i = 0; i < 3; i++ )
//                {
//                    for( unsigned int j = 0; j <= 5; j++ )
//                    {
//                        for( unsigned int k = 0; k <= 5; k++ )
//                        {
//                            BOOST_CHECK_CLOSE_FRACTION(
//                                        numericalTransformedCosineCoefficientsOfBody2Partials.at( i )( j, k ),
//                                        transformedCosineCoefficientsOfBody2Partials.at( i )( j, k ), 2.0E-7 );
//                            BOOST_CHECK_CLOSE_FRACTION(
//                                        numericalTransformedSineCoefficientsOfBody2Partials.at( i )( j, k ),
//                                        transformedSineCoefficientsOfBody2Partials.at( i )( j, k ), 2.0E-7 );
//                        }
//                    }
//                }
//            }

//            {
//                std::map< int, std::map< int, std::map< int, std::map< int, Eigen::Matrix2d > > > > fullCoefficientsWrtBody2Coefficients;
//                sh1ExtendedTorque->getFullCoefficientsWrtBody2Coefficients( fullCoefficientsWrtBody2Coefficients );
//                std::map< int, std::map< int, std::map< int, std::map< int, Eigen::Matrix2d > > > > numericalFullCoefficientsWrtBody2Coefficients;

//                boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualPotentialField =
//                        sh1ExtendedGravity->getEffectiveMutualPotentialField( );
//                std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse =
//                        mutualPotentialField->getCoefficientCombinationsToUse( );

//                int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2;

//                Eigen::MatrixXd nominalTransformedCosineCoefficients = mutualPotentialField->getTransformedCosineCoefficientsOfBody2( );
//                Eigen::MatrixXd nominalTransformedSineCoefficients = mutualPotentialField->getTransformedCosineCoefficientsOfBody2( );

//                double coefficientPerturbation = 1.0E-10;

//                for( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
//                {
//                    degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
//                    orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
//                    degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
//                    orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );

//                    numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ orderOfBody1 ][ degreeOfBody2 ][ orderOfBody2 ] =
//                            getEffectiveMutualPotentialCoefficientWrtBody2Coefficient(
//                                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2,
//                                nominalTransformedCosineCoefficients,
//                                nominalTransformedSineCoefficients,
//                                mutualPotentialField,
//                                coefficientPerturbation );
//                    numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ -orderOfBody1 ][ degreeOfBody2 ][ orderOfBody2 ] =
//                            getEffectiveMutualPotentialCoefficientWrtBody2Coefficient(
//                                degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2,
//                                nominalTransformedCosineCoefficients,
//                                nominalTransformedSineCoefficients,
//                                mutualPotentialField,
//                                coefficientPerturbation );
//                    numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ orderOfBody1 ][ degreeOfBody2 ][ -orderOfBody2 ] =
//                            getEffectiveMutualPotentialCoefficientWrtBody2Coefficient(
//                                degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2,
//                                nominalTransformedCosineCoefficients,
//                                nominalTransformedSineCoefficients,
//                                mutualPotentialField,
//                                coefficientPerturbation );
//                    numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ -orderOfBody1 ][ degreeOfBody2 ][ -orderOfBody2 ] =
//                            getEffectiveMutualPotentialCoefficientWrtBody2Coefficient(
//                                degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2,
//                                nominalTransformedCosineCoefficients,
//                                nominalTransformedSineCoefficients,
//                                mutualPotentialField,
//                                coefficientPerturbation );

//                    //std::cout<<degreeOfBody1<<" "<<orderOfBody1<<" "<<degreeOfBody2<<" "<<orderOfBody2 <<std::endl;
//                    double tolerance = 1.0E-5;
//                    if( degreeOfBody1 == 2 && orderOfBody1 == 1 )
//                    {
//                        tolerance *= 100.0;
//                    }

//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                                numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ orderOfBody1 ][ degreeOfBody2 ][ orderOfBody2 ],
//                                fullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ orderOfBody1 ][ degreeOfBody2 ][ orderOfBody2 ], tolerance );

//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                                numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ -orderOfBody1 ][ degreeOfBody2 ][ orderOfBody2 ],
//                            fullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ -orderOfBody1 ][ degreeOfBody2 ][ orderOfBody2 ], tolerance );

//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                                numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ orderOfBody1 ][ degreeOfBody2 ][ -orderOfBody2 ],
//                            fullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ orderOfBody1 ][ degreeOfBody2 ][ -orderOfBody2 ], tolerance );

//                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                                numericalFullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ -orderOfBody1 ][ degreeOfBody2 ][ -orderOfBody2 ],
//                            fullCoefficientsWrtBody2Coefficients[ degreeOfBody1 ][ -orderOfBody1 ][ degreeOfBody2 ][ -orderOfBody2 ], tolerance );
//                }
//            }

//            {
//                std::map< int, std::map< int, std::map< int, std::map< int, Eigen::Matrix< double, 1, 2 > > > > > potentialComponentsWrtFullCoefficients;
//                sh1ExtendedTorque->getPotentialComponentsWrtFullCoefficients( potentialComponentsWrtFullCoefficients );

//                boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualPotentialField =
//                        sh1ExtendedGravity->getEffectiveMutualPotentialField( );
//                std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse =
//                        mutualPotentialField->getCoefficientCombinationsToUse( );

//                int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2;

//                for( unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
//                {
//                    degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
//                    orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
//                    degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
//                    orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );
//                }

//            }
        }
    }
}
BOOST_AUTO_TEST_SUITE_END( )

}

}

