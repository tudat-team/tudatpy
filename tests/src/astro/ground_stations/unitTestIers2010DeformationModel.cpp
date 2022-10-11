#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iomanip>

#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"

namespace tudat
{
namespace unit_tests
{


using namespace tudat;
using namespace tudat::ephemerides;
using namespace tudat::utilities;
using namespace tudat::spice_interface;
using namespace tudat::basic_astrodynamics;
using namespace tudat::basic_astrodynamics::iers_2010_parameters;

Eigen::Vector6d calculateFundamentalArgumentsIersCode( const double universalTimeSinceJ2000 )
{
    int daysSinceJ2000 = std::floor( universalTimeSinceJ2000 / physical_constants::JULIAN_DAY + 0.5 );
    double FHR = ( universalTimeSinceJ2000 / physical_constants::JULIAN_DAY + 0.5 - static_cast< double >( daysSinceJ2000 ) )/24.0;

    double T = universalTimeSinceJ2000 / ( physical_constants::JULIAN_YEAR * 100.0 );

    double S = 218.31664563
            + (481267.88194
               + (-0.0014663889
                  + (0.00000185139)*T)*T)*T ;

    double TAU = basic_mathematics::computeModulo( FHR*15
                                                   + 280.4606184
                                                   + (36000.7700536
                                                      + (0.00038793
                                                         + (-0.0000000258)*T)*T)*T
                                                   + (-S), 360.0 );

    double PR = basic_mathematics::computeModulo( (1.396971278
                                                   + (0.000308889
                                                      + (0.000000021
                                                         + (0.000000007)*T)*T)*T)*T, 360.0 );

    S = basic_mathematics::computeModulo( S + PR, 360.0 );

    double H = basic_mathematics::computeModulo( 280.46645
                                                 + (36000.7697489
                                                    + (0.00030322222
                                                       + (0.000000020
                                                          + (-0.00000000654)*T)*T)*T)*T, 360.0 );

    double P = basic_mathematics::computeModulo( 83.35324312
                                                 + (4069.01363525
                                                    + (-0.01032172222
                                                       + (-0.0000124991
                                                          + (0.00000005263)*T)*T)*T)*T, 360.0 );

    double ZNS = basic_mathematics::computeModulo(  234.95544499
                                                    + (1934.13626197
                                                       + (-0.00207561111
                                                          + (-0.00000213944
                                                             + (0.00000001650)*T)*T)*T)*T, 360.0 );

    double PS = basic_mathematics::computeModulo(  282.93734098
                                                   + (1.71945766667
                                                      + (0.00045688889
                                                         + (-0.00000001778
                                                            + (-0.00000000334)*T)*T)*T)*T, 360.0 );

    Eigen::Vector6d doodsonArguments =
            ( Eigen::Vector6d( )<<TAU,S, H, P, ZNS, PS ).finished( );
    return doodsonArguments * mathematical_constants::PI / 180.0;
}

BOOST_AUTO_TEST_SUITE( test_iers_2010_solid_earth_tide_deformation )


BOOST_AUTO_TEST_CASE( test_Iers2012DeformationModel )
{


    double lunarMassRatio = 0.0123000371;
    double solarMassRatio = 332946.0482;
    std::function< double( ) > lunarMassFunction = boost::lambda::constant( lunarMassRatio );
    std::function< double( ) > solarMassFunction = boost::lambda::constant( solarMassRatio );
    std::vector< std::function< double( ) > > massFunctions;
    massFunctions.resize( 2 );
    massFunctions[ 0 ] = lunarMassFunction;
    massFunctions[ 1 ] = solarMassFunction;

    double earthEquatorialRadius = 6378136.6;

    Eigen::Vector6d earthState = Eigen::Vector6d::Zero( );
    Eigen::Vector6d siteState = Eigen::Vector6d::Zero( );
    Eigen::Vector6d sunState = Eigen::Vector6d::Zero( );
    Eigen::Vector6d moonState = Eigen::Vector6d::Zero( );

    siteState( 0 ) = 4075578.3850;
    siteState( 1 ) = 931852.890;
    siteState( 2 ) = 4801570.154;
    sunState( 0 ) = 137859926952.015;
    sunState( 1 ) = 54228127881.4350;
    sunState( 2 ) = 23509422341.6960;
    moonState( 0 ) = -179996231.920342;
    moonState( 1 ) = -312468450.131567;
    moonState( 2 ) = -169288918.592160;

    Eigen::Matrix3d orientationMatrix = Eigen::Matrix3d::Identity( );
    Eigen::Quaterniond orientationQuaternion = Eigen::Quaterniond( orientationMatrix );

    std::function< Eigen::Quaterniond( const double ) > rotationEphemeris =
            std::bind( &RotationalEphemeris::getRotationToTargetFrame,
                       std::make_shared< ConstantRotationalEphemeris >( orientationQuaternion ), std::placeholders::_1 );

    std::function< Eigen::Vector6d( const double ) > earthEphemeris = std::bind(
                &Ephemeris::getCartesianState,
                std::make_shared< ConstantEphemeris >( earthState ), std::placeholders::_1 );
    std::function< Eigen::Vector6d( const double ) > sunEphemeris = std::bind(
                &Ephemeris::getCartesianState,
                std::make_shared< ConstantEphemeris >( sunState ), std::placeholders::_1  );
    std::function< Eigen::Vector6d( const double ) > moonEphemeris = std::bind(
                &Ephemeris::getCartesianState,
                std::make_shared< ConstantEphemeris >( moonState ), std::placeholders::_1 );

    std::vector< std::function< Eigen::Vector6d( const double ) > > ephemerides;
    ephemerides.resize( 2 );
    ephemerides[ 0 ] = moonEphemeris;
    ephemerides[ 1 ] = sunEphemeris;

    double evaluationTime = 0.1059411362080767 * 100.0 * 365.25 * 86400.0;


    std::map< int, std::pair< double, double > > nominalDisplacementLoveNumbers;
    nominalDisplacementLoveNumbers[ 2 ] = std::make_pair( 0.0, 0.0 );
    nominalDisplacementLoveNumbers[ 3 ] = std::make_pair( 0.0, 0.0 );


    std::vector< double > firstStepLatitudeDependenceTerms;
    firstStepLatitudeDependenceTerms.resize( 2 );
    firstStepLatitudeDependenceTerms[ 0 ] = 0.0;
    firstStepLatitudeDependenceTerms[ 1 ] = 0.0;

    std::vector< bool > areFirstStepCorrectionsCalculated;
    areFirstStepCorrectionsCalculated.resize( 6 );
    for( int i = 0; i < 6; i++ )
    {
        areFirstStepCorrectionsCalculated[ i ] = 0;
    }

    std::vector< double > correctionLoveAndShidaNumbers;
    correctionLoveAndShidaNumbers.resize( 6 );
    for( int i = 0; i < 6; i++ )
    {
        correctionLoveAndShidaNumbers[ i ] = 0.0;
    }


    std::shared_ptr< Iers2010EarthDeformation > deformationModel;

    Eigen::Vector3d expectedDeformation, calculatedDeformation;

    // Test IERS Eqs. 7.10 component of displacement
    {
        areFirstStepCorrectionsCalculated[ 2 ] = 1;
        areFirstStepCorrectionsCalculated[ 3 ] = 1;

        correctionLoveAndShidaNumbers[ 2 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_DIURNAL_LOVE_NUMBER;
        correctionLoveAndShidaNumbers[ 3 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_DIURNAL_SHIDA_NUMBER;

        deformationModel = std::make_shared< Iers2010EarthDeformation >
                ( earthEphemeris, ephemerides, rotationEphemeris,
                  boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                  nominalDisplacementLoveNumbers ,firstStepLatitudeDependenceTerms,
                  areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers ,"","" );

        expectedDeformation << -0.2836337012840008001E-03, 0.1125342324347507444E-03, -0.2471186224343683169E-03;

        calculatedDeformation = deformationModel->calculateDisplacement(
                    evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-6 );

        areFirstStepCorrectionsCalculated[ 2 ] = 0;
        areFirstStepCorrectionsCalculated[ 3 ] = 0;

        correctionLoveAndShidaNumbers[ 2 ] = 0.0;
        correctionLoveAndShidaNumbers[ 3 ] = 0.0;
    }

    // Test IERS Eqs. 7.11 component of displacement
    {
        areFirstStepCorrectionsCalculated[ 4 ] = 1;
        areFirstStepCorrectionsCalculated[ 5 ] = 1;

        correctionLoveAndShidaNumbers[ 4 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_SEMIDIURNAL_LOVE_NUMBER;
        correctionLoveAndShidaNumbers[ 5 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_SEMIDIURNAL_SHIDA_NUMBER;

        deformationModel = std::make_shared< Iers2010EarthDeformation >
                ( earthEphemeris, ephemerides, rotationEphemeris,
                  boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                  nominalDisplacementLoveNumbers ,firstStepLatitudeDependenceTerms,
                  areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers ,"","" );

        expectedDeformation << -0.2801334805106874015E-03, 0.2939522229284325029E-04,-0.6051677912316721561E-04;
        calculatedDeformation =  deformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-6 );

        areFirstStepCorrectionsCalculated[ 4 ] = 0;
        areFirstStepCorrectionsCalculated[ 5 ] = 0;

        correctionLoveAndShidaNumbers[ 4 ] = 0.0;
        correctionLoveAndShidaNumbers[ 5 ] = 0.0;
    }

    // Test IERS Eqs. 7.8 and 7.9 component of displacement
    {
        areFirstStepCorrectionsCalculated[ 0 ] = 1;
        areFirstStepCorrectionsCalculated[ 1 ] = 1;

        correctionLoveAndShidaNumbers[ 0 ] = iers_2010_parameters::DEGREE_TWO_DIURNAL_TOROIDAL_LOVE_NUMBER;
        correctionLoveAndShidaNumbers[ 1 ] = iers_2010_parameters::DEGREE_TWO_SEMIDIURNAL_TOROIDAL_LOVE_NUMBER;

        deformationModel = std::make_shared< Iers2010EarthDeformation >
                ( earthEphemeris, ephemerides, rotationEphemeris,
                  boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                  nominalDisplacementLoveNumbers ,firstStepLatitudeDependenceTerms,
                  areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers ,"","" );

        expectedDeformation << 0.2367189532359759044E-03, 0.5181609907284959182E-03, -0.3014881422940427977E-03;
        calculatedDeformation =  deformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-6 );

        areFirstStepCorrectionsCalculated[ 0 ] = 0;
        areFirstStepCorrectionsCalculated[ 1 ] = 0;

        correctionLoveAndShidaNumbers[ 0 ] = 0.0;
        correctionLoveAndShidaNumbers[ 1 ] = 0.0;
    }


    // Test IERS Eqs. 7.12 component of displacement
    {
        deformationModel = std::make_shared< Iers2010EarthDeformation >
                ( earthEphemeris, ephemerides, rotationEphemeris,
                  boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                  nominalDisplacementLoveNumbers ,firstStepLatitudeDependenceTerms,
                  areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers,
                  tudat::paths::getEarthDeformationDataFilesPath( ) + "/diurnalDisplacementFrequencyDependence2.txt" ,"",
                  std::bind( &calculateFundamentalArgumentsIersCode, std::placeholders::_1 ) );
        expectedDeformation <<0.4193085327321284701E-02, 0.1456681241014607395E-02, 0.5123366597450316508E-02;
        calculatedDeformation =  deformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-12 );
    }

    // Test IERS Eqs. 7.13 component of displacement
    {
        deformationModel = std::make_shared< Iers2010EarthDeformation >
                ( earthEphemeris, ephemerides, rotationEphemeris,
                  boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                  nominalDisplacementLoveNumbers ,firstStepLatitudeDependenceTerms,
                  areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers, "",
                  tudat::paths::getEarthDeformationDataFilesPath( ) + "/longPeriodDisplacementFrequencyDependence.txt" ,
                  std::bind( &calculateFundamentalArgumentsIersCode, std::placeholders::_1 ) );

        expectedDeformation <<-0.9780962849562107762E-04,-0.2236349699932734273E-04,0.3561945821351565926E-03;
        calculatedDeformation =  deformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-12 );
    }

    // Test full IERS chapter 7 model (case 1)
    {
        nominalDisplacementLoveNumbers[ 2 ] = std::make_pair( PRINCIPAL_DEGREE_TWO_LOVE_NUMBER, PRINCIPAL_DEGREE_TWO_SHIDA_NUMBER );
        nominalDisplacementLoveNumbers[ 3 ] = std::make_pair( PRINCIPAL_DEGREE_THREE_LOVE_NUMBER, PRINCIPAL_DEGREE_THREE_SHIDA_NUMBER );

        firstStepLatitudeDependenceTerms[ 0 ] = DEGREE_TWO_LATITUDE_LOVE_NUMBER;
        firstStepLatitudeDependenceTerms[ 1 ] = DEGREE_TWO_LATITUDE_SHIDA_NUMBER;

        for( int i = 0; i < 6; i++ )
        {
            areFirstStepCorrectionsCalculated[ i ] = 1;
        }

        correctionLoveAndShidaNumbers[ 0 ] = iers_2010_parameters::DEGREE_TWO_DIURNAL_TOROIDAL_LOVE_NUMBER;
        correctionLoveAndShidaNumbers[ 1 ] = iers_2010_parameters::DEGREE_TWO_SEMIDIURNAL_TOROIDAL_LOVE_NUMBER;
        correctionLoveAndShidaNumbers[ 2 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_DIURNAL_LOVE_NUMBER;
        correctionLoveAndShidaNumbers[ 3 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_DIURNAL_SHIDA_NUMBER;
        correctionLoveAndShidaNumbers[ 4 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_SEMIDIURNAL_LOVE_NUMBER;
        correctionLoveAndShidaNumbers[ 5 ] = iers_2010_parameters::IMAGINARY_DEGREE_TWO_SEMIDIURNAL_SHIDA_NUMBER;

        deformationModel = std::make_shared< Iers2010EarthDeformation >
                ( earthEphemeris, ephemerides, rotationEphemeris,
                  boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                  nominalDisplacementLoveNumbers ,firstStepLatitudeDependenceTerms,
                  areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers ,
                  tudat::paths::getEarthDeformationDataFilesPath( ) + "/diurnalDisplacementFrequencyDependence2.txt",
                  tudat::paths::getEarthDeformationDataFilesPath( ) + "/longPeriodDisplacementFrequencyDependence.txt",
                  std::bind( &calculateFundamentalArgumentsIersCode, std::placeholders::_1 ) );

        evaluationTime = 292852800.0;
        expectedDeformation <<0.7700420357108125891E-01,0.6304056321824967613E-01,0.5516568152597246810E-01;
        calculatedDeformation =  deformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-5 );
    }

    // Test full IERS chapter 7 model (case 2)
    {
        for( int test = 0; test < 2; test++ )
        {
        siteState( 0 ) = 1112189.660;
        siteState( 1 ) = -4842955.026;
        siteState( 2 ) = 3985352.284;
        sunState( 0 ) = -54537460436.2357;
        sunState( 1 ) = 130244288385.279;
        sunState( 2 ) = 56463429031.5996;
        moonState( 0 ) = 300396716.912;
        moonState( 1 ) = 243238281.451;
        moonState( 2 ) = 120548075.939;

        earthEphemeris = std::bind( &Ephemeris::getCartesianState,
                                    std::make_shared< ConstantEphemeris >( earthState ), std::placeholders::_1 );
        sunEphemeris = std::bind( &Ephemeris::getCartesianState,
                                  std::make_shared< ConstantEphemeris >( sunState ), std::placeholders::_1 );
        moonEphemeris = std::bind( &Ephemeris::getCartesianState,
                                   std::make_shared< ConstantEphemeris >( moonState ), std::placeholders::_1 );

        ephemerides[ 0 ] = moonEphemeris;
        ephemerides[ 1 ] = sunEphemeris;

        if( test == 0 )
        {
            deformationModel = std::make_shared< Iers2010EarthDeformation >
                    ( earthEphemeris, ephemerides, rotationEphemeris,
                      boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                      nominalDisplacementLoveNumbers, firstStepLatitudeDependenceTerms,
                      areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers ,
                      tudat::paths::getEarthDeformationDataFilesPath( ) + "/diurnalDisplacementFrequencyDependence2.txt",
                      tudat::paths::getEarthDeformationDataFilesPath( ) + "/longPeriodDisplacementFrequencyDependence.txt",
                      std::bind( &calculateFundamentalArgumentsIersCode, std::placeholders::_1 ) );
        }
        else
        {
            deformationModel = createDefaultEarthIers2010DeformationModel(
                        std::make_shared< ConstantEphemeris >( earthState ),
                        std::make_shared< ConstantEphemeris >( moonState ),
                        std::make_shared< ConstantEphemeris >( sunState ),
                        std::make_shared< ConstantRotationalEphemeris >( orientationQuaternion ),
                        boost::lambda::constant( 1.0 ),
                        lunarMassFunction,
                        solarMassFunction );
        }
        evaluationTime = 395409600.0;
        expectedDeformation <<-0.2036831479592075833E-01,0.5658254776225972449E-01,-0.7597679676871742227E-01;
        calculatedDeformation =  deformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-5 );
        }
    }

    // Check equivalance of basic and IERS model
    {
        for( int i = 0; i < 6; i++ )
        {
            areFirstStepCorrectionsCalculated[ i ] = 0;
        }
        firstStepLatitudeDependenceTerms[ 0 ] = 0.0;
        firstStepLatitudeDependenceTerms[ 1 ] = 0.0;

        std::shared_ptr< BasicTidalBodyDeformation > basicDeformationModel =
                std::make_shared< BasicTidalBodyDeformation >(
                    earthEphemeris, ephemerides, rotationEphemeris,
                                      boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                                      nominalDisplacementLoveNumbers );

        deformationModel = std::make_shared< Iers2010EarthDeformation >
                ( earthEphemeris, ephemerides, rotationEphemeris,
                  boost::lambda::constant( 1.0 ), massFunctions, earthEquatorialRadius,
                  nominalDisplacementLoveNumbers, firstStepLatitudeDependenceTerms,
                  areFirstStepCorrectionsCalculated, correctionLoveAndShidaNumbers, "", "" );

        expectedDeformation = deformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );
        calculatedDeformation = basicDeformationModel->calculateDisplacement( evaluationTime, siteState.segment( 0, 3 ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedDeformation, expectedDeformation, 1.0E-12 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}


