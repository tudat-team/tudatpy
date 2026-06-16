/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <memory>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/io/basicInputOutput.h"

#include "tudat/math/basic/sphericalHarmonicTransformations.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/simulation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_mathematics;

BOOST_AUTO_TEST_SUITE( test_spherical_harmonic_transformations )

BOOST_AUTO_TEST_CASE( testAnalyticalTransformations )
{
    Eigen::MatrixXd nominalCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd nominalSineCoefficients = Eigen::Matrix3d::Zero( );

    Eigen::MatrixXd nominalNormalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd nominalNormalizedSineCoefficients = Eigen::Matrix3d::Zero( );
    nominalCosineCoefficients( 0, 0 ) = 1.0;

    Eigen::MatrixXd transformedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd transformedSineCoefficients = Eigen::Matrix3d::Zero( );

    Eigen::MatrixXd transformedNormalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd transformedNormalizedSineCoefficients = Eigen::Matrix3d::Zero( );

    Eigen::MatrixXd transformedRenormalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd transformedRenormalizedSineCoefficients = Eigen::Matrix3d::Zero( );

    SphericalHarmonicTransformationCache sphericalHarmonicTransformationCache( 2, 2 );

    double anglePsi = 0.0;
    double anglePhi = 0.0;
    double angleTheta = mathematical_constants::PI / 2.0;

    sphericalHarmonicTransformationCache.updateFrom313EulerAngles( anglePhi, angleTheta, anglePsi );

    double perturbationMagnitude = 0.1;

    {
        Eigen::MatrixXd normalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::MatrixXd normalizedSineCoefficients = Eigen::Matrix3d::Zero( );
        nominalCosineCoefficients( 1, 0 ) = perturbationMagnitude;

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                nominalCosineCoefficients, nominalSineCoefficients, transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients( nominalCosineCoefficients,
                                                                               nominalSineCoefficients,
                                                                               nominalNormalizedCosineCoefficients,
                                                                               nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree( nominalNormalizedCosineCoefficients,
                                                                            nominalNormalizedSineCoefficients,
                                                                            transformedNormalizedCosineCoefficients,
                                                                            transformedNormalizedSineCoefficients,
                                                                            1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients( transformedNormalizedCosineCoefficients,
                                                                               transformedNormalizedSineCoefficients,
                                                                               transformedRenormalizedCosineCoefficients,
                                                                               transformedRenormalizedSineCoefficients );

        nominalCosineCoefficients( 1, 0 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 1 && j == 1 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) - perturbationMagnitude ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) - perturbationMagnitude ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
            }
        }

        //        std::cout<<"Transformed: "<<std::endl<<transformedNormalizedCosineCoefficients<<std::endl<<std::endl<<
        //                   transformedNormalizedSineCoefficients<<std::endl<<std::endl;
        //        std::cout<<"Original: "<<std::endl<<nominalNormalizedCosineCoefficients<<std::endl<<std::endl<<
        //                   nominalNormalizedSineCoefficients<<std::endl<<std::endl<<std::endl;
    }

    {
        nominalCosineCoefficients( 1, 1 ) = perturbationMagnitude;
        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                nominalCosineCoefficients, nominalSineCoefficients, transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients( nominalCosineCoefficients,
                                                                               nominalSineCoefficients,
                                                                               nominalNormalizedCosineCoefficients,
                                                                               nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree( nominalNormalizedCosineCoefficients,
                                                                            nominalNormalizedSineCoefficients,
                                                                            transformedNormalizedCosineCoefficients,
                                                                            transformedNormalizedSineCoefficients,
                                                                            1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients( transformedNormalizedCosineCoefficients,
                                                                               transformedNormalizedSineCoefficients,
                                                                               transformedRenormalizedCosineCoefficients,
                                                                               transformedRenormalizedSineCoefficients );

        nominalCosineCoefficients( 1, 1 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 1 && j == 1 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - perturbationMagnitude ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - perturbationMagnitude ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
            }
        }
    }

    {
        nominalSineCoefficients( 1, 1 ) = perturbationMagnitude;
        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                nominalCosineCoefficients, nominalSineCoefficients, transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients( nominalCosineCoefficients,
                                                                               nominalSineCoefficients,
                                                                               nominalNormalizedCosineCoefficients,
                                                                               nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree( nominalNormalizedCosineCoefficients,
                                                                            nominalNormalizedSineCoefficients,
                                                                            transformedNormalizedCosineCoefficients,
                                                                            transformedNormalizedSineCoefficients,
                                                                            1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients( transformedNormalizedCosineCoefficients,
                                                                               transformedNormalizedSineCoefficients,
                                                                               transformedRenormalizedCosineCoefficients,
                                                                               transformedRenormalizedSineCoefficients );

        //        std::cout<<"Transformed: "<<std::endl<<transformedNormalizedCosineCoefficients<<std::endl<<std::endl<<
        //                   transformedNormalizedSineCoefficients<<std::endl<<std::endl;
        //        std::cout<<"Original: "<<std::endl<<nominalNormalizedCosineCoefficients<<std::endl<<std::endl<<
        //                   nominalNormalizedSineCoefficients<<std::endl<<std::endl<<std::endl;

        //        std::cout<<"Transformed: "<<std::endl<<transformedCosineCoefficients<<std::endl<<std::endl<<
        //                   transformedSineCoefficients<<std::endl<<std::endl;
        //        std::cout<<"Original: "<<std::endl<<nominalCosineCoefficients<<std::endl<<std::endl<<
        //                   nominalSineCoefficients<<std::endl<<std::endl<<std::endl;

        nominalSineCoefficients( 1, 1 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 1 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) + perturbationMagnitude ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) + perturbationMagnitude ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
            }
        }
    }

    {
        nominalCosineCoefficients( 2, 0 ) = perturbationMagnitude;

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                nominalCosineCoefficients, nominalSineCoefficients, transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients( nominalCosineCoefficients,
                                                                               nominalSineCoefficients,
                                                                               nominalNormalizedCosineCoefficients,
                                                                               nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree( nominalNormalizedCosineCoefficients,
                                                                            nominalNormalizedSineCoefficients,
                                                                            transformedNormalizedCosineCoefficients,
                                                                            transformedNormalizedSineCoefficients,
                                                                            1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients( transformedNormalizedCosineCoefficients,
                                                                               transformedNormalizedSineCoefficients,
                                                                               transformedRenormalizedCosineCoefficients,
                                                                               transformedRenormalizedSineCoefficients );

        nominalSineCoefficients( 2, 0 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 2 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) + perturbationMagnitude / 2.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) + perturbationMagnitude / 2.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }

                else if( i == 2 && j == 2 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) + perturbationMagnitude / 4.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) + perturbationMagnitude / 4.0 ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ),
                                       std::numeric_limits< double >::epsilon( ) );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testSphericalHarmonicTransformationFromAcceleration )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings( "Earth", TUDAT_NAN, TUDAT_NAN );
    bodySettings[ "Earth2" ] = getDefaultSingleBodySettings( "Earth", TUDAT_NAN, TUDAT_NAN );

    Eigen::Quaterniond rotationToEarthFixed = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    bodySettings[ "Earth" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000", "IAU_Earth", rotationToEarthFixed, 0.0, 2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    Eigen::Quaterniond rotationToEarth2Fixed =
            Eigen::Quaterniond( Eigen::AngleAxisd( 0.1, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitX( ) ) *
                                Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitZ( ) ) );
    bodySettings[ "Earth2" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000", "IAU_Mars", rotationToEarth2Fixed, 0.0, 1.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    Eigen::MatrixXd nominalCosineCoefficientsFull =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings[ "Earth2" ]->gravityFieldSettings )
                    ->getCosineCoefficients( );
    Eigen::MatrixXd nominalSineCoefficientsFull =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings[ "Earth2" ]->gravityFieldSettings )
                    ->getSineCoefficients( );

    Eigen::MatrixXd nominalCosineCoefficients = nominalCosineCoefficientsFull.block( 0, 0, 7, 7 );
    Eigen::MatrixXd nominalSineCoefficients = nominalSineCoefficientsFull.block( 0, 0, 7, 7 );

    //    nominalCosineCoefficients( 2, 0 ) = nominalCosineCoefficientsFull( 2, 0 );
    //    nominalCosineCoefficients( 2, 1 ) = nominalCosineCoefficientsFull( 2, 1 );
    //    nominalCosineCoefficients( 2, 2 ) = nominalCosineCoefficientsFull( 2, 2 );

    //    nominalSineCoefficients( 2, 1 ) = nominalSineCoefficientsFull( 2, 1 );
    //    nominalSineCoefficients( 2, 2 ) = nominalSineCoefficientsFull( 2, 2 );

    bodySettings[ "Earth" ]->gravityFieldSettings =
            std::make_shared< SphericalHarmonicsGravityFieldSettings >( tudat::spice_interface::getBodyGravitationalParameter( "Earth" ),
                                                                        tudat::spice_interface::getAverageRadius( "Earth" ),
                                                                        nominalCosineCoefficients,
                                                                        nominalSineCoefficients,
                                                                        "IAU_Earth" );
    bodySettings[ "Earth2" ]->gravityFieldSettings =
            std::make_shared< SphericalHarmonicsGravityFieldSettings >( tudat::spice_interface::getBodyGravitationalParameter( "Earth" ),
                                                                        tudat::spice_interface::getAverageRadius( "Earth" ),
                                                                        nominalCosineCoefficients,
                                                                        nominalSineCoefficients,
                                                                        "IAU_Mars" );

    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
    bodySettings[ "Earth2" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );

    SystemOfBodies bodyMap = createSystemOfBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap.addBody( std::make_shared< simulation_setup::Body >( ), "Asterix" );

    bodyMap.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );
    bodyMap.at( "Earth2" )->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );

    Eigen::Quaterniond rotationToEarthFixedFrame = bodyMap.at( "Earth" )->getCurrentRotationToLocalFrame( );
    Eigen::Quaterniond rotationToEarth2FixedFrame = bodyMap.at( "Earth2" )->getCurrentRotationToLocalFrame( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap, accelerationModelMap2;
    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );
        accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
    }

    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth2" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth2" );
        accelerationModelMap2 = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
    }

    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationFromEarth =
            accelerationModelMap.at( "Asterix" ).at( "Earth" ).at( 0 );
    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationFromEarth2 =
            accelerationModelMap2.at( "Asterix" ).at( "Earth2" ).at( 0 );

    SphericalHarmonicTransformationCache sphericalHarmonicTransformationCacheEarth( 6, 6 );
    sphericalHarmonicTransformationCacheEarth.updateFromQuaternion( rotationToEarthFixedFrame );
    SphericalHarmonicTransformationCache sphericalHarmonicTransformationCacheEarth2( 6, 6 );
    sphericalHarmonicTransformationCacheEarth2.updateFromQuaternion( rotationToEarth2FixedFrame );

    Eigen::MatrixXd earthCosineCoefficients, earthSineCoefficients;
    Eigen::MatrixXd earth2CosineCoefficients, earth2SineCoefficients;

    earthCosineCoefficients.setZero( 7, 7 );
    earthCosineCoefficients( 0, 0 ) = 1.0;
    earthSineCoefficients.setZero( 7, 7 );

    earth2CosineCoefficients.setZero( 7, 7 );
    earth2CosineCoefficients( 0, 0 ) = 1.0;
    earth2SineCoefficients.setZero( 7, 7 );

    sphericalHarmonicTransformationCacheEarth.transformCoefficientsAtDegree(
            nominalCosineCoefficients, nominalSineCoefficients, earthCosineCoefficients, earthSineCoefficients, true );
    sphericalHarmonicTransformationCacheEarth2.transformCoefficientsAtDegree(
            nominalCosineCoefficients, nominalSineCoefficients, earth2CosineCoefficients, earth2SineCoefficients, true );

    earthCosineCoefficients( 0, 0 ) = 0.0;
    earth2CosineCoefficients( 0, 0 ) = 0.0;

    std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >( bodyMap.at( "Earth" )->getGravityFieldModel( ) )
            ->setCosineCoefficients( earthCosineCoefficients );
    std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >( bodyMap.at( "Earth" )->getGravityFieldModel( ) )
            ->setSineCoefficients( earthSineCoefficients );

    std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >( bodyMap.at( "Earth2" )->getGravityFieldModel( ) )
            ->setCosineCoefficients( earth2CosineCoefficients );
    std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >( bodyMap.at( "Earth2" )->getGravityFieldModel( ) )
            ->setSineCoefficients( earth2SineCoefficients );

    bodyMap.at( "Earth" )->setStateFromEphemeris( 0.0 );
    bodyMap.at( "Earth2" )->setStateFromEphemeris( 0.0 );

    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    bodyMap.at( "Asterix" )->setState( asterixInitialState );

    accelerationFromEarth->updateMembers( 0.0 );
    accelerationFromEarth2->updateMembers( 0.0 );

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( accelerationFromEarth->getAcceleration( )( i ) - accelerationFromEarth2->getAcceleration( )( i ) ),
                           1.0E-17 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
