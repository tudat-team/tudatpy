/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Easy calculation. Newton's Law of Gravity Tutorial,
 *          http://easycalculation.com/physics/classical-physics/learn-newtons-law.php, last
 *          accessed: 12th February, 2012.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_gravitational_torque )

//! Test to check degree two gravitational torque
BOOST_AUTO_TEST_CASE( testDegreeTwoGravitationalTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    //! Test torque for various inertia tensors
    for( unsigned int testCase = 0; testCase < 7; testCase++ )
    {


        // Create body objects.
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );

        // Define degree two coefficients
        Eigen::Matrix3d cosineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d sineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d inertiaTensorDeviation = Eigen::Matrix3d::Zero( );
        if( testCase > 0 )
        {
            cosineCoefficients( 0, 0 ) = 1.0;

            if( testCase == 2 )
            {
                cosineCoefficients( 2, 0 ) = 0.01;
            }
            else if( testCase == 3 )
            {
                cosineCoefficients( 2, 2 ) = 0.01;
            }
            else if( testCase == 4 )
            {
                sineCoefficients( 2, 2 ) = 0.01;
            }
            else if( testCase == 5 )
            {
                cosineCoefficients( 2, 1 ) = 0.01;
            }
            else if( testCase == 6 )
            {
                sineCoefficients( 2, 1 ) = 0.01;
            }

            bodySettings[ "Moon" ]->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                        spice_interface::getBodyGravitationalParameter( "Moon" ), spice_interface::getAverageRadius( "Moon" ),
                        cosineCoefficients, sineCoefficients, "IAU_Moon" );
        }

        if( testCase > 0 )
        {
            double c20InertiaContribution  = cosineCoefficients( 2, 0 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) / 3.0;
            inertiaTensorDeviation( 0, 0 ) += c20InertiaContribution;
            inertiaTensorDeviation( 1, 1 ) += c20InertiaContribution;
            inertiaTensorDeviation( 2, 2 ) -= 2.0 * c20InertiaContribution;

            double c21InertiaContribution  = cosineCoefficients( 2, 1 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            inertiaTensorDeviation( 0, 2 ) -= c21InertiaContribution;
            inertiaTensorDeviation( 2, 0 ) -= c21InertiaContribution;

            double c22InertiaContribution  = 2.0 * cosineCoefficients( 2, 2 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
            inertiaTensorDeviation( 0, 0 ) -= c22InertiaContribution;
            inertiaTensorDeviation( 1, 1 ) += c22InertiaContribution;

            double s21InertiaContribution  = sineCoefficients( 2, 1 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            inertiaTensorDeviation( 1, 2 ) -= s21InertiaContribution;
            inertiaTensorDeviation( 2, 1 ) -= s21InertiaContribution;

            double s22InertiaContribution  = 2.0 * sineCoefficients( 2, 2 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
            inertiaTensorDeviation( 0, 1 ) -= s22InertiaContribution;
            inertiaTensorDeviation( 1, 0 ) -= s22InertiaContribution;

            inertiaTensorDeviation *= spice_interface::getAverageRadius( "Moon" ) * spice_interface::getAverageRadius( "Moon" ) *
                    spice_interface::getBodyGravitationalParameter( "Moon" ) / physical_constants::GRAVITATIONAL_CONSTANT;
        }

        // Create bodies
        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Crate torque model
        SelectedTorqueMap selectedTorqueModelMap;
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back(
                    std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                    bodyMap, selectedTorqueModelMap, { "Moon" } );
        std::shared_ptr< TorqueModel > secondDegreeGravitationalTorque =
                torqueModelMap.at( "Moon" ).at( "Earth" ).at( 0 );

        // Update environment to current time
        double evaluationTime = tudat::physical_constants::JULIAN_DAY / 2.0;

        bodyMap.at( "Moon" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setBodyInertiaTensorFromGravityField( 0.4 );

        bodyMap.at( "Earth" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setBodyInertiaTensorFromGravityField( 0.4 );

        if( testCase == 0 )
        {
            // Test id gravity field coefficients are correctly reconstructed
            inertiaTensorDeviation = bodyMap.at( "Moon" )->getBodyInertiaTensor( );
            std::shared_ptr< SphericalHarmonicsGravityField > moonGravityField =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                        bodyMap.at( "Moon" )->getGravityFieldModel( ) );


            double recomputedC20 = std::sqrt( 0.2 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) * (
                        0.5 * ( inertiaTensorDeviation( 0, 0 ) + inertiaTensorDeviation( 1, 1 ) ) -
                        inertiaTensorDeviation( 2, 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedC20, moonGravityField->getCosineCoefficients( )( 2, 0 ), 1.0E-12 );

            double recomputedC21 = -std::sqrt( 3.0 / 5.0 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) * (
                        inertiaTensorDeviation( 0, 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedC21, moonGravityField->getCosineCoefficients( )( 2, 1 ), 1.0E-12 );

            double recomputedC22 = std::sqrt( 3.0 / 20.0 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) * (
                        inertiaTensorDeviation( 1, 1 ) - inertiaTensorDeviation( 0, 0 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedC22, moonGravityField->getCosineCoefficients( )( 2, 2 ), 1.0E-12 );

            double recomputedS21 = -std::sqrt( 3.0 / 5.0 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) * (
                        inertiaTensorDeviation( 1, 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedS21, moonGravityField->getSineCoefficients( )( 2, 1 ), 1.0E-12 );

            double recomputedS22 = -std::sqrt( 12.0 / 5.0 ) /
                    ( 2.0 * moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) * (
                        inertiaTensorDeviation( 0, 1 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedS22, moonGravityField->getSineCoefficients( )( 2, 2 ), 1.0E-12 );
        }

        // Update torque model to current time.
        secondDegreeGravitationalTorque->updateMembers( evaluationTime );

        // Compute torque manually, and test against computed result
        Eigen::Vector3d earthRelativePosition =
                bodyMap.at( "Moon" )->getCurrentRotationToLocalFrame( ) *
                ( bodyMap.at( "Earth" )->getPosition( ) - bodyMap.at( "Moon" )->getPosition( ) );
        Eigen::Vector3d manualTorque = 3.0 * earthRelativePosition.cross( inertiaTensorDeviation * earthRelativePosition ) *
                bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) /
                std::pow( earthRelativePosition.norm( ), 5.0 );

        Eigen::Vector3d currentTorque = secondDegreeGravitationalTorque->getTorque( );
        Eigen::Vector3d torqueError = ( currentTorque - manualTorque );

        std::cout<<"Current torques "<<currentTorque.transpose( )<<std::endl<<
                   torqueError<<std::endl;

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( torqueError( i ) ), 1.0E-14 * currentTorque.norm( ) );
        }
    }
}

//! Test to check spherical harmonic torque.
BOOST_AUTO_TEST_CASE( testSphericalGravitationalTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    //! Test for zero and non-zero spherical harmonic coefficients
    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );

        // Create body objects.
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );

        // Set effective point-mass gravity
        if( testCase == 0 )
        {
            Eigen::Matrix3d cosineCoefficients = Eigen::Matrix3d::Zero( );
            Eigen::Matrix3d sineCoefficients = Eigen::Matrix3d::Zero( );
            cosineCoefficients( 0, 0 ) = 1.0;

            bodySettings[ "Moon" ]->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                        spice_interface::getBodyGravitationalParameter( "Moon" ), spice_interface::getAverageRadius( "Moon" ),
                        cosineCoefficients, sineCoefficients, "IAU_Moon" );
        }

        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Create two torque models
        SelectedTorqueMap selectedTorqueModelMap;
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back(
                    std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back(
                    std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ) );
        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                    bodyMap, selectedTorqueModelMap, { "Moon" } );
        std::shared_ptr< TorqueModel > secondDegreeGravitationalTorque =
                torqueModelMap.at( "Moon" ).at( "Earth" ).at( 0 );
        std::shared_ptr< TorqueModel > sphercialHarmonicGravitationalTorque =
                torqueModelMap.at( "Moon" ).at( "Earth" ).at( 1 );

        // Update Moon to current time.
        double evaluationTime = tudat::physical_constants::JULIAN_DAY / 2.0;
        bodyMap.at( "Moon" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setBodyInertiaTensorFromGravityField( 0.0 );

        {
            // Test reconstructed spherical harmonic coefficients
            Eigen::Matrix3d moonCosineCoefficients =
                    std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >(
                        bodyMap.at( "Moon" )->getGravityFieldModel( ) )->getCosineCoefficients( ).block( 0, 0, 3, 3 );
            Eigen::Matrix3d moonSineCoefficients =
                    std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >(
                        bodyMap.at( "Moon" )->getGravityFieldModel( ) )->getSineCoefficients( ).block( 0, 0, 3, 3 );
            Eigen::MatrixXd moonReconstructedCosineCoefficients = Eigen::Matrix3d::Zero( ),
                    moonReconstructedSineCoefficients = Eigen::Matrix3d::Zero( );

            double reconstructedScaledMeanMomentOfInertia;
            gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                        bodyMap.at( "Moon" )->getBodyInertiaTensor( ),
                        bodyMap.at( "Moon" )->getGravityFieldModel( )->getGravitationalParameter( ),
                        std::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                            bodyMap.at( "Moon" )->getGravityFieldModel( ) )->getReferenceRadius( ), true,
                        moonReconstructedCosineCoefficients, moonReconstructedSineCoefficients,
                        reconstructedScaledMeanMomentOfInertia );

            for( unsigned int j = 0; j < 3; j++ )
            {
                for( unsigned int k = 0; k <= j; k++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( moonReconstructedCosineCoefficients( j, k ) -
                                                  moonCosineCoefficients( j, k ) ), 1.0E-20 );

                    BOOST_CHECK_SMALL( std::fabs( moonReconstructedSineCoefficients( j, k ) -
                                                  moonSineCoefficients( j, k ) ), 1.0E-20 );
                }
            }
        }

        // Update Earth to current time.
        bodyMap.at( "Earth" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setBodyInertiaTensorFromGravityField( 0.0 );

        // Update and compute torque values
        secondDegreeGravitationalTorque->updateMembers( evaluationTime );
        Eigen::Vector3d currentExplicitTorque = secondDegreeGravitationalTorque->getTorque( );
        sphercialHarmonicGravitationalTorque->updateMembers( evaluationTime );
        Eigen::Vector3d currentSphericalHarmonicTorque = sphercialHarmonicGravitationalTorque->getTorque( );

        // Test difference between models (should be equivalent)
        Eigen::Vector3d torqueError = currentExplicitTorque - currentSphericalHarmonicTorque;
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( torqueError( i ) ), 1.0E-14 * currentExplicitTorque.norm( ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
