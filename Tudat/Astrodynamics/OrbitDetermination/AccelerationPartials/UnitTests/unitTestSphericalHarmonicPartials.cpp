/*    Copyright (c) 2010-2012 Delft University of Technology./
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
#include <string>
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/gravitationalParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/numericalAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicPartialFunctions.h"
#include "Tudat/SimulationSetup/createAccelerationPartials.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/defaultBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination::partial_derivatives;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electro_magnetism;

BOOST_AUTO_TEST_SUITE( test_spherical_harmonic_partials )

BOOST_AUTO_TEST_CASE( testSphericalHarmonicPartials )
{
    // Short-cuts.
    using namespace gravitation;


    const double gravitationalParameter = 3.986004418e14;
    const double planetaryRadius = 6378137.0;

    const Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );
    const Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );
    Eigen::Vector3d nominalSphericalPosition = coordinate_conversions::
            convertCartesianToSpherical( position );
    nominalSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 -
            nominalSphericalPosition( 1 );

    boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache
            = boost::make_shared< basic_mathematics::SphericalHarmonicsCache >( 6, 6 );
    sphericalHarmonicsCache->update(
                nominalSphericalPosition( 0 ), std::sin( nominalSphericalPosition( 1 ) ), nominalSphericalPosition( 2 ),
                planetaryRadius );
    boost::shared_ptr< basic_mathematics::LegendreCache > legendreCache = sphericalHarmonicsCache->getLegendreCache( );

    double currentLongitude = sphericalHarmonicsCache->getCurrentLongitude( );
    double currentPolynomialArgument = legendreCache->getCurrentPolynomialParameter( );

    Eigen::MatrixXd upPerturbedLegendrePolynomials = Eigen::MatrixXd( 6, 6 );
    Eigen::MatrixXd downPerturbedLegendrePolynomials = Eigen::MatrixXd( 6, 6 );

    Eigen::MatrixXd upPerturbedLegendrePolynomialPartials = Eigen::MatrixXd( 6, 6 );
    Eigen::MatrixXd downPerturbedLegendrePolynomialPartials = Eigen::MatrixXd( 6, 6 );

    Eigen::MatrixXd analyticalLegendrePolynomialPartials = Eigen::MatrixXd( 6, 6 );
    Eigen::MatrixXd analyticalLegendrePolynomialSecondPartials = Eigen::MatrixXd( 6, 6 );

    legendreCache->setComputeSecondDerivatives( 1 );
    legendreCache->update( currentPolynomialArgument + 0.1  );

    legendreCache->update( currentPolynomialArgument );

    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            analyticalLegendrePolynomialPartials( i, j ) = legendreCache->getLegendrePolynomialDerivative( i, j );
            analyticalLegendrePolynomialSecondPartials( i, j ) = legendreCache->getLegendrePolynomialSecondDerivative( i, j );
        }

    }

    double polynomialArgumentPerturbation = 1.0E-6;
    {
        legendreCache->update( currentPolynomialArgument + polynomialArgumentPerturbation );
        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                upPerturbedLegendrePolynomials( i, j ) = legendreCache->getLegendrePolynomial( i, j );
                upPerturbedLegendrePolynomialPartials( i, j ) = legendreCache->getLegendrePolynomialDerivative( i, j );
            }
        }

        legendreCache->update( currentPolynomialArgument - polynomialArgumentPerturbation );
        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                downPerturbedLegendrePolynomials( i, j ) = legendreCache->getLegendrePolynomial( i, j );
                downPerturbedLegendrePolynomialPartials( i, j ) = legendreCache->getLegendrePolynomialDerivative( i, j );
            }
        }
    }

    Eigen::MatrixXd numericalLegendrePolynomialPartials =
            ( upPerturbedLegendrePolynomials - downPerturbedLegendrePolynomials ) /
            ( 2.0 * polynomialArgumentPerturbation );
    Eigen::MatrixXd numericalLegendrePolynomialSecondPartials =
            ( upPerturbedLegendrePolynomialPartials - downPerturbedLegendrePolynomialPartials ) /
            ( 2.0 * polynomialArgumentPerturbation );

    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( numericalLegendrePolynomialPartials( i, j ) - analyticalLegendrePolynomialPartials( i, j ) ), 1.0E-8 );
            BOOST_CHECK_SMALL(
                        std::fabs( numericalLegendrePolynomialSecondPartials( i, j ) - analyticalLegendrePolynomialSecondPartials( i, j ) ), 1.0E-8 );
        }
    }


    std::vector< std::vector< Eigen::Vector3d > > sphericalPotentialGradients;

    std::vector< std::vector< Eigen::Matrix3d > > upPerturbedSphericalPotentialGradients;
    std::vector< std::vector< Eigen::Matrix3d > > downPerturbedSphericalPotentialGradients;

    std::vector< std::vector< Eigen::Matrix3d > > numericalSphericalPotentialHessian;
    std::vector< std::vector< Eigen::Matrix3d > > analyticalSphericalPotentialHessian;

    Eigen::Matrix3d normalization;
    normalization <<  nominalSphericalPosition( 0 ) * nominalSphericalPosition( 0 ), nominalSphericalPosition( 0 ), nominalSphericalPosition( 0 ),
            nominalSphericalPosition( 0 ), 1.0 , 1.0,
            nominalSphericalPosition( 0 ), 1.0, 1.0;
    for( unsigned int i = 0; i < 6; i++ )
    {
        std::vector< Eigen::Vector3d > singleTermPotentialGradients;
        singleTermPotentialGradients.resize( 6 );
        sphericalPotentialGradients.push_back( singleTermPotentialGradients );

        std::vector< Eigen::Matrix3d > singleTermPotentialPartials;
        singleTermPotentialPartials.resize( 6 );
        upPerturbedSphericalPotentialGradients.push_back( singleTermPotentialPartials );
        downPerturbedSphericalPotentialGradients.push_back( singleTermPotentialPartials );

        numericalSphericalPotentialHessian.push_back( singleTermPotentialPartials );
        analyticalSphericalPotentialHessian.push_back( singleTermPotentialPartials );


    }

    sphericalHarmonicsCache->update( position.norm( ), currentPolynomialArgument, currentLongitude, planetaryRadius );

    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            sphericalPotentialGradients[ i ][ j ] = basic_mathematics::computePotentialGradient(
                        nominalSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                        cosineCoefficients( i, j ), sineCoefficients( i, j ), legendreCache->getLegendrePolynomial( i, j ),
                        legendreCache->getLegendrePolynomialDerivative( i, j ), sphericalHarmonicsCache );
        }
    }


    Eigen::Vector3d perturbedSphericalPosition;
    Eigen::Vector3d sphericalStatePerturbation;
    sphericalStatePerturbation<<10.0, 1.0E-7, 1.0E-8;

    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) += sphericalStatePerturbation( parameter );

        sphericalHarmonicsCache->update(
                    perturbedSphericalPosition( 0 ), std::sin( perturbedSphericalPosition( 1 ) ),
                    perturbedSphericalPosition( 2 ), planetaryRadius );

        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                upPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) = basic_mathematics::computePotentialGradient(
                            perturbedSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                            cosineCoefficients( i, j ), sineCoefficients( i, j ), legendreCache->getLegendrePolynomial( i, j ),
                            legendreCache->getLegendrePolynomialDerivative( i, j ), sphericalHarmonicsCache );
            }
        }

        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) -= sphericalStatePerturbation( parameter );

        sphericalHarmonicsCache->update(
                    perturbedSphericalPosition( 0 ), std::sin( perturbedSphericalPosition( 1 ) ),
                    perturbedSphericalPosition( 2 ), planetaryRadius );

        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                downPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) = basic_mathematics::computePotentialGradient(
                            perturbedSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                            cosineCoefficients( i, j ), sineCoefficients( i, j ), legendreCache->getLegendrePolynomial( i, j ),
                            legendreCache->getLegendrePolynomialDerivative( i, j ), sphericalHarmonicsCache );
            }
        }

        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                numericalSphericalPotentialHessian[ i ][ j ].block( 0, parameter, 3, 1 ) =
                        ( ( upPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) -
                            downPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) ) /
                          ( 2.0 * sphericalStatePerturbation( parameter ) ) );

            }

        }
    }

    for( unsigned int i = 1; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            computePotentialSphericalHessian(
                        nominalSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                        cosineCoefficients( i, j ), sineCoefficients( i, j ), sphericalHarmonicsCache,
                        analyticalSphericalPotentialHessian[ i ][ j ] );

            analyticalSphericalPotentialHessian[ i ][ j ] = analyticalSphericalPotentialHessian[ i ][ j ].cwiseProduct(
                        normalization );
            numericalSphericalPotentialHessian[ i ][ j ] = numericalSphericalPotentialHessian[ i ][ j ].cwiseProduct(
                        normalization );
            for( unsigned int k = 0; k < 3; k++ )
            {
                for( unsigned int l = 0; l < 3; l++ )
                {
                    BOOST_CHECK_SMALL( analyticalSphericalPotentialHessian[ i ][ j ]( k, l ) -
                            numericalSphericalPotentialHessian[ i ][ j ]( k, l ), 1.0E-5 );
                }

            }

        }
    }

    Eigen::Matrix3d cumulativeSphericalHessian =  computeCumulativeSphericalHessian(
                nominalSphericalPosition, planetaryRadius, gravitationalParameter, cosineCoefficients, sineCoefficients,
                sphericalHarmonicsCache );

    Eigen::Matrix3d nominalGradientTransformationMatrix =
            coordinate_conversions::getSphericalToCartesianGradientMatrix( position );

    Eigen::Vector3d upPerturbedTotalGradient;
    Eigen::Vector3d downPerturbedTotalGradient;
    Eigen::Vector3d perturbedCartesianPosition;

    Eigen::Matrix3d numericalTotalSphericalGradient;

    sphericalStatePerturbation( 0 ) *= 10.0;
    sphericalStatePerturbation( 1 ) *= 100.0;
    sphericalStatePerturbation( 2 ) *= 1000.0;

    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) += sphericalStatePerturbation( parameter );

        perturbedSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - perturbedSphericalPosition( 1 );
        perturbedCartesianPosition = coordinate_conversions::convertSphericalToCartesian(
                    perturbedSphericalPosition );

        upPerturbedTotalGradient =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache );

        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) -= sphericalStatePerturbation( parameter );

        perturbedSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - perturbedSphericalPosition( 1 );
        perturbedCartesianPosition = coordinate_conversions::convertSphericalToCartesian(
                    perturbedSphericalPosition );

        downPerturbedTotalGradient =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache );

        numericalTotalSphericalGradient.block( 0, parameter, 3, 1 ) =
                ( upPerturbedTotalGradient - downPerturbedTotalGradient ) / ( 2.0 * sphericalStatePerturbation( parameter ) );

    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                cumulativeSphericalHessian, numericalTotalSphericalGradient, 1.0E-6 );

    Eigen::Vector3d cartesianStatePerturbation;
    cartesianStatePerturbation<<10.0, 10.0, 10.0;

    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) += cartesianStatePerturbation( parameter );

        upPerturbedTotalGradient =
                nominalGradientTransformationMatrix *
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache );

        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) -= cartesianStatePerturbation( parameter );

        downPerturbedTotalGradient =
                nominalGradientTransformationMatrix *
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache );

        numericalTotalSphericalGradient.block( 0, parameter, 3, 1 ) =
                ( upPerturbedTotalGradient - downPerturbedTotalGradient ) / ( 2.0 * cartesianStatePerturbation( parameter ) );

    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( nominalGradientTransformationMatrix * cumulativeSphericalHessian * nominalGradientTransformationMatrix.transpose( ) ),
                numericalTotalSphericalGradient, 1.0E-6 );


    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) += cartesianStatePerturbation( parameter );

        upPerturbedTotalGradient =
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache );

        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) -= cartesianStatePerturbation( parameter );

        downPerturbedTotalGradient =
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache );

        numericalTotalSphericalGradient.block( 0, parameter, 3, 1 ) =
                ( upPerturbedTotalGradient - downPerturbedTotalGradient ) / ( 2.0 * cartesianStatePerturbation( parameter ) );

    }


    Eigen::Matrix3d totalGradientCartesianPartial =
            computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
                position, planetaryRadius, gravitationalParameter, cosineCoefficients, sineCoefficients,
                sphericalHarmonicsCache );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                totalGradientCartesianPartial, numericalTotalSphericalGradient, 1.0E-6 );
}

BOOST_AUTO_TEST_CASE( testSphericalHarmonicAccelerationpartial )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    // Create empty bodies, earth and vehicle.
    boost::shared_ptr< Body > earth = boost::make_shared< Body >( );
    boost::shared_ptr< Body > vehicle = boost::make_shared< Body >( );

    const double gravitationalParameter = 3.986004418e14;
    const double planetaryRadius = 6378137.0;

    Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );


    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = earth;
    bodyMap[ "Vehicle" ] = vehicle;

    boost::shared_ptr< gravitation::SphericalHarmonicsGravityField > earthGravityField =
            boost::make_shared< gravitation::SphericalHarmonicsGravityField  >(
                gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients, "IAU_Earth" );

    earth->setGravityFieldModel( earthGravityField );

    boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            boost::make_shared< ephemerides::SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000" , "IAU_Earth", 0.0 ),
                2.0 * mathematical_constants::PI / 86400.0,
                1.0E7, basic_astrodynamics::JULIAN_DAY_ON_J2000,
                "ECLIPJ2000" , "IAU_Earth" );

    earth->setRotationalEphemeris( simpleRotationalEphemeris );

    // Set current state of vehicle and earth.
    earth->setState( basic_mathematics::Vector6d::Zero( ) );

    // Set Keplerian elements for Asterix.
    basic_mathematics::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    basic_mathematics::Vector6d asterixInitialState = orbital_element_conversions::convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, gravitationalParameter );


    vehicle->setState( asterixInitialState );



    // Create acceleration due to vehicle on earth.
    boost::shared_ptr< SphericalHarmonicAccelerationSettings > accelerationSettings =
            boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 );
    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > gravitationalAcceleration =
            boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( vehicle, earth, accelerationSettings, "Vehicle", "Earth" ) );

    gravitationalAcceleration->updateMembers( 0.0 );
    Eigen::Vector3d nominalAcceleration = gravitationalAcceleration->getAcceleration( );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation<<10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation<< 1.0E-3, 1.0E-3, 1.0E-3;

    // Create state access/modification functions for bodies.
    boost::function< void( basic_mathematics::Vector6d ) > earthStateSetFunction =
            boost::bind( &Body::setState, earth, _1  );
    boost::function< void( basic_mathematics::Vector6d ) > vehicleStateSetFunction =
            boost::bind( &Body::setState, vehicle, _1  );
    boost::function< basic_mathematics::Vector6d ( ) > earthStateGetFunction =
            boost::bind( &Body::getState, earth );
    boost::function< basic_mathematics::Vector6d ( ) > vehicleStateGetFunction =
            boost::bind( &Body::getState, vehicle );



    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", constant_rotation_rate ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", rotation_pole_position ) );

    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 5, 4, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 1, 5, 4, "Earth", spherical_harmonics_sine_coefficient_block ) );

    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodyMap );

    // Create acceleration partial object.
    boost::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartial =
            boost::dynamic_pointer_cast< SphericalHarmonicsGravityPartial > (
                createAnalyticalAccelerationPartial(
                    gravitationalAcceleration, std::make_pair( "Vehicle", vehicle ), std::make_pair( "Earth", earth ),
                    bodyMap, parameterSet ) );

    double testTime = 1.0E6;
    earth->setCurrentRotationToLocalFrameFromEphemeris( testTime );

    accelerationPartial->update( testTime );

    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );





    // Calculate numerical partials.
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), positionPerturbation, 0 );
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), velocityPerturbation, 3 );

    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3 );



    // Calculate numerical partials.
    std::map< int, boost::shared_ptr< EstimatableParameter< double > > > doubleParameters =
            parameterSet->getDoubleParameters( );
    std::map< int, boost::shared_ptr< EstimatableParameter< double > > >::iterator doubleParametersIterator =
            doubleParameters.begin( );

    Eigen::Vector3d testPartialWrtEarthGravitationalParameter = calculateAccelerationWrtParameterPartials(
                doubleParametersIterator->second, gravitationalAcceleration, 1.0E12 );
    Eigen::Vector3d partialWrtEarthGravitationalParameter = accelerationPartial->wrtParameter(
                doubleParametersIterator->second );
    doubleParametersIterator++;

    Eigen::Vector3d partialWrtEarthRotationRate = accelerationPartial->wrtParameter(
                doubleParametersIterator->second );
    Eigen::Vector3d testPartialWrtEarthRotationRate = calculateAccelerationWrtParameterPartials(
                doubleParametersIterator->second, gravitationalAcceleration, 1.0E-12, &emptyFunction, testTime, boost::bind(
                    &Body::setCurrentRotationToLocalFrameFromEphemeris, earth, _1 ) );




    std::map< int, boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parameterSet->getVectorParameters( );
    std::map< int, boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >::iterator vectorParametersIterator =
            vectorParameters.begin( );
    Eigen::MatrixXd partialWrtPolePosition = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtPosition = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, Eigen::Vector2d::Constant( 1.0E-6 ),
                &emptyFunction, testTime, boost::bind(
                    &Body::setCurrentRotationToLocalFrameFromEphemeris, earth, _1 ) );
    vectorParametersIterator++;

    Eigen::MatrixXd partialWrtCosineCoefficients = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtCosineCoefficients = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, vectorParametersIterator->second->getParameterValue( ) * 1.0E-2 );
    vectorParametersIterator++;

    Eigen::MatrixXd partialWrtSineCoefficients = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtSineCoefficients = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, vectorParametersIterator->second->getParameterValue( ) * 1.0E-2 );



    //    Eigen::Vector3d testPartialWrtPolePosition = calculateAccelerationWrtParameterPartials(
    //                earthPolePosition, accelerationModel, ( Eigen::VectorXd( 2 ) << 1.0E-4, 1.0E-4 ).finished( ), &emptyFunction, 1.0E6, boost::bind(
    //                    &CelestialBody::setCurrentRotationToLocalFrameFromEphemeris, earth, _1 ) );
    //    Eigen::MatrixXd testPartialWrtCosineCoefficients = calculateAccelerationWrtParameterPartials(
    //                cosineCoefficients, accelerationModel, cosineCoefficients->getParameterValue( ) * 1.0E-2 );
    //    Eigen::MatrixXd testPartialWrtSineCoefficients = calculateAccelerationWrtParameterPartials(
    //                sineCoefficients, accelerationModel, sineCoefficients->getParameterValue( ) * 1.0E-2 );

    //    Eigen::Vector3d partialWrtEarthPolePosition = accelerationPartial->wrtParameter(
    //                earthPolePosition );
    //    Eigen::MatrixXd partialWrtCosineCoefficients = accelerationPartial->wrtParameter(
    //                cosineCoefficients );
    //    Eigen::MatrixXd partialWrtSineCoefficients = accelerationPartial->wrtParameter(
    //                sineCoefficients );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition, partialWrtVehiclePosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity, partialWrtVehicleVelocity, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition, partialWrtEarthPosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity, partialWrtEarthVelocity, 1.0E-6 );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthGravitationalParameter, partialWrtEarthGravitationalParameter, 1.0E-12 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthRotationRate, partialWrtEarthRotationRate, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPosition, partialWrtPolePosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtCosineCoefficients, partialWrtCosineCoefficients, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSineCoefficients, partialWrtSineCoefficients, 1.0E-6 );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





